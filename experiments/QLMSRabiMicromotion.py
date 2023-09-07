import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.system.subsequences import TickleDDS
import LAX_exp.experiments.SidebandCooling as SidebandCooling

# LABRAD IMPORTS
import labrad
from os import environ
from EGGS_labrad.config.dc_config import dc_config


class QLMSRabiMicromotion(SidebandCooling.SidebandCooling):
    """
    Experiment: QLMS Rabi - Micromotion

    # todo: redocument
    Quantum Logic Mass Spectroscopy - Rabi Flopping
    Cool the ions to the ground state of motion via sideband cooling,
    then apply a tickle to create a coherent state to be read out via RSB/BSB comparison.
    """
    name = 'QLMSRabiMicromotion'


    def build_experiment(self):
        # run regular sideband cooling build
        super().build_experiment()

        # QLMS configuration
        self.setattr_argument("freq_qlms_rabi_khz_list",                        Scannable(
                                                                                    default=[
                                                                                        ExplicitScan([1088]),
                                                                                        CenterScan(1088, 10, 0.5, randomize=True)
                                                                                    ],
                                                                                    global_min=0, global_max=10000, global_step=1,
                                                                                    unit="kHz", scale=1, ndecimals=3
                                                                                ), group=self.name)

        # shimming voltages
        self.dc_micromotion_channeldict =                                       dc_config.channeldict
        self.setattr_argument("dc_micromotion_channel_1",                       EnumerationValue(list(self.dc_micromotion_channeldict.keys()), default='H Shim'), group=self.name)
        self.setattr_argument("dc_micromotion_voltages_1_v_list",               Scannable(
                                                                                    default=[
                                                                                        ExplicitScan([35.5]),
                                                                                        CenterScan(35.5, 10., 1., randomize=True)
                                                                                    ],
                                                                                    global_min=0, global_max=400, global_step=1,
                                                                                    unit="V", scale=1, ndecimals=1
                                                                                ), group=self.name)

        self.setattr_argument("dc_micromotion_channel_2",                       EnumerationValue(list(self.dc_micromotion_channeldict.keys()), default='V Shim'), group=self.name)
        self.setattr_argument("dc_micromotion_voltages_2_v_list",               Scannable(
                                                                                    default=[
                                                                                        ExplicitScan([35.5]),
                                                                                        CenterScan(57., 10., 1., randomize=True)
                                                                                    ],
                                                                                    global_min=0, global_max=400, global_step=1,
                                                                                    unit="V", scale=1, ndecimals=1
                                                                                ), group=self.name)

        # subsequences
        self.tickle_subsequence =                                               TickleDDS(self)

        # get relevant devices
        self.setattr_device('dds_modulation')
        # tmp remove
        self.setattr_device('ttl8')
        # tmp remove

    def prepare_experiment(self):
        # run preparations for sideband cooling
        super().prepare_experiment()

        # get voltage parameters
        self.dc_micromotion_channel_1_num =                                     self.dc_micromotion_channeldict[self.dc_micromotion_channel_1]['num']
        self.dc_micromotion_voltages_1_v_list =                                 np.array(list(self.dc_micromotion_voltages_1_v_list))

        self.dc_micromotion_channel_2_num =                                     self.dc_micromotion_channeldict[self.dc_micromotion_channel_2]['num']
        self.dc_micromotion_voltages_2_v_list =                                 np.array(list(self.dc_micromotion_voltages_2_v_list))

        # convert QLMS modulation to machine units
        self.freq_qlms_rabi_ftw_list =                                          np.array([
                                                                                    self.dds_modulation.frequency_to_ftw(freq_khz * kHz)
                                                                                    for freq_khz in self.freq_qlms_rabi_khz_list
                                                                                ])


        # connect to labrad
        self.cxn =                                                              labrad.connect(environ['LABRADHOST'],
                                                                                               port=7682, tls_mode='off',
                                                                                               username='', password='lab')
        self.dc =                                                               self.cxn.dc_server

        # create an array of values for the experiment to sweep:
        # [readout AOM FTW, DDS tickle frequency, H shim voltage (volts), V shim voltage (volts)]
        self.config_qlms_rabi_micromotion_list =                                np.stack(np.meshgrid(self.freq_qlms_rabi_ftw_list,
                                                                                                     self.freq_readout_ftw_list,
                                                                                                     self.dc_micromotion_voltages_1_v_list,
                                                                                                     self.dc_micromotion_voltages_2_v_list),
                                                                                         -1).reshape(-1, 4)
        np.random.shuffle(self.config_qlms_rabi_micromotion_list)

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_qlms_rabi_micromotion_list),
                5)


    # LABRAD FUNCTIONS
    @rpc
    def voltage_set(self, channel: TInt32, voltage_v: TFloat) -> TNone:
        """
        Set the channel to the desired voltage.
        """
        # set desired voltage
        voltage_set_v = self.dc.voltage_fast(channel, voltage_v)

    @rpc(flags={"async"})
    def prepareDevicesLabrad(self):
        """
        Prepare LabRAD devices for the experiment via RPC.
        """
        # set up amo8
        self.dc.polling(False)
        self.dc.alarm(False)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self):
        self.core.break_realtime()

        # set up labrad devices via RPC
        self.prepareDevicesLabrad()

        # record subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.sidebandcool_subsequence.record_dma()
        self.tickle_subsequence.record_dma()

        # record custom readout sequence
        # note: this is necessary since DMA sequences will preserve urukul attenuation register
        with self.core_dma.record('_SBC_READOUT'):
            # set readout waveform for qubit
            self.qubit.set_profile(0)
            self.qubit.set_att_mu(self.att_readout_mu)

            # transfer population to D-5/2 state
            self.rabiflop_subsequence.run()

            # read out fluorescence
            self.readout_subsequence.run()


    @kernel(flags={"fast-math"})
    def run_main(self):
        self.core.reset()

        # get custom readout handle
        _handle_sbc_readout = self.core_dma.get_handle('_SBC_READOUT')
        self.core.break_realtime()

        for trial_num in range(self.repetitions):

            # sweep experiment configurations
            for config_vals in self.config_qlms_rabi_micromotion_list:

                # extract values from config list
                freq_qlms_ftw =     np.int32(config_vals[0])
                freq_readout_ftw =  np.int32(config_vals[1])
                voltage_v_1 =       config_vals[2]
                voltage_v_2 =       config_vals[3]
                self.core.break_realtime()
                # tmp remove
                self.ttl8.off()
                delay_mu(10000)
                # tmp remove

                # set DC shimming voltage
                # note: we add delays here to ensure that voltages have time to settle
                self.voltage_set(self.dc_micromotion_channel_1_num, voltage_v_1)
                self.core.break_realtime()
                delay_mu(250000)
                self.voltage_set(self.dc_micromotion_channel_2_num, voltage_v_2)
                self.core.break_realtime()
                delay_mu(250000)
                # todo: do we need to add extra wait time for voltages to settle?

                # set QLMS modulation and 729nm readout frequencies
                with parallel:
                    self.dds_modulation.set_mu(freq_qlms_ftw, asf=self.dds_modulation.ampl_modulation_asf, profile=0)
                    self.qubit.set_mu(freq_readout_ftw, asf=self.ampl_readout_pipulse_asf, profile=0)
                    self.core.break_realtime()

                # initialize ion in S-1/2 state
                self.initialize_subsequence.run_dma()

                # sideband cool
                self.sidebandcool_subsequence.run_dma()

                # QLMS tickle
                self.tickle_subsequence.run_dma()

                # custom SBC readout
                self.core_dma.playback_handle(_handle_sbc_readout)

                # update dataset
                with parallel:
                    self.update_results(freq_readout_ftw,
                                        self.readout_subsequence.fetch_count(),
                                        freq_qlms_ftw,
                                        voltage_v_1,
                                        voltage_v_2)
                    self.core.break_realtime()

            # rescue ion as needed
            self.rescue_subsequence.run(trial_num)

            # support graceful termination
            with parallel:
                self.check_termination()
                self.core.break_realtime()
