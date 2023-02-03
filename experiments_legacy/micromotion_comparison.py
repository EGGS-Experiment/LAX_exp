import labrad
import numpy as np
from os import environ
from time import sleep

from artiq.experiment import *
from EGGS_labrad.config.dc_config import dc_config

_DMA_HANDLE_SEQUENCE = "micromotion_comparison_sequence"


class MicromotionComparison(EnvExperiment):
    """
    Micromotion Comparison

    Compare the ratio of the Rabi frequencies of the carrier and the micromotion sideband.
    """

    # kernel_invariants = {}
    global_parameters = [
        "pmt_input_channel",
        "pmt_gating_edge",

        "time_profileswitch_delay_us",
        "time_repump_qubit_us",
        "time_doppler_cooling_us",
        "time_readout_us",
        "time_redist_us",

        "dds_board_num",
        "dds_board_qubit_num",
        "dds_probe_channel",
        "dds_pump_channel",
        "dds_repump_cooling_channel",
        "dds_repump_qubit_channel",
        "dds_qubit_channel",

        "freq_redist_mhz",
        "freq_pump_cooling_mhz",
        "freq_pump_readout_mhz",
        "freq_pump_rescue_mhz",
        "freq_repump_cooling_mhz",
        "freq_repump_qubit_mhz",

        "ampl_redist_pct",
        "ampl_pump_cooling_pct",
        "ampl_pump_readout_pct",
        "ampl_pump_rescue_pct",
        "ampl_repump_cooling_pct",
        "ampl_repump_qubit_pct",
        "ampl_qubit_pct"
    ]

    def build(self):
        """
        Set devices and arguments for the experiment.
        """
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # experiment runs
        self.setattr_argument("repetitions",                            NumberValue(default=50, ndecimals=0, step=1, min=1, max=10000))

        # additional cooling
        self.setattr_argument("repetitions_per_cooling",                NumberValue(default=1, ndecimals=0, step=1, min=1, max=10000))
        self.setattr_argument("additional_cooling_time_s",              NumberValue(default=1, ndecimals=5, step=0.1, min=0, max=10000))

        # qubit parameters
        self.setattr_argument("freq_carrier_mhz",                       NumberValue(default=104.335, ndecimals=6, step=1, min=1, max=10000))
        self.setattr_argument("freq_micromotion_mhz",                   NumberValue(default=94.635, ndecimals=6, step=1, min=1, max=10000))
        self.setattr_argument("time_pipulse_us",                        NumberValue(default=250, ndecimals=5, step=1, min=1, max=10000))

        # shimming parameters
        self.dc_micromotion_channeldict =                               dc_config.channeldict
        self.setattr_argument("dc_micromotion_channel_1",               EnumerationValue(list(self.dc_micromotion_channeldict.keys()), default='V Shim'))
        self.setattr_argument("dc_micromotion_voltages_1_v_list",       Scannable(
                                                                            default=CenterScan(100.0, 10.0, 0.1, randomize=True),
                                                                            global_min=0, global_max=1000, global_step=1,
                                                                            unit="V", scale=1, ndecimals=2
                                                                        ))
        self.setattr_argument("dc_micromotion_channel_2",               EnumerationValue(list(self.dc_micromotion_channeldict.keys()), default='W Endcap'))
        self.setattr_argument("dc_micromotion_voltages_2_v_list",       Scannable(
                                                                            default=CenterScan(250.0, 5.0, 0.1, randomize=True),
                                                                            global_min=0, global_max=1000, global_step=1,
                                                                            unit="V", scale=1, ndecimals=2
                                                                        ))

        # get global parameters
        for param_name in self.global_parameters:
            self.setattr_dataset(param_name, archive=True)


    def prepare(self):
        """
        Set up the dataset and prepare things such that
        the kernel functions have minimal overhead.
        """
        # PMT devices
        self.pmt_counter =                                              self.get_device("ttl{:d}_counter".format(self.pmt_input_channel))
        self.pmt_gating_edge =                                          getattr(self.pmt_counter, 'gate_{:s}_mu'.format(self.pmt_gating_edge))

        # convert time values to machine units
        self.time_profileswitch_delay_mu =                              self.core.seconds_to_mu(self.time_profileswitch_delay_us * us)
        self.time_rfswitch_delay_mu =                                   self.core.seconds_to_mu(2 * us)
        self.time_repump_qubit_mu =                                     self.core.seconds_to_mu(self.time_repump_qubit_us * us)
        self.time_cooling_mu =                                          self.core.seconds_to_mu(self.time_doppler_cooling_us * us)
        self.time_readout_mu =                                          self.core.seconds_to_mu(self.time_readout_us * us)
        self.time_redist_mu =                                           self.core.seconds_to_mu(self.time_redist_us * us)

        # DDS devices
        self.dds_board =                                                self.get_device("urukul{:d}_cpld".format(self.dds_board_num))
        self.dds_qubit_board =                                          self.get_device("urukul{:d}_cpld".format(self.dds_board_qubit_num))

        self.dds_probe =                                                self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_probe_channel))
        self.dds_pump =                                                 self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_pump_channel))
        self.dds_repump_cooling =                                       self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_repump_cooling_channel))
        self.dds_repump_qubit =                                         self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_repump_qubit_channel))
        self.dds_qubit =                                                self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_qubit_num, self.dds_qubit_channel))

        # RF switches
        self.dds_repump_qubit_switch =                                  self.get_device("ttl21")

        # convert frequency to ftw
        self.freq_redist_ftw =                                          self.dds_qubit.frequency_to_ftw(self.freq_redist_mhz * MHz)
        self.freq_pump_cooling_ftw =                                    self.dds_qubit.frequency_to_ftw(self.freq_pump_cooling_mhz * MHz)
        self.freq_pump_readout_ftw =                                    self.dds_qubit.frequency_to_ftw(self.freq_pump_readout_mhz * MHz)
        self.freq_pump_rescue_ftw =                                     self.dds_qubit.frequency_to_ftw(self.freq_pump_rescue_mhz * MHz)
        self.freq_repump_cooling_ftw =                                  self.dds_qubit.frequency_to_ftw(self.freq_repump_cooling_mhz * MHz)
        self.freq_repump_qubit_ftw =                                    self.dds_qubit.frequency_to_ftw(self.freq_repump_qubit_mhz * MHz)

        # convert amplitude to asf
        self.ampl_redist_asf =                                          self.dds_qubit.amplitude_to_asf(self.ampl_redist_pct / 100)
        self.ampl_pump_cooling_asf =                                    self.dds_qubit.amplitude_to_asf(self.ampl_pump_cooling_pct / 100)
        self.ampl_pump_readout_asf =                                    self.dds_qubit.amplitude_to_asf(self.ampl_pump_readout_pct / 100)
        self.ampl_pump_rescue_asf =                                     self.dds_qubit.amplitude_to_asf(self.ampl_pump_rescue_pct / 100)
        self.ampl_repump_cooling_asf =                                  self.dds_qubit.amplitude_to_asf(self.ampl_repump_cooling_pct / 100)
        self.ampl_repump_qubit_asf =                                    self.dds_qubit.amplitude_to_asf(self.ampl_repump_qubit_pct / 100)
        self.ampl_qubit_asf =                                           self.dds_qubit.amplitude_to_asf(self.ampl_qubit_pct / 100)


        # convert values for rabi flopping
        self.freq_carrier_ftw =                                         self.dds_qubit.frequency_to_ftw(self.freq_carrier_mhz * MHz)
        self.freq_micromotion_ftw =                                     self.dds_qubit.frequency_to_ftw(self.freq_micromotion_mhz * MHz)
        self.time_pipulse_mu =                                          self.core.seconds_to_mu(self.time_pipulse_us * us)

        # get voltage parameters
        self.dc_micromotion_voltages_1_v_list =                         np.array(list(self.dc_micromotion_voltages_1_v_list))
        self.dc_micromotion_channel_1 =                                 self.dc_micromotion_channeldict[self.dc_micromotion_channel_1]['num']

        self.dc_micromotion_voltages_2_v_list =                         np.array(list(self.dc_micromotion_voltages_2_v_list))
        self.dc_micromotion_channel_2 =                                 self.dc_micromotion_channeldict[self.dc_micromotion_channel_2]['num']


        # connect to labrad
        self.cxn =                                                      labrad.connect(environ['LABRADHOST'], port=7682, tls_mode='off', username='', password='lab')
        self.dc =                                                       self.cxn.dc_server

        # set up datasets
        self._iter_dataset =                                            0
        self.set_dataset("micromotion_comparison",                      np.zeros((
                                                                            len(self.dc_micromotion_voltages_1_v_list) * len(self.dc_micromotion_voltages_2_v_list) * self.repetitions,
                                                                            4
                                                                        )))
        self.setattr_dataset("micromotion_comparison")


    @kernel(flags={"fast-math"})
    def run(self):
        """
        Run the experimental sequence.
        """
        self.core.reset()

        # prepare devices
        self.prepareDevices()

        # record dma and get handle
        self.DMArecord()
        handle_sequence = self.core_dma.get_handle(_DMA_HANDLE_SEQUENCE)
        self.core.break_realtime()

        # MAIN SEQUENCE
        for trial_num in range(self.repetitions):

            # sweep voltage 1
            for voltage_1_v in self.dc_micromotion_voltages_1_v_list:

                # set channel 1 voltage
                self.voltage_set(self.dc_micromotion_channel_1, voltage_1_v)
                self.core.break_realtime()

                # sweep voltage 2
                for voltage_2_v in self.dc_micromotion_voltages_2_v_list:

                    # set channel 2 voltage
                    self.voltage_set(self.dc_micromotion_channel_2, voltage_2_v)
                    self.core.break_realtime()

                    # rabi flop on carrier
                    self.dds_qubit_board.set_profile(0)
                    self.core_dma.playback_handle(handle_sequence)

                    # get carrier counts
                    counts_carrier = self.pmt_counter.fetch_count()
                    self.core.break_realtime()

                    # rabi flop on sideband
                    self.dds_qubit_board.set_profile(1)
                    self.core_dma.playback_handle(handle_sequence)

                    # get sideband counts
                    counts_sideband = self.pmt_counter.fetch_count()
                    self.core.break_realtime()

                    # update dataset
                    with parallel:
                        self.update_dataset(voltage_1_v, voltage_2_v, counts_carrier, counts_sideband)
                        self.core.break_realtime()

            # add post repetition cooling
            if (trial_num > 0) and (trial_num % self.repetitions_per_cooling == 0):
                # set rescue waveform
                with parallel:
                    self.dds_board.set_profile(2)
                    delay_mu(self.time_profileswitch_delay_mu)

                # start rescuing
                self.dds_board.io_update.pulse_mu(8)
                self.dds_board.cfg_switches(0b0110)
                delay(self.additional_cooling_time_s)
                self.dds_board.cfg_switches(0b0100)

        # reset board profiles
        self.dds_pump.set_mu(self.freq_pump_rescue_ftw, asf=self.ampl_pump_rescue_asf, profile=0)
        self.core.break_realtime()
        self.dds_board.set_profile(0)
        self.dds_qubit_board.set_profile(0)

        # reset AOMs after experiment
        self.dds_board.cfg_switches(0b1110)
        self.dds_qubit.cfg_sw(False)
        self.dds_repump_qubit_switch.on()


    @kernel(flags={"fast-math"})
    def DMArecord(self):
        """
        Record onto core DMA the AOM sequence for a single data point.
        """
        # complete rabi flopping sequence
        with self.core_dma.record(_DMA_HANDLE_SEQUENCE):
            with sequential:
                # enable 854 rf switch
                # with parallel:
                self.dds_repump_qubit_switch.on()
                delay_mu(self.time_rfswitch_delay_mu)
                self.dds_board.cfg_switches(0b1100)

                # qubit repump (854) pulse
                delay_mu(self.time_repump_qubit_mu)
                # with parallel:
                self.dds_repump_qubit_switch.off()
                self.dds_board.cfg_switches(0b0100)
                delay_mu(self.time_rfswitch_delay_mu)

                # set cooling waveform
                with parallel:
                    self.dds_board.set_profile(0)
                    delay_mu(self.time_profileswitch_delay_mu)

                # do doppler cooling
                self.dds_board.cfg_switches(0b0110)
                delay_mu(self.time_cooling_mu)
                self.dds_board.cfg_switches(0b0100)

                # do spin depolarization using probe
                self.dds_board.cfg_switches(0b0101)
                delay_mu(self.time_redist_mu)
                self.dds_board.cfg_switches(0b0100)

                # rabi flop
                self.dds_qubit.cfg_sw(True)
                delay_mu(self.time_pipulse_mu)
                self.dds_qubit.cfg_sw(False)

                # set readout waveform
                with parallel:
                    self.dds_board.set_profile(1)
                    delay_mu(self.time_profileswitch_delay_mu)

                # readout pulse
                self.dds_board.cfg_switches(0b0110)
                self.pmt_gating_edge(self.time_readout_mu)
                self.dds_board.cfg_switches(0b0100)


    @kernel(flags={"fast-math"})
    def prepareDevices(self):
        """
        Prepare devices for the experiment.
        """
        self.core.break_realtime()

        # set AOM DDS waveforms; profile 0 is cooling, profile 1 is readout
        self.dds_probe.set_mu(self.freq_redist_ftw, asf=self.ampl_redist_asf, profile=0)
        self.dds_probe.set_mu(self.freq_redist_ftw, asf=self.ampl_redist_asf, profile=1)
        self.dds_probe.set_mu(self.freq_redist_ftw, asf=self.ampl_redist_asf, profile=2)
        self.core.break_realtime()

        self.dds_pump.set_mu(self.freq_pump_cooling_ftw, asf=self.ampl_pump_cooling_asf, profile=0)
        self.dds_pump.set_mu(self.freq_pump_readout_ftw, asf=self.ampl_pump_readout_asf, profile=1)
        self.dds_pump.set_mu(self.freq_pump_rescue_ftw, asf=self.ampl_pump_rescue_asf, profile=2)
        self.core.break_realtime()

        self.dds_repump_cooling.set_mu(self.freq_repump_cooling_ftw, asf=self.ampl_repump_cooling_asf, profile=0)
        self.dds_repump_cooling.set_mu(self.freq_repump_cooling_ftw, asf=self.ampl_repump_cooling_asf, profile=1)
        self.dds_repump_cooling.set_mu(self.freq_repump_cooling_ftw, asf=self.ampl_repump_cooling_asf, profile=2)
        self.core.break_realtime()

        self.dds_repump_qubit.set_mu(self.freq_repump_qubit_ftw, asf=self.ampl_repump_qubit_asf, profile=0)
        self.dds_repump_qubit.set_mu(self.freq_repump_qubit_ftw, asf=self.ampl_repump_qubit_asf, profile=1)
        self.dds_repump_qubit.set_mu(self.freq_repump_qubit_ftw, asf=self.ampl_repump_qubit_asf, profile=2)
        self.core.break_realtime()

        self.dds_qubit.set_mu(self.freq_carrier_ftw, asf=self.ampl_qubit_asf, profile=0)
        self.dds_qubit.set_mu(self.freq_micromotion_ftw, asf=self.ampl_qubit_asf, profile=1)
        self.core.break_realtime()

        # set rf switches
        self.dds_repump_qubit_switch.off()


    #@rpc
    def voltage_set(self, channel, voltage_v):
        """
        Set the DC voltage channel to the desired voltage.
        """
        voltage_set_v = self.dc.voltage(channel, voltage_v)
        #sleep(0.1)
        #print('\tvoltage set: {}'.format(voltage_set_v))


    @rpc(flags={"async"})
    def update_dataset(self, voltage_1, voltage_2, pmt_counts_carrier, pmt_counts_sideband):
        """
        Records values via rpc to minimize kernel overhead.
        """
        data_tmp = np.array([voltage_1, voltage_2, pmt_counts_carrier, pmt_counts_sideband])
        self.mutate_dataset('micromotion_comparison', self._iter_dataset, data_tmp)
        self._iter_dataset += 1


    def analyze(self):
        """
        Analyze the results from the experiment.
        """
        pass
