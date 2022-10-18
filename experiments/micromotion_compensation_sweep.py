import numpy as np
from time import sleep
from artiq.experiment import *

_DMA_HANDLE_MM_COMP = "micromotion_compensation_sweep_sequence"

import labrad
from os import environ
from EGGS_labrad.config.dc_config import dc_config


class MicromotionCompensationSweep(EnvExperiment):
    """
    Micromotion Compensation Sweep
    Adjusts the endcap voltages and compares the Rabi frequencies
        of the micromotion sideband with that of the carrier.
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
        "freq_repump_cooling_mhz",
        "freq_repump_qubit_mhz",
        "freq_qubit_mhz",

        "ampl_redist_pct",
        "ampl_pump_cooling_pct",
        "ampl_pump_readout_pct",
        "ampl_repump_cooling_pct",
        "ampl_repump_qubit_pct",
        "ampl_qubit_pct",

        "pmt_discrimination"
    ]

    def build(self):
        """
        Set devices and arguments for the experiment.
        """
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # experiment runs
        self.setattr_argument("repetitions",                        NumberValue(default=20, ndecimals=0, step=1, min=1, max=10000))

        # qubit time scan
        self.setattr_argument("time_rabi_us",                       NumberValue(default=8, ndecimals=5, step=1, min=1, max=1000000))
        self.setattr_argument("time_dc_set_s",                      NumberValue(default=5, ndecimals=3, step=1, min=0, max=100))

        # AOM values
        self.setattr_argument("freq_micromotion_mhz",               NumberValue(default=19, ndecimals=5, step=1, min=1, max=10000))

        # voltage values
        #self.setattr_argument("dc_micromotion_channels",            NumberValue(default=19, ndecimals=5, step=1, min=1, max=10000))
        self.dc_micromotion_channeldict =                           dc_config.channeldict
        self.setattr_argument("dc_micromotion_channels",            EnumerationValue(list(self.dc_micromotion_channeldict.keys())))
        self.setattr_argument("dc_micromotion_voltages_v",          Scannable(
                                                                        default=RangeScan(30, 60, 31),
                                                                        global_min=0, global_max=1000, global_step=1,
                                                                        unit="V", scale=1, ndecimals=4
                                                                    ))


        # get global parameters
        for param_name in self.global_parameters:
            self.setattr_dataset(param_name, archive=True)

        # connect to labrad
        self.cxn = labrad.connect(environ['LABRADHOST'], port=7682, tls_mode='off', username='', password='lab')
        self.dc = self.cxn.dc_server


    def prepare(self):
        """
        Set up the dataset and prepare things such that
        the kernel functions have minimal overhead.
        """
        # PMT devices
        self.pmt_counter =                      self.get_device("ttl_counter{:d}".format(self.pmt_input_channel))
        self.pmt_gating_edge =                  getattr(self.pmt_counter, 'gate_{:s}_mu'.format(self.pmt_gating_edge))

        # convert time values to machine units
        self.time_profileswitch_delay_mu =      self.core.seconds_to_mu(self.time_profileswitch_delay_us * us)
        self.time_repump_qubit_mu =             self.core.seconds_to_mu(self.time_repump_qubit_us * us)
        self.time_cooling_mu =                  self.core.seconds_to_mu(self.time_doppler_cooling_us * us)
        self.time_readout_mu =                  self.core.seconds_to_mu(self.time_readout_us * us)
        self.time_redist_mu =                   self.core.seconds_to_mu(self.time_redist_us * us)

        # rabi flopping timing
        self.time_rabi_mu =                     self.core.seconds_to_mu(self.time_rabi_us * us)

        # DDS devices
        self.dds_board =                        self.get_device("urukul{:d}_cpld".format(self.dds_board_num))
        self.dds_qubit_board =                  self.get_device("urukul{:d}_cpld".format(self.dds_board_qubit_num))

        self.dds_probe =                        self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_probe_channel))
        self.dds_pump =                         self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_pump_channel))
        self.dds_repump_cooling =               self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_repump_cooling_channel))
        self.dds_repump_qubit =                 self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_repump_qubit_channel))
        self.dds_qubit =                        self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_qubit_num, self.dds_qubit_channel))

        # convert frequency to ftw
        self.freq_redist_ftw =                  self.dds_qubit.frequency_to_ftw(self.freq_redist_mhz * MHz)
        self.freq_pump_cooling_ftw =            self.dds_qubit.frequency_to_ftw(self.freq_pump_cooling_mhz * MHz)
        self.freq_pump_readout_ftw =            self.dds_qubit.frequency_to_ftw(self.freq_pump_readout_mhz * MHz)
        self.freq_repump_cooling_ftw =          self.dds_qubit.frequency_to_ftw(self.freq_repump_cooling_mhz * MHz)
        self.freq_repump_qubit_ftw =            self.dds_qubit.frequency_to_ftw(self.freq_repump_qubit_mhz * MHz)
        self.freq_qubit_ftw =                   self.dds_qubit.frequency_to_ftw(self.freq_qubit_mhz * MHz)
        self.freq_bsb_ftw =                     self.dds_qubit.frequency_to_ftw((0.5 * self.freq_micromotion_mhz + self.freq_qubit_mhz) * MHz)

        # convert amplitude to asf
        self.ampl_redist_asf =                  self.dds_qubit.amplitude_to_asf(self.ampl_redist_pct / 100)
        self.ampl_pump_cooling_asf =            self.dds_qubit.amplitude_to_asf(self.ampl_pump_cooling_pct / 100)
        self.ampl_pump_readout_asf =            self.dds_qubit.amplitude_to_asf(self.ampl_pump_readout_pct / 100)
        self.ampl_repump_cooling_asf =          self.dds_qubit.amplitude_to_asf(self.ampl_repump_cooling_pct / 100)
        self.ampl_repump_qubit_asf =            self.dds_qubit.amplitude_to_asf(self.ampl_repump_qubit_pct / 100)
        self.ampl_qubit_asf =                   self.dds_qubit.amplitude_to_asf(self.ampl_qubit_pct / 100)

        # get voltage parameters
        self.dc_micromotion_channels =          self.dc_micromotion_channeldict["dc_micromotion_channels"]
        self.dc_micromotion_voltages_v =        list(self.dc_micromotion_voltages_v)

        # set up datasets
        self.set_dataset("micromotion_compensation", np.zeros([self.repetitions, 2]))
        self.setattr_dataset("micromotion_compensation")


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
        handle_sequence = self.core_dma.get_handle(_DMA_HANDLE_MM_COMP)
        self.core.break_realtime()

        # MAIN SEQUENCE
        for voltage_v in self.dc_micromotion_voltages_v:

            # set voltage
            self.voltage_set(self.dc_micromotion_channels, voltage_v)

            # repeat experiment
            for trial_num in range(self.repetitions):

                # set carrier waveform
                with parallel:
                    self.dds_board.set_profile(0)
                    delay_mu(self.time_profileswitch_delay_mu)

                # run sequence
                self.core_dma.playback_handle(handle_sequence)

                # set bsb waveform
                with parallel:
                    self.dds_board.set_profile(1)
                    delay_mu(self.time_profileswitch_delay_mu)

                # run sequence
                self.core_dma.playback_handle(handle_sequence)

                # add data to dataset
                self.mutate_dataset("micromotion_compensation", trial_num, [self.pmt_counter.fetch_count(), self.pmt_counter.fetch_count()])

        # reset board profiles
        self.dds_board.set_profile(0)
        self.dds_qubit_board.set_profile(0)

        # reset after experiment
        self.dds_board.cfg_switches(0b1110)
        self.dds_qubit.cfg_sw(0)


    @kernel(flags={"fast-math"})
    def DMArecord(self):
        """
        Record onto core DMA the AOM sequence for a single data point.
        """
        # compensation sequence
        with self.core_dma.record(_DMA_HANDLE_MM_COMP):
            # qubit repump pulse
            self.dds_board.cfg_switches(0b1100)
            delay_mu(self.time_repump_qubit_mu)
            self.dds_board.cfg_switches(0b0100)

            # set cooling waveform
            with parallel:
                self.dds_board.set_profile(0)
                delay_mu(self.time_profileswitch_delay_mu)

            # do cooling
            self.dds_board.cfg_switches(0b0110)
            delay_mu(self.time_cooling_mu)
            self.dds_board.cfg_switches(0b0100)

            # do spin depolarization using probe
            self.dds_board.cfg_switches(0b0101)
            delay_mu(self.time_redist_mu)
            self.dds_board.cfg_switches(0b0100)

            # set readout waveform
            with parallel:
                self.dds_board.set_profile(1)
                delay_mu(self.time_profileswitch_delay_mu)

            # rabi flopping w/qubit laser
            self.dds_qubit.cfg_sw(1)
            delay_mu(self.time_rabi_mu)
            self.dds_qubit.cfg_sw(0)

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
        self.core.break_realtime()

        self.dds_pump.set_mu(self.freq_pump_cooling_ftw, asf=self.ampl_pump_cooling_asf, profile=0)
        self.dds_pump.set_mu(self.freq_pump_readout_ftw, asf=self.ampl_pump_readout_asf, profile=1)
        self.core.break_realtime()

        self.dds_repump_cooling.set_mu(self.freq_repump_cooling_ftw, asf=self.ampl_repump_cooling_asf, profile=0)
        self.dds_repump_cooling.set_mu(self.freq_repump_cooling_ftw, asf=self.ampl_repump_cooling_asf, profile=1)
        self.core.break_realtime()

        self.dds_repump_qubit.set_mu(self.freq_repump_qubit_ftw, asf=self.ampl_repump_qubit_asf, profile=0)
        self.dds_repump_qubit.set_mu(self.freq_repump_qubit_ftw, asf=self.ampl_repump_qubit_asf, profile=1)
        self.core.break_realtime()

        # profile 0 is carrier, profile 1 is bsb
        self.dds_qubit.set_mu(self.freq_qubit_ftw, asf=self.ampl_qubit_asf, profile=0)
        self.dds_qubit.set_mu(self.freq_bsb_ftw, asf=self.ampl_qubit_asf, profile=1)
        self.core.break_realtime()


    @rpc(flags={"async"})
    def voltage_set(self, channel, voltage_v):
        """
        Set the channel to the desired voltage.
        """
        voltage_set_v = self.dc.voltage(channel, voltage_v)
        sleep(self.time_dc_set_s)


    def analyze(self):
        """
        Analyze the results from the experiment.
        """
        print("carrier counts: {:f} +/- {:f}".format(np.mean(self.micromotion_compensation[:, 0]), np.std(self.micromotion_compensation[:, 0])))
        print("bsb counts: {:f} +/- {:f}".format(np.mean(self.micromotion_compensation[:, 1]), np.std(self.micromotion_compensation[:, 1])))
