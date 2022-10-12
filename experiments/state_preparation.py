import numpy as np
from artiq.experiment import *

_DMA_HANDLE_RESET = "state_preparation_reset"
_DMA_HANDLE_LOOP = "state_preparation_loop"
_DMA_HANDLE_SHELVING = "state_preparation_shelving"
_DMA_HANDLE_READOUT = "state_preparation_readout"


class StatePreparation(EnvExperiment):
    """
    State Preparation
    Measures the proportion of population in the correct Zeeman level.
    """

    # kernel_invariants = {}
    global_parameters = [
        "pmt_input_channel",
        "pmt_gating_edge",
        "time_profileswitch_delay_us",
        "time_repump_qubit_us",
        "time_doppler_cooling_us",
        "time_readout_us",
        "dds_board_num",
        "dds_pump_channel",
        "dds_repump_cooling_channel",
        "freq_redist_mhz",
        "freq_pump_cooling_mhz",
        "freq_pump_readout_mhz",
        "freq_repump_cooling_mhz",
        "ampl_pump_cooling_pct",
        "ampl_pump_readout_pct",
        "ampl_repump_cooling_pct",
        "pmt_discrimination"
    ]

    def build(self):
        """
        Set devices and arguments for the experiment.
        """
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # experiment runs
        self.setattr_argument("repetitions",                    NumberValue(default=400, ndecimals=0, step=1, min=1, max=10000))
        self.setattr_argument("num_loops",                      NumberValue(default=1, ndecimals=0, step=1, min=1, max=20))

        # qubit (729) waveforms
        self.setattr_argument("freq_qubit_loop_mhz",            NumberValue(default=92, ndecimals=5, step=1, min=10, max=200))
        self.setattr_argument("freq_qubit_shelving_mhz",        NumberValue(default=92, ndecimals=5, step=1, min=10, max=200))

        self.setattr_argument("ampl_qubit_loop_pct",            NumberValue(default=50, ndecimals=5, step=1, min=10, max=100))
        self.setattr_argument("ampl_qubit_shelving_pct",        NumberValue(default=50, ndecimals=5, step=1, min=10, max=100))

        # timing
        self.setattr_argument("time_loop_us",                   NumberValue(default=10, ndecimals=5, step=1, min=1, max=1000000))
        self.setattr_argument("time_shelving_us",               NumberValue(default=5, ndecimals=5, step=1, min=1, max=1000000))

        # get global parameters
        for param_name in self.global_parameters:
            self.setattr_dataset(param_name, archive=True)


    def prepare(self):
        """
        Set up the dataset and prepare things such that
        the kernel functions have minimal overhead.
        """
        # PMT devices
        self.pmt_counter =                                      self.get_device("ttl_counter{:d}".format(self.pmt_input_channel))
        self.pmt_gating_edge =                                  getattr(self.pmt_counter, 'gate_{:s}_mu'.format(self.pmt_gating_edge))

        # convert time values to machine units
        self.time_profileswitch_delay_mu =                      self.core.seconds_to_mu(self.time_profileswitch_delay_us * us)
        self.time_repump_qubit_mu =                             self.core.seconds_to_mu(self.time_repump_qubit_us * us)
        self.time_cooling_mu =                                  self.core.seconds_to_mu(self.time_doppler_cooling_us * us)
        self.time_readout_mu =                                  self.core.seconds_to_mu(self.time_readout_us * us)

        # DDS devices
        self.dds_board =                                        self.get_device("urukul{:d}_cpld".format(self.dds_board_num))
        self.dds_qubit_board =                                  self.get_device("urukul{:d}_cpld".format(self.dds_board_qubit_num))

        self.dds_pump =                                         self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_pump_channel))
        self.dds_repump_cooling =                               self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_repump_cooling_channel))
        self.dds_repump_qubit =                                 self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_repump_qubit_channel))
        self.dds_qubit =                                        self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_qubit_num, self.dds_qubit_channel))

        # convert frequency to ftw
        self.freq_pump_cooling_ftw =                            self.dds_qubit.frequency_to_ftw(self.freq_pump_cooling_mhz * MHz)
        self.freq_pump_readout_ftw =                            self.dds_qubit.frequency_to_ftw(self.freq_pump_readout_mhz * MHz)
        self.freq_repump_cooling_ftw =                          self.dds_qubit.frequency_to_ftw(self.freq_repump_cooling_mhz * MHz)
        self.freq_repump_qubit_ftw =                            self.dds_qubit.frequency_to_ftw(self.freq_repump_qubit_mhz * MHz)

        # convert amplitude to asf
        self.ampl_pump_cooling_asf =                            self.dds_qubit.amplitude_to_asf(self.ampl_pump_cooling_pct / 100)
        self.ampl_pump_readout_asf =                            self.dds_qubit.amplitude_to_asf(self.ampl_pump_readout_pct / 100)
        self.ampl_repump_cooling_asf =                          self.dds_qubit.amplitude_to_asf(self.ampl_repump_cooling_pct / 100)
        self.ampl_repump_qubit_asf =                            self.dds_qubit.amplitude_to_asf(self.ampl_repump_qubit_pct / 100)

        # novel parameters
        self.freq_qubit_loop_ftw =                              self.dds_qubit.frequency_to_ftw(self.freq_qubit_loop_mhz)
        self.freq_qubit_shelving_ftw =                          self.dds_qubit.frequency_to_ftw(self.freq_qubit_shelving_mhz)

        self.ampl_qubit_loop_asf =                              self.dds_qubit.amplitude_to_asf(self.ampl_qubit_loop_pct / 100)
        self.ampl_qubit_shelving_asf =                          self.dds_qubit.amplitude_to_asf(self.ampl_qubit_shelving_pct / 100)

        self.time_loop_mu =                                     self.core.seconds_to_mu(self.time_loop_us * us)
        self.time_shelving_mu =                                 self.core.seconds_to_mu(self.time_shelving_us * us)

        # set up datasets
        self.set_dataset("state_preparation", np.zeros(self.repetitions))
        self.setattr_dataset("state_preparation")


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
        handle_reset = self.core_dma.get_handle(_DMA_HANDLE_RESET)
        handle_loop = self.core_dma.get_handle(_DMA_HANDLE_LOOP)
        handle_shelving = self.core_dma.get_handle(_DMA_HANDLE_SHELVING)
        handle_readout = self.core_dma.get_handle(_DMA_HANDLE_READOUT)
        self.core.break_realtime()

        # MAIN SEQUENCE
        for trial_num in range(self.repetitions):

            # reset ion using 854 then doppler cool
            self.core_dma.playback_handle(handle_reset)

            # run X spin polarization loops
            for i in range(self.num_loops):
                self.core_dma.playback_handle(handle_loop)

            # shelve
            self.core_dma.playback_handle(handle_shelving)

            # readout
            self.core_dma.playback_handle(handle_readout)

            # update dataset
            with parallel:
                self.mutate_dataset("state_preparation", trial_num, self.pmt_counter.fetch_count())
                self.core.break_realtime()

        # reset after experiment
        self.dds_board.cfg_switches(0b1110)
        self.dds_qubit.cfg_sw(0)


    @kernel(flags={"fast-math"})
    def DMArecord(self):
        """
        Record onto core DMA the AOM sequence for a single data point.
        """
        # reset sequence
        with self.core_dma.record(_DMA_HANDLE_RESET):
            with sequential:
                # set loop waveform
                with parallel:
                    self.dds_qubit_board.set_profile(0)
                    delay_mu(self.time_profileswitch_delay_mu)

                # qubit repump (854) pulse
                self.dds_board.cfg_switches(0b1100)
                delay_mu(self.time_repump_qubit_mu)
                self.dds_board.cfg_switches(0b0100)

                # set cooling waveform
                with parallel:
                    self.dds_board.set_profile(0)
                    delay_mu(self.time_profileswitch_delay_mu)

                # doppler cooling
                self.dds_board.cfg_switches(0b0110)
                delay_mu(self.time_cooling_mu)
                self.dds_board.cfg_switches(0b0100)

        # spin polarization loop sequence
        with self.core_dma.record(_DMA_HANDLE_LOOP):
            with sequential:
                # 729 pulse
                self.dds_qubit.cfg_sw(1)
                delay_mu(self.time_loop_mu)
                self.dds_qubit.cfg_sw(0)

                # 854 repump pulse
                self.dds_board.cfg_switches(0b1100)
                delay_mu(self.time_repump_qubit_mu)
                self.dds_board.cfg_switches(0b0100)

        # shelving sequence
        with self.core_dma.record(_DMA_HANDLE_SHELVING):
            with sequential:
                # set 729 waveform for shelving
                self.dds_qubit.set_mu(self.freq_qubit_shelving_ftw, asf=self.ampl_qubit_shelving_asf)

                # set shelving waveform
                with parallel:
                    self.dds_qubit_board.set_profile(1)
                    delay_mu(self.time_profileswitch_delay_mu)

                # shelving pulse
                self.dds_.cfg_switches(0b0110)
                delay_mu(self.time_shelving_mu)
                self.dds_board.cfg_switches(0b0100)

        # readout sequence
        with self.core_dma.record(_DMA_HANDLE_READOUT):
            with sequential:
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
        self.dds_pump.set_mu(self.freq_pump_cooling_ftw, asf=self.ampl_pump_cooling_asf, profile=0)
        self.dds_pump.set_mu(self.freq_pump_readout_ftw, asf=self.ampl_pump_readout_asf, profile=1)
        self.core.break_realtime()

        self.dds_repump_cooling.set_mu(self.freq_repump_cooling_ftw, asf=self.ampl_repump_cooling_asf, profile=0)
        self.dds_repump_cooling.set_mu(self.freq_repump_cooling_ftw, asf=self.ampl_repump_cooling_asf, profile=1)
        self.core.break_realtime()

        self.dds_repump_qubit.set_mu(self.freq_repump_qubit_ftw, asf=self.ampl_repump_qubit_asf, profile=0)
        self.dds_repump_qubit.set_mu(self.freq_repump_qubit_ftw, asf=self.ampl_repump_qubit_asf, profile=1)
        self.core.break_realtime()

        self.dds_qubit.set_mu(self.freq_qubit_loop_ftw, asf=self.ampl_qubit_loop_asf, profile=0)
        self.dds_qubit.set_mu(self.freq_qubit_shelving_ftw, asf=self.ampl_qubit_shelving_asf, profile=1)
        self.core.break_realtime()


    def analyze(self):
        """
        Analyze the results from the experiment.
        """
        print("counts: {:.3f} +/- {:.3f}".format(np.mean(self.state_preparation), np.std(self.state_preparation)))
