import numpy as np
from artiq.experiment import *

_DMA_HANDLE_INITIALIZATION =            "rabi_flopping_initialization"
_DMA_HANDLE_READOUT =                   "rabi_flopping_readout"

_URUKUL0_PROFILE_DEFAULT =              0
_URUKUL0_PROFILE_QUBIT =                1
_URUKUL0_PROFILE_OFF =                  2

_URUKUL1_PROFILE_DEFAULT =              0
_URUKUL1_PROFILE_DOPPLERCOOLING =       1
_URUKUL1_PROFILE_READOUT =              2
_URUKUL1_PROFILE_RESCUE =               3
_URUKUL1_PROFILE_DEPLETION =            4
_URUKUL1_PROFILE_SPINPOLARIZATION =     5


class RabiFloppingFrequencyShiftTest(EnvExperiment):
    """
    Rabi Flopping Freq Shift Test

    Measures ion fluorescence vs 729nm pulse time and frequency.
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

        # todo: add to dataset manager
        "freq_probe_shutter_mhz",
        "freq_pump_shutter_mhz",
        "freq_repump_cooling_shutter_mhz",
        "freq_repump_qubit_shutter_mhz",
        "freq_qubit_shutter_mhz",

        "ampl_redist_pct",
        "ampl_pump_cooling_pct",
        "ampl_pump_readout_pct",
        "ampl_pump_rescue_pct",
        "ampl_repump_cooling_pct",
        "ampl_repump_qubit_pct",
        "ampl_qubit_pct",

        # todo: add to dataset manager
        "ampl_shutter_pct"
    ]

    def build(self):
        """
        Set devices and arguments for the experiment.
        """
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # experiment runs
        self.setattr_argument("repetitions",                            NumberValue(default=100, ndecimals=0, step=1, min=1, max=10000))

        # additional cooling
        self.setattr_argument("repetitions_per_cooling",                NumberValue(default=1, ndecimals=0, step=1, min=1, max=10000))
        self.setattr_argument("additional_cooling_time_s",              NumberValue(default=1, ndecimals=5, step=0.1, min=0, max=10000))

        # qubit parameters
        self.setattr_argument("time_rabi_us_list",                      Scannable(
                                                                            default=RangeScan(0, 400, 401, randomize=True),
                                                                            global_min=1, global_max=100000, global_step=1,
                                                                            unit="us", scale=1, ndecimals=5
                                                                        ))

        # AOM values
        self.setattr_argument("freq_qubit_mhz",                         NumberValue(default=104.335, ndecimals=5, step=1, min=1, max=10000))

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
        self.time_repump_qubit_mu =                                     self.core.seconds_to_mu(self.time_repump_qubit_us * us)
        self.time_cooling_mu =                                          self.core.seconds_to_mu(self.time_doppler_cooling_us * us)
        self.time_readout_mu =                                          self.core.seconds_to_mu(self.time_readout_us * us)
        self.time_redist_mu =                                           self.core.seconds_to_mu(self.time_redist_us * us)

        # rabi flopping timing
        self.time_rabi_mu_list =                                        [self.core.seconds_to_mu(time_us * us) for time_us in self.time_rabi_us_list]
        max_time_us =                                                   np.max(list(self.time_rabi_us_list))
        self.time_delay_mu_list =                                       [self.core.seconds_to_mu((max_time_us - time_us) * us) for time_us in self.time_rabi_us_list]
        self.num_time_points_list =                                     list(range(len(self.time_rabi_mu_list)))

        # DDS devices
        self.dds_board =                                                self.get_device("urukul{:d}_cpld".format(self.dds_board_num))
        self.dds_qubit_board =                                          self.get_device("urukul{:d}_cpld".format(self.dds_board_qubit_num))

        self.dds_probe =                                                self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_probe_channel))
        self.dds_pump =                                                 self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_pump_channel))
        self.dds_repump_cooling =                                       self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_repump_cooling_channel))
        self.dds_repump_qubit =                                         self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_repump_qubit_channel))
        self.dds_qubit =                                                self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_qubit_num, self.dds_qubit_channel))

        # convert frequency to ftw
        self.freq_redist_ftw =                                          self.dds_qubit.frequency_to_ftw(self.freq_redist_mhz * MHz)
        self.freq_pump_cooling_ftw =                                    self.dds_qubit.frequency_to_ftw(self.freq_pump_cooling_mhz * MHz)
        self.freq_pump_readout_ftw =                                    self.dds_qubit.frequency_to_ftw(self.freq_pump_readout_mhz * MHz)
        self.freq_pump_rescue_ftw =                                     self.dds_qubit.frequency_to_ftw(self.freq_pump_rescue_mhz * MHz)
        self.freq_repump_cooling_ftw =                                  self.dds_qubit.frequency_to_ftw(self.freq_repump_cooling_mhz * MHz)
        self.freq_repump_qubit_ftw =                                    self.dds_qubit.frequency_to_ftw(self.freq_repump_qubit_mhz * MHz)
        self.freq_qubit_ftw =                                           self.dds_qubit.frequency_to_ftw(self.freq_qubit_mhz * MHz)

        self.freq_probe_shutter_ftw =                                   self.dds_qubit.frequency_to_ftw(self.freq_probe_shutter_mhz * MHz)
        self.freq_pump_shutter_ftw =                                    self.dds_qubit.frequency_to_ftw(self.freq_pump_shutter_mhz * MHz)
        self.freq_repump_cooling_shutter_ftw =                          self.dds_qubit.frequency_to_ftw(self.freq_repump_cooling_shutter_mhz * MHz)
        self.freq_repump_qubit_shutter_ftw =                            self.dds_qubit.frequency_to_ftw(self.freq_repump_qubit_shutter_mhz * MHz)
        self.freq_qubit_shutter_ftw =                                   self.dds_qubit.frequency_to_ftw(self.freq_qubit_shutter_mhz * MHz)

        # convert amplitude to asf
        self.ampl_redist_asf =                                          self.dds_qubit.amplitude_to_asf(self.ampl_redist_pct / 100)
        self.ampl_pump_cooling_asf =                                    self.dds_qubit.amplitude_to_asf(self.ampl_pump_cooling_pct / 100)
        self.ampl_pump_readout_asf =                                    self.dds_qubit.amplitude_to_asf(self.ampl_pump_readout_pct / 100)
        self.ampl_pump_rescue_asf =                                     self.dds_qubit.amplitude_to_asf(self.ampl_pump_rescue_pct / 100)
        self.ampl_repump_cooling_asf =                                  self.dds_qubit.amplitude_to_asf(self.ampl_repump_cooling_pct / 100)
        self.ampl_repump_qubit_asf =                                    self.dds_qubit.amplitude_to_asf(self.ampl_repump_qubit_pct / 100)
        self.ampl_qubit_asf =                                           self.dds_qubit.amplitude_to_asf(self.ampl_qubit_pct / 100)

        self.ampl_shutter_asf =                                         self.dds_qubit.amplitude_to_asf(self.ampl_shutter_pct / 100)


        # set up datasets
        self.set_dataset("rabi_flopping", [])
        self.setattr_dataset("rabi_flopping")
        self.set_dataset("rabi_flopping_processed", np.zeros([len(self.time_rabi_mu_list), 2]))
        self.setattr_dataset("rabi_flopping_processed")


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
        handle_initialize = self.core_dma.get_handle(_DMA_HANDLE_INITIALIZATION)
        handle_readout = self.core_dma.get_handle(_DMA_HANDLE_READOUT)
        self.core.break_realtime()

        # MAIN SEQUENCE
        for trial_num in range(self.repetitions):

            # sweep time
            for i in self.num_time_points_list:

                # get timings
                time_delay_mu = self.time_delay_mu_list[i]
                time_rabi_mu = self.time_rabi_mu_list[i]

                # run repump and cooling
                self.core_dma.playback_handle(handle_initialize)

                # prepare qubit waveform
                self.dds_qubit_board.set_profile(_URUKUL0_PROFILE_QUBIT)
                delay_mu(self.time_profileswitch_delay_mu)

                # wait given time
                delay_mu(time_delay_mu)

                # rabi flopping w/qubit laser
                self.dds_qubit.cfg_sw(True)
                delay_mu(time_rabi_mu)
                self.dds_qubit.cfg_sw(False)

                # turn off qubit waveform
                self.dds_qubit_board.set_profile(_URUKUL0_PROFILE_OFF)
                delay_mu(self.time_profileswitch_delay_mu)

                # do readout
                self.core_dma.playback_handle(handle_readout)

                # update dataset
                with parallel:
                    self.update_dataset(time_rabi_mu, self.pmt_counter.fetch_count())
                    self.core.break_realtime()

            # add post repetition cooling
            if (trial_num > 0) and (trial_num % self.repetitions_per_cooling == 0):
                # set rescue waveform
                self.dds_board.set_profile(_URUKUL1_PROFILE_RESCUE)
                delay_mu(self.time_profileswitch_delay_mu)
                self.dds_board.io_update.pulse_mu(8)

                # start rescuing
                self.dds_board.cfg_switches(0b0110)
                delay(self.additional_cooling_time_s)
                self.dds_board.cfg_switches(0b0100)

        # reset after experiment
        self.dds_board.cfg_switches(0b1110)
        self.dds_board.set_profile(_URUKUL1_PROFILE_DEFAULT)
        self.dds_board.io_update.pulse_mu(8)

        self.dds_qubit.cfg_sw(False)
        self.dds_qubit_board.set_profile(_URUKUL0_PROFILE_DEFAULT)
        self.dds_qubit_board.io_update.pulse_mu(8)



    @kernel(flags={"fast-math"})
    def DMArecord(self):
        """
        Record onto core DMA the AOM sequence for a single data point.
        """
        # initialization sequence
        with self.core_dma.record(_DMA_HANDLE_INITIALIZATION):
            with sequential:
                # deplete qubit state
                # set depletion waveform
                self.dds_board.set_profile(_URUKUL1_PROFILE_DEPLETION)
                delay_mu(self.time_profileswitch_delay_mu)
                self.dds_board.cfg_switches(0b1100)
                delay_mu(self.time_repump_qubit_mu)
                self.dds_board.cfg_switches(0b0100)

                # doppler cooling
                # set doppler cooling waveform
                self.dds_board.set_profile(_URUKUL1_PROFILE_DOPPLERCOOLING)
                delay_mu(self.time_profileswitch_delay_mu)
                self.dds_board.cfg_switches(0b0110)
                delay_mu(self.time_cooling_mu)
                self.dds_board.cfg_switches(0b0100)

                # spin polarization
                # set spin polarization waveform
                self.dds_board.set_profile(_URUKUL1_PROFILE_SPINPOLARIZATION)
                delay_mu(self.time_profileswitch_delay_mu)
                self.dds_board.cfg_switches(0b0101)
                delay_mu(self.time_redist_mu)
                self.dds_board.cfg_switches(0b0100)

        # readout sequence
        with self.core_dma.record(_DMA_HANDLE_READOUT):
            with sequential:
                # readout
                # set readout waveform
                self.dds_board.set_profile(_URUKUL1_PROFILE_READOUT)
                delay_mu(self.time_profileswitch_delay_mu)
                self.dds_board.cfg_switches(0b0110)
                self.pmt_gating_edge(self.time_readout_mu)
                self.dds_board.cfg_switches(0b0100)


    @kernel(flags={"fast-math"})
    def prepareDevices(self):
        """
        Prepare devices for the experiment.
        """
        self.core.break_realtime()

        # set main AOM DDS waveforms
        self.dds_probe.set_mu(self.freq_probe_shutter_ftw,                  asf=self.ampl_shutter_asf,              profile=_URUKUL1_PROFILE_DOPPLERCOOLING)
        self.dds_probe.set_mu(self.freq_probe_shutter_ftw,                  asf=self.ampl_shutter_asf,              profile=_URUKUL1_PROFILE_READOUT)
        self.dds_probe.set_mu(self.freq_probe_shutter_ftw,                  asf=self.ampl_shutter_asf,              profile=_URUKUL1_PROFILE_RESCUE)
        self.dds_probe.set_mu(self.freq_probe_shutter_ftw,                  asf=self.ampl_shutter_asf,              profile=_URUKUL1_PROFILE_DEPLETION)
        self.dds_probe.set_mu(self.freq_redist_ftw,                         asf=self.ampl_redist_asf,               profile=_URUKUL1_PROFILE_SPINPOLARIZATION)
        self.core.break_realtime()

        self.dds_pump.set_mu(self.freq_pump_cooling_ftw,                    asf=self.ampl_pump_cooling_asf,         profile=_URUKUL1_PROFILE_DOPPLERCOOLING)
        self.dds_pump.set_mu(self.freq_pump_readout_ftw,                    asf=self.ampl_pump_readout_asf,         profile=_URUKUL1_PROFILE_READOUT)
        self.dds_pump.set_mu(self.freq_pump_rescue_ftw,                     asf=self.ampl_pump_rescue_asf,          profile=_URUKUL1_PROFILE_RESCUE)
        self.dds_pump.set_mu(self.freq_pump_shutter_ftw,                    asf=self.ampl_shutter_asf,              profile=_URUKUL1_PROFILE_DEPLETION)
        self.dds_pump.set_mu(self.freq_pump_shutter_ftw,                    asf=self.ampl_shutter_asf,              profile=_URUKUL1_PROFILE_SPINPOLARIZATION)
        self.core.break_realtime()

        self.dds_repump_cooling.set_mu(self.freq_repump_cooling_ftw,        asf=self.ampl_repump_cooling_asf,       profile=_URUKUL1_PROFILE_DOPPLERCOOLING)
        self.dds_repump_cooling.set_mu(self.freq_repump_cooling_ftw,        asf=self.ampl_repump_cooling_asf,       profile=_URUKUL1_PROFILE_READOUT)
        self.dds_repump_cooling.set_mu(self.freq_repump_cooling_ftw,        asf=self.ampl_repump_cooling_asf,       profile=_URUKUL1_PROFILE_RESCUE)
        self.dds_repump_cooling.set_mu(self.freq_repump_cooling_ftw,        asf=self.ampl_repump_cooling_asf,       profile=_URUKUL1_PROFILE_DEPLETION)
        self.dds_repump_cooling.set_mu(self.freq_repump_cooling_ftw,        asf=self.ampl_repump_cooling_asf,       profile=_URUKUL1_PROFILE_SPINPOLARIZATION)
        self.core.break_realtime()

        self.dds_repump_qubit.set_mu(self.freq_repump_qubit_shutter_ftw,    asf=self.ampl_shutter_asf,              profile=_URUKUL1_PROFILE_DOPPLERCOOLING)
        self.dds_repump_qubit.set_mu(self.freq_repump_qubit_shutter_ftw,    asf=self.ampl_shutter_asf,              profile=_URUKUL1_PROFILE_READOUT)
        self.dds_repump_qubit.set_mu(self.freq_repump_qubit_shutter_ftw,    asf=self.ampl_shutter_asf,              profile=_URUKUL1_PROFILE_RESCUE)
        self.dds_repump_qubit.set_mu(self.freq_repump_qubit_ftw,            asf=self.ampl_repump_qubit_asf,         profile=_URUKUL1_PROFILE_DEPLETION)
        self.dds_repump_qubit.set_mu(self.freq_repump_qubit_shutter_ftw,    asf=self.ampl_shutter_asf,              profile=_URUKUL1_PROFILE_SPINPOLARIZATION)
        self.core.break_realtime()

        # set qubit AOM DDS waveforms
        self.dds_qubit.set_mu(self.freq_qubit_ftw,                          asf=self.ampl_qubit_asf,                profile=_URUKUL0_PROFILE_QUBIT)
        self.dds_qubit.set_mu(self.freq_qubit_shutter_ftw,                  asf=self.ampl_shutter_asf,              profile=_URUKUL0_PROFILE_OFF)
        self.core.break_realtime()


    @rpc(flags={"async"})
    def update_dataset(self, time_mu, pmt_counts):
        """
        Records values via rpc to minimize kernel overhead.
        """
        self.append_to_dataset('rabi_flopping', [self.core.mu_to_seconds(time_mu), pmt_counts])


    def analyze(self):
        """
        Analyze the results from the experiment.
        """
        # tmp remove
        self.pmt_discrimination = 17

        # turn dataset into numpy array for ease of use
        self.rabi_flopping = np.array(self.rabi_flopping)

        # get sorted x-values (time, seconds)
        time_list_s = sorted(set(self.rabi_flopping[:, 0]))

        # collate results
        collated_results = {
            time: []
            for time in time_list_s
        }
        for time_s, pmt_counts in self.rabi_flopping:
            collated_results[time_s].append(pmt_counts)

        # process counts for mean and std and put into processed dataset
        for i, (time_s, count_list) in enumerate(collated_results.items()):
            binned_count_list = np.heaviside(np.array(count_list) - self.pmt_discrimination, 1)
            self.rabi_flopping_processed[i] = np.array([time_s, np.mean(binned_count_list)])
