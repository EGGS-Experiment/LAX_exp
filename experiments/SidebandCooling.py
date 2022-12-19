import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import InitializeQubit, RabiFlop, Readout


class SidebandCooling2(LAXExperiment, Experiment):
    """
    Sideband Cooling 2
    Measures temperature after a given number of RSB pulses.
    """

    name = 'Rabi Flopping 2'

    def build_experiment(self):
        self.setattr_argument("calibration",                                BooleanValue(default=False))
        self.setattr_argument("sideband_cycles",                            NumberValue(default=100, ndecimals=0, step=1, min=1, max=10000))
        self.setattr_argument("cycles_per_spin_polarization",               NumberValue(default=20, ndecimals=0, step=1, min=1, max=10000))

        # sideband cooling
        self.setattr_argument("time_repump_sideband_cooling_us",            NumberValue(default=20, ndecimals=5, step=1, min=1, max=1000000))
        self.setattr_argument("time_min_sideband_cooling_us_list",          PYONValue([50, 75, 80, 91]))
        self.setattr_argument("time_max_sideband_cooling_us_list",          PYONValue([250, 271, 239, 241]))
        self.setattr_argument("freq_sideband_cooling_mhz_list",             PYONValue([104.012, 103.012, 105.012, 107.711]))
        self.setattr_argument("ampl_sideband_cooling_pct",                  NumberValue(default=50, ndecimals=5, step=1, min=10, max=100))

        # readout
        self.setattr_argument("freq_rsb_scan_mhz",                          Scannable(
                                                                                default=CenterScan(104.012, 0.04, 0.001),
                                                                                global_min=30, global_max=200, global_step=1,
                                                                                unit="MHz", scale=1, ndecimals=5
                                                                            ))

        self.setattr_argument("freq_bsb_scan_mhz",                          Scannable(
                                                                                default=CenterScan(105.214, 0.04, 0.001),
                                                                                global_min=30, global_max=200, global_step=1,
                                                                                unit="MHz", scale=1, ndecimals=5
                                                                            ))

        self.setattr_argument("time_readout_pipulse_us",                    NumberValue(default=250, ndecimals=5, step=1, min=1, max=10000))
        #self.setattr_argument("ampl_readout_pipulse_pct",                  NumberValue(default=50, ndecimals=5, step=1, min=1, max=100))

    def prepare_experiment(self):
        # combine and shuffle readout frequencies
        self.freq_qubit_scan_ftw =                                  [mhz_to_ftw(freq_mhz * MHz)
                                                                     for freq_mhz in (list(self.freq_rsb_scan_mhz) + list(self.freq_bsb_scan_mhz))]
        shuffle(self.freq_qubit_scan_ftw)

        # calculate number of spin polarizations
        num_spin_depolarizations = int(self.sideband_cycles / self.cycles_per_spin_polarization)
        if (num_spin_depolarizations < 1):
            num_spin_depolarizations = 1

        # sideband cooling timing
        self.time_sideband_cooling_list_mu =                        np.array([
                                                                        self.core.seconds_to_mu(time_us * us)
                                                                        for time_us in np.linspace(
                                                                            self.time_min_sideband_cooling_us_list,
                                                                            self.time_max_sideband_cooling_us_list,
                                                                            self.sideband_cycles
                                                                        )
                                                                    ])
        self.time_sideband_cooling_list_mu =                        np.array_split(self.time_sideband_cooling_list_mu, num_spin_depolarizations)

        # other sideband cooling parameters
        self.time_repump_sideband_cooling_mu =                      us_to_mu(self.time_repump_sideband_cooling_us)
        self.freq_sideband_cooling_ftw_list =                       [self.dds_qubit.frequency_to_ftw(freq_mhz * MHz) for freq_mhz in self.freq_sideband_cooling_mhz_list]
        self.ampl_sideband_cooling_asf =                            self.dds_qubit.amplitude_to_asf(self.ampl_sideband_cooling_pct / 100)
        self.iter_sideband_cooling_modes_list =                     list(range(1, 1 + len(self.freq_sideband_cooling_ftw_list)))

        # readout pi-pulse
        self.time_readout_pipulse_mu =                              self.core.seconds_to_mu(self.time_readout_pipulse_us * us)
        self.ampl_qubit_asf =                                       self.dds_qubit.amplitude_to_asf(self.ampl_qubit_pct / 100)
        #self.ampl_readout_pipulse_asf =                            self.dds_qubit.amplitude_to_asf(self.ampl_readout_pipulse_pct / 100)

        # calibration setup
        self.calibration_qubit_status =                             not self.calibration

        # get devices
        self.setattr_device('qubit')
        self.setattr_device('pmt')

        # prepare sequences
        self.initialize_subsequence =                               InitializeQubit(self)
        self.readout_subsequence =                                  Readout(self)

        # dataset
        self.set_dataset('results',                                 np.zeros((self.repetitions * len(self.time_rabiflop_mu_list), 2)))
        self.setattr_dataset('results')


    # PREPARE MAIN SEQUENCE
    @kernel
    def run_initialize(self):
        self.core.reset()

        # set qubit beam parameters
        self.qubit.set_mu(self.freq_rabiflop_mhz, asf=self.qubit.ampl_qubit_asf)

        # record subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.core.break_realtime()

        self.readout_subsequence.record_dma()
        self.core.break_realtime()

        # load subsequences from DMA
        self.initialize_subsequence.load_dma()
        self.core.break_realtime()

        self.readout_subsequence.load_dma()
        self.core.break_realtime()


    # MAIN LOOP
    @kernel
    def run_main(self):
        for trial_num in range(self.repetitions):
            # sweep time
            for time_rabi_pair_mu in self.time_rabiflop_mu_list:

                # initialize ion in S-1/2 state
                self.initialize_subsequence.run_dma()

                # wait given time
                delay_mu(time_rabi_pair_mu[0])

                # rabi flopping w/qubit laser
                self.qubit.cfg_sw(True)
                delay_mu(time_rabi_pair_mu[1])
                self.qubit.cfg_sw(False)

                # do readout
                self.readout_subsequence.run_dma()

                # update dataset
                with parallel:
                    self.update_dataset(time_rabi_pair_mu[1], self.pmt_counter.fetch_count())
                    self.core.break_realtime()
