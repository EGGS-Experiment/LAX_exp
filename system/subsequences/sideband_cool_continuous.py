from artiq.experiment import *

import numpy as np
from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence


class SidebandCoolContinuous(LAXSubsequence):
    """
    Subsequence: Sideband Cool - Continuous

    Cool the ion to the ground state using a continuous RSB pulse on the S-1/2 to D-5/2 transition.
    """
    name = 'sideband_cool_continuous'

    def build_subsequence(self):
        # get devices
        self.setattr_device('probe')
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')
        self.setattr_device('repump_qubit')
        self.setattr_device('qubit')

        # sideband cooling configuration
        self.setattr_argument("calibration_continuous",                 BooleanValue(default=False), group='sideband_cooling.continuous')
        self.setattr_argument("sideband_cycles_continuous",             NumberValue(default=1, ndecimals=0, step=1, min=1, max=10000), group='sideband_cooling.continuous')
        self.setattr_argument("time_sideband_cooling_us",               NumberValue(default=16000, ndecimals=3, step=100, min=0.001, max=1000000), group='sideband_cooling.continuous')
        self.setattr_argument("pct_per_spin_polarization",              NumberValue(default=20, ndecimals=3, step=1, min=0.01, max=100), group='sideband_cooling.continuous')

        # sideband cooling modes
        self.setattr_argument("freq_sideband_cooling_mhz_pct_list",     PYONValue({102.7765: 100}), group='sideband_cooling.continuous')

        # sideband cooling powers
        self.setattr_argument("att_sidebandcooling_continuous_db",      NumberValue(default=8, ndecimals=1, step=0.5, min=8, max=31.5), group='sideband_cooling.continuous')
        self.setattr_argument("ampl_quench_pct",                        NumberValue(default=10, ndecimals=2, step=1, min=5, max=50),    group='sideband_cooling.continuous')

    def prepare_subsequence(self):
        # ensure mode percentages add up to 100%
        mode_total_pct =                                                np.sum(self.freq_sideband_cooling_mhz_pct_list.values())
        # assert mode_total_pct == 100,                                   "Error: total sideband cooling mode percentages exceed 100%."

        # DDS parameters
        self.ampl_qubit_asf =                                           self.get_parameter('ampl_qubit_pct',
                                                                                           group='beams.ampl_pct', override=True,
                                                                                           conversion_function=pct_to_asf)

        # get timing parameters
        self.time_repump_qubit_mu =                                     self.get_parameter('time_repump_qubit_us',
                                                                                           group='timing', override=True,
                                                                                           conversion_function=seconds_to_mu, units=us)
        self.time_spinpol_mu =                                          self.get_parameter('time_spinpol_us',
                                                                                           group='timing', override=True,
                                                                                           conversion_function=seconds_to_mu, units=us)
        # get waveform parameters
        self.freq_repump_qubit_ftw =                                    self.get_parameter('freq_repump_qubit_mhz',
                                                                                           group='beams.freq_mhz', override=False,
                                                                                           conversion_function=hz_to_ftw, units=MHz)


        ### SIDEBAND COOLING ###

        # CONFIG
        self.qubit_func =                                               self.qubit.off if self.calibration_continuous is True else self.qubit.on
        self.freq_sideband_cooling_ftw_list =                           np.array([hz_to_ftw(freq_mhz * MHz)
                                                                                  for freq_mhz in self.freq_sideband_cooling_mhz_pct_list.keys()])
        self.iter_sideband_cooling_modes_list =                         np.array(list(range(1, 1 + len(self.freq_sideband_cooling_mhz_pct_list))))

        # POWER
        self.ampl_quench_asf =                                          pct_to_asf(self.ampl_quench_pct)
        self.att_sidebandcooling_mu =                                   att_to_mu(self.att_sidebandcooling_continuous_db * dB)

        # TIMING
        self.time_sideband_cooling_mu =                                 self.core.seconds_to_mu(self.time_sideband_cooling_us * us)
        # create list of cycle times
        cycle_time_us =                                                 self.time_sideband_cooling_us / self.sideband_cycles_continuous
        cycle_timings_us_list =                                         np.linspace(0, self.time_sideband_cooling_us, self.sideband_cycles_continuous + 1)[:-1]

        # create list of SBC profile times within a cycle
        cycle_profile_times_pct =                                       np.array([0] + list(self.freq_sideband_cooling_mhz_pct_list.values()))
        cycle_profile_timings_us =                                      np.cumsum(cycle_profile_times_pct / 100 * cycle_time_us)[:-1]

        # create final timing list
        self.delay_sideband_cooling_cycle_us_list =                     np.array([time_cycle_start_us + cycle_profile_timings_us
                                                                                 for time_cycle_start_us in cycle_timings_us_list])
        self.delay_sideband_cooling_cycle_mu_list =                     np.array([self.core.seconds_to_mu(time_us_list * us)
                                                                                  for time_us_list in self.delay_sideband_cooling_cycle_us_list])

        # spin polarization timings
        self.time_spinpolarization_mu_list =                            self.core.seconds_to_mu(np.arange(0, 1, self.pct_per_spin_polarization / 100)
                                                                                                * (self.time_sideband_cooling_us * us))
        self.time_spinpolarization_mu_list =                            self.time_spinpolarization_mu_list[1:]


    @kernel(flags={"fast-math"})
    def initialize_subsequence(self):
        # set quench waveform
        self.repump_qubit.set_mu(self.freq_repump_qubit_ftw, asf=self.ampl_quench_asf, profile=3)

        # set sideband cooling profiles for 729nm qubit laser
            # profile 0: reserved for readout
            # profile 1 & greater: sideband cooling
        for i in self.iter_sideband_cooling_modes_list:
            self.qubit.set_mu(self.freq_sideband_cooling_ftw_list[i - 1], asf=self.ampl_qubit_asf, profile=i)
            self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def run(self):
        # prepare beams for sideband cooling
        with parallel:
            # set sideband cooling profiles for regular beams
            self.repump_qubit.set_profile(3)

            # set sideband cooling attenuation for qubit beam
            self.qubit.set_att_mu(self.att_sidebandcooling_mu)

        # do spin polarization before SBC per Guggemos' thesis
        self.probe.on()
        self.repump_cooling.on()
        delay_mu(self.time_spinpol_mu)
        self.probe.off()
        self.repump_cooling.off()

        # turn on sideband cooling beams
        with parallel:
            # turn on qubit repump
            self.repump_qubit.on()

            # turn on qubit beam
            self.qubit_func()           # we use qubit_func() instead of self.qubit.on()
                                        # to allow for variable behavior due to calibration

        # intersperse state preparation with normal SBC
        time_start_mu = now_mu()
        with parallel:

            # interleave sideband cooling of different modes every cycle
            for time_delay_mu_list in self.delay_sideband_cooling_cycle_mu_list:
                # for i in self.iter_sideband_cooling_modes_list:
                for i in range(len(time_delay_mu_list)):
                    # switch profile at set time
                    at_mu(time_start_mu + time_delay_mu_list[i])
                    self.qubit.set_profile(self.iter_sideband_cooling_modes_list[i])

            # spin polarization
            for time_spinpol_mu in self.time_spinpolarization_mu_list:
                # do spin polarization at set time
                at_mu(time_start_mu + time_spinpol_mu)
                self.probe.on()
                self.repump_cooling.on()
                delay_mu(self.time_spinpol_mu)
                self.probe.off()
                self.repump_cooling.off()

        # stop sideband cooling
        at_mu(time_start_mu + self.time_sideband_cooling_mu)
        with parallel:
            self.repump_qubit.set_profile(1)
            self.qubit.off()

        # repump qubit after sideband cooling
        delay_mu(self.time_repump_qubit_mu)
        self.repump_qubit.off()
