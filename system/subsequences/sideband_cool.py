from artiq.experiment import *

import numpy as np
from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence


class SidebandCool(LAXSubsequence):
    """
    Subsequence: Sideband Cool

    Cool the ion to the ground state using RSB pulses on the S-1/2 to D-5/2 transition.
    """
    name = 'sideband_cool'

    def build_subsequence(self):
        # get devices
        self.setattr_device('probe')
        self.setattr_device('pump')
        self.setattr_device('repump_qubit')
        self.setattr_device('qubit')

        # sideband cooling configuration
        self.setattr_argument('calibration',                            BooleanValue(default=False), group='sideband_cooling')
        self.setattr_argument('sideband_cycles',                        NumberValue(default=40, ndecimals=0, step=1, min=1, max=10000), group='sideband_cooling')
        self.setattr_argument('cycles_per_spin_polarization',           NumberValue(default=11, ndecimals=0, step=1, min=1, max=10000), group='sideband_cooling')

        # sideband cooling timing
        self.setattr_argument('time_min_sideband_cooling_us_list',      PYONValue([50, 75, 80, 91]), group='sideband_cooling')
        self.setattr_argument('time_max_sideband_cooling_us_list',      PYONValue([250, 271, 239, 241]), group='sideband_cooling')

        # sideband cooling waveform
        self.setattr_argument('freq_sideband_cooling_mhz_list',         PYONValue([104.012, 103.012, 105.012, 107.711]), group='sideband_cooling')
        self.setattr_argument('ampl_sideband_cooling_pct',              NumberValue(default=50, ndecimals=5, step=1, min=10, max=100), group='sideband_cooling')

    def prepare_subsequence(self):
        # ensure input has correct dimensions and uses < 7 modes (due to max of 8 profiles per urukul channel)
        num_min_times = len(list(self.time_min_sideband_cooling_us_list))
        num_max_times = len(list(self.time_max_sideband_cooling_us_list))
        num_modes =     len(list(self.freq_sideband_cooling_mhz_list))
        assert num_modes < 7, "Exceeded maximum number of cooling frequencies."
        assert num_min_times == num_max_times == num_modes, "Number of modes and timings are not equal."

        # timing parameters
        self.time_repump_qubit_mu =         self.get_parameter('time_repump_qubit_us', group='timing', override=True, conversion_function=seconds_to_mu, units=us)
        self.time_spinpol_mu =              self.get_parameter('time_spinpol_us', group='timing', override=True, conversion_function=seconds_to_mu, units=us)

        # calibration setup
        self.qubit_func =                                               self.qubit.off if self.calibration is True else self.qubit.on

        # calculate number of spin polarizations
        # todo: do number of spin polarizations more accurately
        num_spin_polarizations =                                        int(self.sideband_cycles / self.cycles_per_spin_polarization + 1)

        # sequence sideband cooling pulses for each mode
        self.time_sideband_cooling_list_mu =                            np.array(
                                                                            self.core.seconds_to_mu(
                                                                                np.linspace(
                                                                                    self.time_min_sideband_cooling_us_list,
                                                                                    self.time_max_sideband_cooling_us_list,
                                                                                    self.sideband_cycles
                                                                                ) * us
                                                                            )
                                                                        )
        self.time_sideband_cooling_list_mu =                            np.array_split(self.time_sideband_cooling_list_mu, num_spin_polarizations)

        # sideband cooling waveforms
        self.freq_sideband_cooling_ftw_list =                           np.array([hz_to_ftw(freq_mhz * MHz) for freq_mhz in self.freq_sideband_cooling_mhz_list])
        self.ampl_sideband_cooling_asf =                                pct_to_asf(self.ampl_sideband_cooling_pct)
        self.iter_sideband_cooling_modes_list =                         np.array(range(1, 1 + len(self.freq_sideband_cooling_ftw_list)))

    @kernel(flags={"fast-math"})
    def initialize_subsequence(self):
        # set sideband cooling profiles for 729nm qubit laser
        # profile 0: reserved for readout
        # profile 1 & greater: sideband cooling
        for i in self.iter_sideband_cooling_modes_list:
            self.qubit.set_mu(self.freq_sideband_cooling_ftw_list[i - 1], asf=self.ampl_sideband_cooling_asf, profile=i)
            self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def run(self):
        # run sideband cooling with interspersed spin polarization
        for time_list_mu in self.time_sideband_cooling_list_mu:

            # spin polarization/state preparation (397 probe beam)
            self.probe.on()
            delay_mu(self.time_spinpol_mu)
            self.probe.off()

            # sweep pi-pulse times
            for time_modes_mu in time_list_mu:

                # sweep over modes
                for i in self.iter_sideband_cooling_modes_list:

                    # set the waveform for the i-th mode
                    self.qubit.set_profile(i)

                    # qubit pi-pulse
                    self.qubit_func()               # we use qubit_func() instead of self.qubit.on()
                                                    # to allow for variable behavior due to calibration
                    delay_mu(time_modes_mu[i - 1])
                    self.qubit.off()

                    # qubit repump
                    self.repump_qubit.on()
                    delay_mu(self.time_repump_qubit_mu)
                    self.repump_qubit.off()
