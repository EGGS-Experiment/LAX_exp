from artiq.experiment import *

import numpy as np
from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence


class SidebandCoolPulsed(LAXSubsequence):
    """
    Subsequence: Sideband Cool - Pulsed

    Cool the ion to the ground state using RSB pulses on the S-1/2 to D-5/2 transition.
    """
    name = 'sideband_cool_pulsed'
    kernel_invariants = {
        "ampl_qubit_asf",
        "time_repump_qubit_mu",
        "time_spinpol_mu",
        "time_sideband_cooling_list_mu",
        "freq_sideband_cooling_ftw_list",
        "att_sidebandcooling_mu",
        "iter_sideband_cooling_modes_list"
    }

    def build_subsequence(self):
        # get devices
        self.setattr_device('probe')
        self.setattr_device('pump')
        self.setattr_device('repump_qubit')
        self.setattr_device('qubit')

        # sideband cooling configuration
        self.setattr_argument('calibration_pulsed',                     BooleanValue(default=False), group='sideband_cooling.pulsed')
        self.setattr_argument('sideband_cycles_pulsed',                 NumberValue(default=80, precision=0, step=1, min=1, max=10000), group='sideband_cooling.pulsed')
        self.setattr_argument("extra_sideband_cycles",                  NumberValue(default=0, precision=0, step=1, min=0, max=10000), group='sideband_cooling.pulsed')
        self.setattr_argument('cycles_per_spin_polarization',           NumberValue(default=15, precision=0, step=1, min=1, max=10000), group='sideband_cooling.pulsed')

        # sideband cooling timing
        self.setattr_argument("time_form_sideband_cooling",             EnumerationValue(['Linear', 'Inverse Square Root'], default='Linear'), group='sideband_cooling.pulsed')
        self.setattr_argument('time_min_sideband_cooling_us_list',      PYONValue([30]), group='sideband_cooling.pulsed')
        self.setattr_argument('time_max_sideband_cooling_us_list',      PYONValue([150]), group='sideband_cooling.pulsed')

        # sideband cooling waveform
        self.setattr_argument('freq_sideband_cooling_mhz_list',         PYONValue([103.77]), group='sideband_cooling.pulsed')
        self.setattr_argument("att_sidebandcooling_pulsed_db",          NumberValue(default=8, precision=1, step=0.5, min=8, max=31.5), group='sideband_cooling.pulsed')

    def prepare_subsequence(self):
        # ensure input has correct dimensions and uses < 7 modes (due to max of 8 profiles per urukul channel)
        num_min_times =                                                 len(list(self.time_min_sideband_cooling_us_list))
        num_max_times =                                                 len(list(self.time_max_sideband_cooling_us_list))
        num_modes =                                                     len(list(self.freq_sideband_cooling_mhz_list))
        assert num_modes < 7,                                           "Exceeded maximum number of cooling frequencies."
        assert num_min_times == num_max_times == num_modes,             "Number of modes and timings are not equal."

        # DDS parameters
        self.ampl_qubit_asf =                                           self.get_parameter('ampl_qubit_pct',
                                                                                           group='beams.ampl_pct', override=True,
                                                                                           conversion_function=pct_to_asf)

        # timing parameters
        self.time_repump_qubit_mu =                                     self.get_parameter('time_repump_qubit_us',
                                                                                           group='timing', override=True,
                                                                                           conversion_function=seconds_to_mu, units=us)
        self.time_spinpol_mu =                                          self.get_parameter('time_spinpol_us',
                                                                                           group='timing', override=True,
                                                                                           conversion_function=seconds_to_mu, units=us)

        # calibration setup
        self.qubit_func =                                               self.qubit.off if self.calibration_pulsed is True else self.qubit.on


        # timing
        self.time_sideband_cooling_list_mu =                            np.array([])

        # calculate cooling timeform: linear
        if self.time_form_sideband_cooling == 'Linear':

            self.time_sideband_cooling_list_mu = np.array([
                self.core.seconds_to_mu(time_us * us)
                for time_us in np.linspace(
                    self.time_min_sideband_cooling_us_list,
                    self.time_max_sideband_cooling_us_list,
                    self.sideband_cycles_pulsed
                )
            ])

        # calculate cooling timeform: inverse square root
        elif self.time_form_sideband_cooling == 'Inverse Square Root':

            # alias variables for compactness of notation
            steps = self.sideband_cycles_pulsed
            (t_min, t_max) = (self.time_min_sideband_cooling_us_list, self.time_max_sideband_cooling_us_list)

            # calculate timeshaping
            timeshape_t0 = np.sqrt((steps - 1) / (np.power(t_min, -2.) - np.power(t_max, -2.)))
            timeshape_n0 = np.power(timeshape_t0 / t_max, 2.) + (steps - 1)

            # calculate timeshape
            self.time_sideband_cooling_list_mu = timeshape_t0 / np.sqrt(timeshape_n0 - np.array([np.arange(steps)] * len(t_min)).transpose())
            self.time_sideband_cooling_list_mu = np.array([
                self.core.seconds_to_mu(time_mode_list_us * us)
                for time_mode_list_us in self.time_sideband_cooling_list_mu
            ])

        # account for errors in timing
        else:
            raise Exception('Unknown error in SBC - Pulsed')


        # extra sideband cooling cycles
        extra_cycles_arr =                                              np.tile(self.time_sideband_cooling_list_mu[1], (self.extra_sideband_cycles, 1))
        self.time_sideband_cooling_list_mu =                            np.concatenate([extra_cycles_arr, self.time_sideband_cooling_list_mu])

        # split up sideband cooling times to intersperse spin polarization
        num_spin_polarizations =                                        int(self.sideband_cycles_pulsed / self.cycles_per_spin_polarization + 0.5)
        self.time_sideband_cooling_list_mu =                            np.array_split(self.time_sideband_cooling_list_mu, num_spin_polarizations)

        # sideband cooling waveforms
        self.freq_sideband_cooling_ftw_list =                           np.array([hz_to_ftw(freq_mhz * MHz) for freq_mhz in self.freq_sideband_cooling_mhz_list])
        self.att_sidebandcooling_mu =                                   att_to_mu(self.att_sidebandcooling_pulsed_db * dB)
        self.iter_sideband_cooling_modes_list =                         np.array(range(1, 1 + len(self.freq_sideband_cooling_ftw_list)))

    @kernel(flags={"fast-math"})
    def initialize_subsequence(self) -> TNone:
        # set sideband cooling profiles for 729nm qubit laser
        # profile 0: reserved for readout
        # profile 1 & greater: sideband cooling
        for i in self.iter_sideband_cooling_modes_list:
            self.qubit.set_mu(self.freq_sideband_cooling_ftw_list[i - 1], asf=self.ampl_qubit_asf, profile=i)
            self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def run(self) -> TNone:
        # set sideband cooling attenuation
        self.qubit.set_att_mu(self.att_sidebandcooling_mu)

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

            # repump qubit after sideband cooling
            self.repump_qubit.on()
            delay_mu(self.time_repump_qubit_mu)
            self.repump_qubit.off()
