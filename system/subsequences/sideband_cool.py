from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence


class SidebandCool(LAXSubsequence):
    """
    Subsequence: Sideband Cooling

    Cool the ion to the ground state using the S-1/2 D-5/2 transition.
    """
    name = 'sideband_cool'

    parameters = {
        'freq_sideband_cooling_list_ftw':   ('freq_sideband_cooling_list_mhz',          mhz_to_ftw),
        'time_sideband_cooling_list_mu':    ('timing.time_sideband_cooling_list_us',    us_to_mu),

        'time_spinpol_mu':                  ('timing.time_spinpol_us',                  us_to_mu),
        'time_repump_qubit_mu':             ('timing.time_repump_qubit_us',             us_to_mu),
    }
    devices = [
        'probe',
        'pump',
        'qubit_repump',
        'qubit'
    ]

    @kernel(flags={"fast-math"})
    def run(self):
        # set cooling waveforms
        with parallel:
            self.probe.cooling()
            self.qubit.rsb()

        # ensure state preparation is properly interspersed
        for time_list_mu in self.time_sideband_cooling_list_mu:

            # spin polarization/redistribute S-1/2 (397)
            self.probe.cfg_sw(True)
            delay_mu(self.time_redist_mu)
            self.probe.cfg_sw(False)

            # sweep pi-pulse times
            for time_sideband_mu in time_list_mu:
                # todo: change frequencies to hit multiple modes
                # qubit pi-pulse
                self.qubit.cfg_sw(True)
                delay_mu(time_sideband_mu)
                self.qubit.cfg_sw(False)

                # qubit repump
                self.qubit_repump.cfg_sw(True)
                delay_mu(self.time_repump_sideband_cooling_mu)
                self.qubit_repump.cfg_sw(False)

        # repump qubit after sideband cooling
        self.qubit_repump.cfg_sw(True)
        delay_mu(self.time_repump_qubit_mu)
        self.qubit_repump.cfg_sw(False)
