from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence


class SidebandCool(LAXSubsequence):
    """
    Subsequence: Sideband Cool

    Cool the ion to the ground state using RSB pulses on the S-1/2 to D-5/2 transition.
    """
    name = 'sideband_cool'

    parameters = {
        'freq_sideband_cooling_list_ftw':   ('freq_sideband_cooling_list_mhz',          mhz_to_ftw),
        'time_sideband_cooling_list_mu':    ('timing.time_sideband_cooling_list_us',    us_to_mu),

        'time_spinpol_mu':                  ('timing.time_spinpol_us',                  us_to_mu),
        'time_repump_qubit_mu':             ('timing.time_repump_qubit_us',             us_to_mu),
    }

    def build_subsequence(self):
        self.setattr_device('probe')
        self.setattr_device('pump')
        self.setattr_device('repump_qubit')
        self.setattr_device('qubit')


    @kernel(flags={"fast-math"})
    def run(self):
        # set cooling waveforms
        self.probe.cooling()

        # ensure state preparation is properly interspersed
        for time_list_mu in self.time_sideband_cooling_list_mu:

            # spin polarization/redistribute S-1/2 (397)
            self.probe.on()
            delay_mu(self.time_redist_mu)
            self.probe.off()

            # sweep pi-pulse times
            for time_sideband_mu in time_list_mu:
                # todo: change frequencies to hit multiple modes
                # qubit pi-pulse
                self.qubit.on()
                delay_mu(time_sideband_mu)
                self.qubit.off()

                # qubit repump
                self.repump_qubit.on()
                delay_mu(self.time_repump_sideband_cooling_mu)
                self.repump_qubit.off()

        # repump qubit after sideband cooling
        self.repump_qubit.on()
        delay_mu(self.time_repump_qubit_mu)
        self.repump_qubit.off()
