from artiq.experiment import *
from LAX_exp.LAX.base_classes import LAXSubsequence, us_to_mu, mhz_to_ftw


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
        'time_doppler_cooling_mu':          ('timing.time_doppler_cooling_us',          us_to_mu),
        'time_repump_sideband_cooling_mu':  ('timing.time_repump_sideband_cooling_us',  us_to_mu)
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
        self.probe.cooling()
        self.qubit.rsb()

        # ensure state preparation is properly interspersed
        for time_list_mu in self.time_sideband_cooling_list_mu:

            # spin polarization/redistribute S-1/2 (397)
            self.probe.on()
            delay_mu(self.time_redist_mu)
            self.probe.off()

            # sweep pi-pulse times
            for time_sideband_mu in time_list_mu:
                # todo: change frequencies
                # qubit pi-pulse
                self.qubit.on()
                delay_mu(time_sideband_mu)
                self.qubit.off()

                # qubit repump
                self.qubit_repump.on()
                delay_mu(self.time_repump_sideband_cooling_mu)
                self.qubit_repump.off()

        # repump qubit after sideband cooling
        self.qubit_repump.on()
        delay_mu(self.time_repump_qubit_mu)
        self.qubit_repump.off()
