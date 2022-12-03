from artiq.experiment import *
from LAX_exp.LAX.base_classes import LAXDevice, mhz_to_ftw, pct_to_asf


class Beam397Probe(LAXDevice):
    """
    Wrapper for the 397nm probe beam (polarized).
        Uses the DDS channel to drive an AOM in double-pass configuration.
    """
    name = "probe"

    device_names = {'beam': 'urukul1_ch0'}
    device_parameters = {
        'freq_redist_ftw': ('beams.freq_mhz.freq_probe_redist_mhz', mhz_to_ftw),
        'ampl_redist_asf': ('beams.ampl_pct.ampl_probe_redist_pct', pct_to_asf)
    }


    @kernel(flags='fast-math')
    def prepare_devices(self):
        # set cooling and readout profiles
        self.core.break_realtime()
        self.set_mu(self.freq_redist_ftw, asf=self.ampl_redist_asf, profile=0)
        self.core.break_realtime()
        self.set_mu(self.freq_redist_ftw, asf=self.ampl_redist_asf, profile=1)
