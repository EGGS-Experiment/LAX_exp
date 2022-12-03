from artiq.experiment import *
from LAX_exp.LAX.base_classes import LAXDevice, mhz_to_ftw, pct_to_asf


class Beam729(LAXDevice):
    """
    Wrapper for the 729nm qubit beam (polarized).
        Uses the DDS channel to drive an AOM in double-pass configuration.
    """
    name = "qubit"

    device_names = {'beam': 'urukul0_ch1'}
    device_parameters = {
        'freq_qubit_carrier_ftw': ('beams.freq_mhz.freq_qubit_carrier_mhz', mhz_to_ftw),
        'freq_qubit_rsb_ftw': ('beams.freq_mhz.freq_qubit_rsb_mhz', mhz_to_ftw),
        'freq_qubit_bsb_ftw': ('beams.freq_mhz.freq_qubit_bsb_mhz', mhz_to_ftw),
        'ampl_qubit_pct': ('beams.ampl_pct.ampl_qubit_pct', pct_to_asf)
    }


    @kernel(flags='fast-math')
    def prepare_devices(self):
        # set qubit profiles
        self.core.break_realtime()
        self.dev.set_mu(self.freq_qubit_carrier_ftw, asf=self.ampl_qubit_asf, profile=0)
        self.core.break_realtime()
        self.dev.set_mu(self.freq_qubit_rsb_ftw, asf=self.ampl_qubit_asf, profile=1)
        self.core.break_realtime()
        self.dev.set_mu(self.freq_qubit_bsb_ftw, asf=self.ampl_qubit_asf, profile=2)
