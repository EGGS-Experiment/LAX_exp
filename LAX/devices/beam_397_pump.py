from artiq.experiment import *
from LAX_exp.LAX.base_classes import LAXDevice, mhz_to_ftw, pct_to_asf


class Beam397Pump(LAXDevice):
    """
    Wrapper for the 397nm pump beam.
        Uses the DDS channel to drive an AOM in double-pass configuration.
    """
    name = "pump"

    device_names = {'beam': 'urukul1_ch1'}
    device_parameters = {
        'freq_cooling_ftw': ('beams.freq_mhz.freq_pump_cooling_mhz', mhz_to_ftw),
        'ampl_cooling_asf': ('beams.ampl_pct.ampl_pump_cooling_pct', pct_to_asf),
        'freq_readout_ftw': ('beams.freq_mhz.freq_pump_readout_mhz', mhz_to_ftw),
        'ampl_readout_asf': ('beams.ampl_pct.ampl_pump_readout_pct', pct_to_asf)
    }

    @kernel(flags='fast-math')
    def prepare_devices(self):
        # set cooling and readout profiles
        self.core.break_realtime()
        self.set_mu(self.freq_cooling_ftw, asf=self.ampl_cooling_asf, profile=0)
        self.core.break_realtime()
        self.set_mu(self.freq_readout_ftw, asf=self.ampl_readout_asf, profile=1)
