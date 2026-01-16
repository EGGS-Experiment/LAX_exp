from LAX_exp.base import LAXSubsequence
from artiq.experiment import *
from artiq.language import Scannable
from LAX_exp.system.subsequences import InitializeQubit, Readout, RescueIon

class LaserScan(LAXSubsequence):
    """
    Sequence: Laser Scan
    """
    name = 'Laser Scan'

    def build_subsequence(self, **kwargs):

        group_name = 'laser_scan_subsequence'

        self.setattr_argument('pulse_freq_range_mhz',
                              Scannable(deafult=[
                                  ExplicitScan([101.0968]),
                                  CenterScan(101.09, 0.1, 0.0001),
                                  RangeScan(101.09, 101.1, 20)
                              ], global_min=90., global_max=110., precision=5,
                                    global_step =1e-6, unit='MHz'),
                              group=group_name)

        self.setattr_argument('time_pulse_us',
                              NumberValue(default=3, min=0.01, max=1000,
                                          step=0.1, precision = 2, unit='us'),
                              group=group_name)

        self.setattr_argument('ampl_pulse_pct',
                              NumberValue(default=50., min=0, max=50.,
                                          step=0.1, precision = 2, unit='%'),
                                          group=group_name)

        self.setattr_argument('att_pulse_dB', NumberValue(default=8.,
                              min=0., max=31.5, step=0.5, precision=1, unit='dB'),
                            group=group_name)

        self.setattr_argument('enable_pulse_shaping', BooleanValue(default=True))

        ### SUBSEQUENCES ###
        # note: must be initialized last b/c it depends on arguments
        self.pulseshape_subsequence =   QubitPulseShape(
            self, ram_profile=self.profile_729_readout, ram_addr_start=0,
            num_samples=500, ampl_max_pct=self.ampl_qubit_pct,
        )
        self.initialize_subsequence =   InitializeQubit(self)
        self.readout_subsequence =      Readout(self)
        self.rescue_subsequence =       RescueIon(self)

        self.setattr_device('qubit')
        self.setattr_device('probe')
        self.setattr_device('pump')
        self.setattr_device('toptica')

    def prepare_subsequence(self):

    def initialize_subsequence(self):

    def run(self):

    def cleanup_subsequence(self):

