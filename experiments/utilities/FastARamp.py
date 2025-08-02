from artiq.experiment import *

from time import sleep
from LAX_exp.language import *


class FastARamp(LAXExperiment, Experiment):
    """
    Utility: Fast A-Ramp

    Eject excess ions in the trap by pulsing the voltage on the A-ramp.
    Does nothing except set a voltage and wait.
    """
    name = 'Fast ARamp'
    kernel_invariants = set()

    def build_experiment(self):
        # ending trap arguments
        self.setattr_argument('normal_aramp_voltage', NumberValue(default=2.8, precision=1, step=0.1, min=0., max=50., scale=1., unit="V"),
                              group='A-Ramp')
        self.setattr_argument("eject_aramp_voltage", NumberValue(default=1, precision=1, step=0.5, min=0., max=50, scale=1., unit='V'),
                              group='A-Ramp')
        self.setattr_argument("time_aramp_s", NumberValue(default=1., precision=1, step=0.1, min=0.1, max=10, scale=1., unit='s'),
                              group='A-Ramp')

        # relevant devices - labrad
        self.setattr_device('trap_dc')

    @property
    def results_shape(self):
        return (1,1)


    '''
    MAIN SEQUENCE
    '''
    def initialize_experiment(self) -> TNone:
        self.trap_dc.set_aramp_voltage(self.normal_aramp_voltage)
        self.trap_dc.aramp_toggle(True)

    def run_main(self) -> TNone:
        # note: add initial "sleep" to ensure everything has settled
        sleep(1.)

        # a-ramp!
        self.trap_dc.set_aramp_voltage(self.eject_aramp_voltage)
        sleep(self.time_aramp_s)
        self.trap_dc.set_aramp_voltage(self.normal_aramp_voltage)
        sleep(self.time_aramp_s)

    def cleanup_experiment(self) -> TNone:
        """
        Restore voltages to normal conditions.
        Note: we keep this here in case of error in run_main.
        """
        self.trap_dc.set_aramp_voltage(self.normal_aramp_voltage)
        self.trap_dc.aramp_toggle(True)
