from LAX_exp.base import LAXExperiment
import labrad
from artiq.experiment import *

class TrapDCTest(LAXExperiment, Experiment):
    """
    Experiment: Trap DC Test

    Test Oven
    """
    name = 'Trap DC Test'

    def build_experiment(self):

        self.setattr_device('trap_dc')


    def prepare_experiment(self):
        pass

    @property
    def results_shape(self):
        return (2,2)

    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self):
        self.trap_dc.set_east_endcap_voltage(19.)
        self.trap_dc.set_west_endcap_voltage(25.)
        self.trap_dc.east_endcap_on()
        self.trap_dc.west_endcap_on()
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def run_main(self):
        delay(3*s)
        self.trap_dc.set_east_endcap_voltage(202.)
        self.trap_dc.set_west_endcap_voltage(287.)
        self.core.break_realtime()
        delay(3*s)

    # ANALYSIS
    def analyze_experiment(self):

        self.cleanup_devices()

    @rpc
    def cleanup_devices(self):
        # reset trap voltages
        self.trap_dc.set_east_endcap_voltage(19.)
        self.trap_dc.set_west_endcap_voltage(25.)
        self.trap_dc.east_endcap_off()
        self.trap_dc.west_endcap_off()

