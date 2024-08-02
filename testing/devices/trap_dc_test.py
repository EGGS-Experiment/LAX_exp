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
        self.trap_dc.set_h_shim_voltage(0.)
        self.trap_dc.set_v_shim_voltage(0.)
        self.trap_dc.set_a_ramp2_voltage(0.)
        self.trap_dc.east_endcap_on()
        self.trap_dc.west_endcap_on()
        self.trap_dc.h_shim_off()
        self.trap_dc.v_shim_off()
        self.trap_dc.a_ramp2_off()
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def run_main(self):

        delay(5*s)
        # self.core.wait_until_mu(now_mu())
        # self.trap_dc.set_east_endcap_voltage(202.)
        # self.trap_dc.set_west_endcap_voltage(287.)
        self.trap_dc.ramp_both_endcaps([202, 287], [100,100])
        self.trap_dc.set_h_shim_voltage(51.5)
        self.trap_dc.set_v_shim_voltage(60.)
        self.trap_dc.set_a_ramp2_voltage(2.)
        # self.trap_dc.east_endcap_on()
        # self.trap_dc.west_endcap_on()
        self.trap_dc.h_shim_on()
        self.trap_dc.v_shim_on()
        self.trap_dc.a_ramp2_on()
        self.core.break_realtime()
        delay(5*s)
        self.core.wait_until_mu(now_mu())

    # ANALYSIS
    def analyze_experiment(self):

        self.cleanup_devices()

    @rpc
    def cleanup_devices(self):
        self.trap_dc.set_east_endcap_voltage(19.)
        self.trap_dc.set_west_endcap_voltage(25.)
        self.trap_dc.set_v_shim_voltage(0.)
        self.trap_dc.set_h_shim_voltage(0.)
        self.trap_dc.set_a_ramp2_voltage(0.)
        self.trap_dc.east_endcap_off()
        self.trap_dc.west_endcap_off()
        self.trap_dc.v_shim_off()
        self.trap_dc.h_shim_off()
        self.trap_dc.a_ramp2_off()


