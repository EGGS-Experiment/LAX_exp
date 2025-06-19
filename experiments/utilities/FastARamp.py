import numpy as np
from artiq.experiment import *

import time
from LAX_exp.language import *


class FastARamp(LAXExperiment, Experiment):
    """
    Utility: Fast A-Ramp

    Eject excess ions in the trap by pulsing the voltage on the A-ramp.
    Does nothing except run the
    """
    name = 'Fast ARamp'

    kernel_invariants = {
        "pmt_sample_num", "pmt_dark_threshold_counts",
        "att_397_mu", "att_866_mu", "att_854_mu",
        "IMAGE_HEIGHT", "IMAGE_WIDTH", "image_region", "data_path", "time_aramp_pulse_s"
    }

    def build_experiment(self):
        # ending trap arguments
        self.setattr_argument("repetitions",            NumberValue(default=5, precision=0, step=1, min=1, max=10000))
        self.setattr_argument('normal_aramp_voltage',   NumberValue(default=2.8, precision=1, step=0.1, min=0., max=50.), group='A-Ramp')
        self.setattr_argument("time_aramp_s",           NumberValue(default=1, precision=0, step=1, min=1, max=100000, unit='s', scale=1.), group='A-Ramp')
        self.setattr_argument("aramp_voltage_list",     Scannable(
                                                                default=[
                                                                    ExplicitScan([14]),
                                                                    RangeScan(18, 24, 20, randomize=True),
                                                                ],
                                                                global_min=0.0, global_max=30.0, global_step=1,
                                                                unit="V", scale=1, precision=2
                                                            ), group='A-Ramp')

        # relevant devices - labrad
        self.setattr_device('trap_dc')

    def prepare_experiment(self):
        """
        Prepare experimental values and precompute/preallocate to reduce kernel overheads.
        """
        self.time_aramp_pulse_s = 2.


    @property
    def results_shape(self):
        return (1,1)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self):
        """
        Hardware initialization immediately before run().
        """
        # note: we add sleep and set voltages AGAIN to reflect updates on GUIs
        time.sleep(3)
        self.trap_dc.set_east_endcap_voltage(self.east_endcap_voltage)
        self.trap_dc.aramp_toggle(True)
        self.core.break_realtime()

    def run_main(self) -> TNone:
        """
        Run till ion is loaded or timeout.
        """
        for i in range(self.repetitions):
            self.aramp_ions()
        self.cleanup_devices()


    @rpc
    def aramp_ions(self) -> TNone:
        """
        Pulse A-ramp to get rid of excess ions.
        """
        # todo: allow for some error conditions
        for aramp_voltage in self.aramp_voltage_list:

            # pulse A-ramp for given period
            print(f"ARAMPING AT VOLTAGE {np.round(aramp_voltage,2)}")
            self.trap_dc.set_aramp_voltage(aramp_voltage)
            time.sleep(self.time_aramp_pulse_s)
            self.trap_dc.set_aramp_voltage(self.final_aramp_voltage)
            time.sleep(self.time_aramp_pulse_s)

            # check termination
            self.check_termination()

        self.aperture.close_aperture()

    @rpc
    def cleanup_devices(self) -> TNone:
        """
        Restore A-Ramp voltage to normal conditions.
        """
        self.trap_dc.set_aramp_voltage(self.normal_aramp_voltage)
        self.trap_dc.aramp_toggle(True)

    @rpc
    def check_termination(self) -> TNone:
        """
        OVERRIDE base check_termination to ensure devices are always cleaned up
        in case of termination.
        """
        if self.scheduler.check_termination():
            self.cleanup_devices()
            raise TerminationRequested
