import extensions
import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
import labrad


class OvenTest(LAXExperiment, Experiment):
    """
    Experiment: Oven Test

    Test Oven
    """
    name = 'Oven Test'

    def build_experiment(self):

        self.setattr_device('oven')


    def prepare_experiment(self):
        pass

    @property
    def results_shape(self):
        return (2,2)

    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self):
        self.oven.on()     # turn on oven
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def run_main(self):
        delay(10*s)

    # ANALYSIS
    def analyze_experiment(self):

        self.cleanup_devices()

    @rpc
    def cleanup_devices(self):
        # turn off oven
        self.oven.off()

