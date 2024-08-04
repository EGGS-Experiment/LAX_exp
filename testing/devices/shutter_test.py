import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
import labrad


class ShutterTest(LAXExperiment, Experiment):
    """
    Experiment: Shutter Test

    Shutter Oven
    """
    name = 'Shutter Test'

    def build_experiment(self):

        self.setattr_device('shutters')

    def prepare_experiment(self):
        pass
    @property
    def results_shape(self):
        return (2,2)

    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self):
        # open shutters
        self.core.break_realtime()
        delay(3*s)
        self.shutters.open_377_shutter()
        self.core.break_realtime()
        self.shutters.open_423_shutter()
        self.core.break_realtime()
        delay(10*s)


    @kernel(flags={"fast-math"})
    def run_main(self):
        delay(3*s)

    # ANALYSIS
    def analyze_experiment(self):

        self.cleanup_devices()

    @rpc
    def cleanup_devices(self):
        # turn off oven
        self.shutters.close_377_shutter()
        self.shutters.close_423_shutter()
        self.shutters.close_labjack()

