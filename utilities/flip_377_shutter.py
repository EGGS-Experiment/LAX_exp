import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
import labrad


class Flip377Shutter(LAXExperiment, Experiment):
    """
    Utility to flip 377 shutter
    """

    name = "Flip 377 Shutter"

    def build_experiment(self):
        self.setattr_argument("open_377_shutter", BooleanValue(False))

        self.setattr_argument("close_377_shutter", BooleanValue(True))

        self.setattr_device('shutters')

    def prepare_experiment(self):
        pass

    @property
    def results_shape(self):
        return (2, 2)


    @kernel(flags={"fast-math"})
    def initialize_experiment(self):
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def run_main(self):
        self.core.break_realtime()
        if self.open_377_shutter:
            self.shutters.open_377_shutter()
        elif self.close_377_shutter:
            self.shutters.close_377_shutter()
        self.core.wait_until_mu(now_mu())
        self.core.break_realtime()

        self.shutters.close_labjack()


    # ANALYSIS
    def analyze_experiment(self):
        pass
