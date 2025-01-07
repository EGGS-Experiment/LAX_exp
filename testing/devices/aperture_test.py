from artiq.experiment import *

import labrad
import numpy as np
from artiq.experiment import *

from LAX_exp.base import LAXExperiment
class ApertureTest(LAXExperiment, Experiment):
    """
    todo: document
    """
    name = "Aperture Test"

    def build_experiment(self):

        # grab oven
        self.setattr_device('aperture')
        self.aperture_wait_time = 1 * s


    def prepare_experiment(self):
        pass

    @property
    def results_shape(self):
        return (2, 2)

    @kernel(flags={"fast-math"})
    def initialize_experiment(self):
        self.core.break_realtime()

    def run_main(self):
        self.open_aperture()

    def analyze(self):
        pass

    def pulse_aperture(self):
        self.aperture.open()
        self.aperture.wait(self.aperture_wait_time)
        self.aperture.close()

    def open_aperture(self):
        self.aperture.open()
        self.aperture.wait(self.aperture_wait_time)

    def close_aperture(self):
        self.aperture.close()
        self.aperture.wait(self.aperture_wait_time)

