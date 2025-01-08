from artiq.experiment import *

import labrad
import numpy as np
import time
from artiq.experiment import *

from LAX_exp.base import LAXExperiment


class FlipperTest(LAXExperiment, Experiment):
    """
    todo: document
    """
    name = "Flipper Test"

    def build_experiment(self):
        # grab flipper
        self.setattr_device('flipper')

    def prepare_experiment(self):
        pass

    @property
    def results_shape(self):
        return (2, 2)

    @kernel(flags={"fast-math"})
    def initialize_experiment(self):
        self.core.break_realtime()

    def run_main(self):
        for i in range(10):
            self.flipper.flip()
            time.sleep(2)

    def analyze(self):
        pass