
from artiq.experiment import *

import labrad
import numpy as np
from artiq.experiment import *

from LAX_exp.base import LAXExperiment

>>>>>>> 524534b4c277a6647269fd41893ccfcf493fad2b
class OvenTest(LAXExperiment, Experiment):
    """
    todo: document
    """
    name = "Oven Test"

    def build_experiment(self):

        # grab oven
        self.setattr_device('oven')


    def prepare_experiment(self):

        pass

    @property
    def results_shape(self):
        return (2, 2)


    def initialize_experiment(self):
        pass

    def run_main(self):
        self.oven.set_oven_voltage(0.)
        self.oven.set_oven_current(0.)
        self.oven.toggle(False)
        print(self.oven.get_oven_current())
        print(self.oven.get_oven_voltage())

    def analyze_experiment(self):
        pass


