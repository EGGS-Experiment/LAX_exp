import numpy as np
from artiq.experiment import *

from LAX_exp.base import LAXExperiment
import labrad


class WavemeterTest(LAXExperiment, Experiment):
    """
    Experiment: Wavemeter Test

    Test Wavemeter
    """
    name = 'Wavemeter Test'

    def build_experiment(self):

        self.setattr_device('wavemeter')

    def prepare_experiment(self):
        pass

    @property
    def results_shape(self):
        return (2,2)

    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self):
        pass

    @kernel(flags={"fast-math"})
    def run_main(self):
        delay(4*s)
        print(self.wavemeter.channels["397nm"][0])
        # self.wavemeter.read_channel_frequency(self.wavemeter.channels["397nm"][0])

    # ANALYSIS
    def analyze_experiment(self):
        pass

