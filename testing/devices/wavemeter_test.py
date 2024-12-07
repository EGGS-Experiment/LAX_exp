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


    def run_main(self):
        print(self.wavemeter.channels["397nm"][0])
        print(self.wavemeter.read_channel_frequency(self.wavemeter.channels["397nm"][0]))
        self.wavemeter.set_channel_frequency(self.wavemeter.channels["423nm"][0],709.077647)
        print(self.wavemeter.read_channel_frequency(self.wavemeter.channels["423nm"][0]))

    # ANALYSIS
    def analyze_experiment(self):
        pass

