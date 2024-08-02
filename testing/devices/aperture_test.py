from LAX_exp.base import LAXExperiment
import labrad
from artiq.experiment import *

class ApertureTest(LAXExperiment, Experiment):
    """
    Experiment: Aperture Test

    Test Oven
    """
    name = 'Aperture Test'

    def build_experiment(self):

        self.setattr_device('aperture')


    def prepare_experiment(self):
        pass

    @property
    def results_shape(self):
        return (2,2)

    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self):
        self.aperture.open_aperture()       # open aperture
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def run_main(self):
        delay(3*s)

    # ANALYSIS
    def analyze_experiment(self):

        self.cleanup_devices()

    @rpc
    def cleanup_devices(self):
        # close aperture
        self.aperture.close_aperture()

