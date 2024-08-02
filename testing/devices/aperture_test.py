from LAX_exp.base import LAXExperiment
import labrad


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
        self.aperture.open_aperture()       # turn on oven
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def run_main(self):
        delay(15*s)

    # ANALYSIS
    def analyze_experiment(self):

        self.cleanup_devices()

    @rpc
    def cleanup_devices(self):
        # turn off oven
        self.aperture.close_aperture()

