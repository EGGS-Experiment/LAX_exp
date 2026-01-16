import numpy as np
from artiq.experiment import *

from LAX_exp.language import *
from LAX_exp.system.subsequences import (
    InitializeQubit, Readout, QubitPulseShape, RescueIon, NoOperation,
    SidebandCoolContinuousRAM, FockStateGenerator
)

from sipyco import pyon



class AutomatedDissociaton(LAXExperiment, Experiment):
    """
    Experiment: Automated Dissociation
    Measures ion fluorescence vs 729nm pulse time and frequency.
    """
    name = 'Automated Dissociation'
    kernel_invariants = {
    }

    def build_experiment(self):
        pass

    def prepare_experiment(self):
        pass

    @property
    def results_shape(self):
        return (self.repetitions * len(self.time_rabiflop_mu_list),
                2)

    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self) -> TNone:
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        self.core.break_realtime()

    '''
    ANALYSIS
    '''

    def analyze_experiment(self):
       pass