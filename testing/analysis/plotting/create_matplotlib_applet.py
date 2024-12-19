from artiq.experiment import *

import numpy as np
from artiq.experiment import *
from artiq.experiment import EnvExperiment
from artiq.experiment import Experiment
from artiq.gui import applets


class MatplotlibTest(EnvExperiment, Experiment):
    """
    Matplotlib Experiment
    """
    name = "Matplotlib Test"

    def build(self):
        self.setattr_device("ccb")

    def prepare(self):
        pass

    def initialize(self):
        pass

    def run(self):

       pass
    def analyze(self):
        # self.ccb.issue("disable_applet", "first_matplotlib")

        ccb_command = '$python -m LAX_exp.applets.plot_matplotlib temp.x' \
                       ' temp.y --subplot-x-labels temp.x_labels --subplot-y-labels temp.y_labels' \
                       ' --subplot-titles temp.titles --title "Test Plot" --x-label "x labels" ' \
                       '--y-label "y labels"'

        ccb_command += ' --num-subplots 2'

        self.ccb.issue("create_applet", "first_matplotlib1",
                       ccb_command)
        # # #
        # self.set_dataset('temp.x', [[1, 2, 3], [6, 7, 8]], broadcast=True)
        # self.set_dataset('temp.y', [[3, 6, 5], [4, 5, 6]], broadcast=True)
        # self.set_dataset('temp.titles', ["Test Title 1", "Test Title 2"], broadcast=True)
        # self.set_dataset('temp.x_labels', "Label", broadcast=True)
        # self.set_dataset('temp.y_labels', ["Test Y Label 1", "Test Y Label 2"], broadcast=True)
