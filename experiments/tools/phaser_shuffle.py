import numpy as np
from artiq.experiment import *
from LAX_exp.system.subsequences import PhaserShuffle


class PhaserShuffle(LAXExperiment):
    """
    Tool: Phaser Shuffle

    todo: document
    """
    name = "Phaser Shuffle"
    kernel_invariants = set()

    def build_experiment(self):
        # set core arguments
        self.setattr_argument("repetitions",        NumberValue(default=10, precision=0, step=1, min=1, max=100000))
        self.setattr_argument("randomize_config",   BooleanValue(default=True))

        # frequency configuration
        self.setattr_argument("phaser_att_db",      NumberValue(default=3, precision=1, step=0.5, min=0., max=31.5))
        self.setattr_argument("phaser_config_list", ([[65., 100., 1000], [72., 100., 1000], [81., 100., 1000],
                                                      [83., 100., 1000], [85., 100., 1000], [91., 100., 1000],
                                                      [210., 100., 1000], [376., 100., 1000], [521., 100., 1000]]),
                              tooltip="[[freq_mhz, ampl_pct, time_us]]")

        # instantiate subsequences
        self.shuffle_subsequence = PhaserShuffle()

    def prepare_experiment(self):
        """
        Prepare kernel values before running.
        """
        # convert configs to numpy array
        # note: if config list dimensions are not all same, then we get an error automatically
        self.phaser_config_list = np.array(self.phaser_config_list)
        self._prepare_argument_checks()

        # convert config units (MHz to Hz, pct to frac, us to mu)
        self.phaser_config_list[:, 0] *= MHz
        self.phaser_config_list[:, 1] /= 100.
        self.phaser_config_list[:, 2] = self.core.seconds_to_mu(self.phaser_config_list[:, 2] * us)

        # if randomize_config is enabled, completely randomize the sweep configuration
        if self.randomize_config: np.random.shuffle(self.phaser_config_list)

    def _prepare_argument_checks(self) -> TNone:
        """
        Check experiment arguments for validity.
        """
        # check that all elements in phaser_config_list are valid
        if not all(50. < freq_mhz < 550. for freq_mhz in self.phaser_config_list[:, 0]):
            raise ValueError("Invalid frequencie(s) in phaser_config_list. Must be in [50., 550.] MHz.")
        if not all(0.01 < ampl_pct < 100. for ampl_pct in self.phaser_config_list[:, 1]):
            raise ValueError("Invalid amplitude(s) in phaser_config_list. Must be in [0.01, 100) %.")
        if not all(1 < time_us < 10000000 for time_us in self.phaser_config_list[:, 2]):
            raise ValueError("Invalid time(s) in phaser_config_list. Must be in [1us, 10s].")

        # todo: ensure we're using an upconverter phaser and phaser_eggs is NOT the same phaser

    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        self.shuffle_subsequence.run_config(self.phaser_config_list, self.repetitions)
        delay_mu(10000)

