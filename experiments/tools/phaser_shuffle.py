import numpy as np
from artiq.experiment import *
from artiq.coredevice.dac34h84 import DAC34H84
from artiq.coredevice.trf372017 import TRF372017

from LAX_exp.system.subsequences import PhaserShuffle


class PhaserShuffle(LAXExperiment):
    """
    Tool: Phaser Shuffle

    todo: document
    """
    name = "Phaser Shuffle"
    kernel_invariants = {
        'idk', 'idk2'
    }

    def build_experiment(self):
        # set core arguments
        self.setattr_argument("repetitions",        NumberValue(default=10, precision=0, step=1, min=1, max=100000))
        self.setattr_argument("randomize_config",   BooleanValue(default=True))

        # frequency configuration
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

        # convert config times from us to mu

        # if randomize_config is enabled, completely randomize the sweep configuration
        if self.randomize_config: np.random.shuffle(self.phaser_config_list)

    def _prepare_argument_checks(self) -> TNone:
        """
        Check experiment arguments for validity.
        """
        # todo: frequencies valid
        # todo: ampls valid
        # todo: times valid
        # todo: ensure we're using an upconverter phaser and phaser_eggs is NOT the same phaser


        # ensure NCO frequency is valid
        if (self.freq_nco_mhz > 400.) or (self.freq_nco_mhz < -400.):
            raise Exception("Invalid phaser NCO frequency. Must be in range [-400, 400].")
        elif (self.freq_nco_mhz > 300.) or (self.freq_nco_mhz < -300.):
            print("Warning: Phaser NCO frequency outside passband of [-300, 300] MHz.")

    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        """
        Initialize and configure hardware elements on the phaser.
        """
        self.core.reset()
        self.shuffle_subsequence.run_config(self.config_list, self.repetitions)
        delay_mu(10000)
