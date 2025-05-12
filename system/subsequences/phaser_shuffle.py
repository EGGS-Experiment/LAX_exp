import numpy as np
from artiq.experiment import *
from LAX_exp.base import LAXSubsequence


class PhaserShuffle(LAXSubsequence):
    """
    Subsequence: Phaser Shuffle

    todo: document
    """
    name = 'phaser_shuffle'
    kernel_invariants = {
        # hardware
        'phaser_dev_name', 'phaser_dev_chan', 'phaser_dev',

        # hardware values
        'freq_band1_center_hz', 'freq_band2_center_hz', 'freq_nco_band1_hz', 'freq_nco_band2_hz',
        'freq_nco_threshold_hz'
    }

    def build_subsequence(self, phaser_att_db: TFloat = 10.):
        # get relevant devices
        self.setattr_device("core")
        self.setattr_device('phaser_eggs')

        self.phaser_dev_name = 'phaser0'
        self.phaser_dev_chan = 0
        self.phaser_att_db = phaser_att_db

    def prepare_subsequence(self):
        """
        Prepare values for speedy evaluation.
        """
        # sanitize inputs
        self._prepare_argument_checks()

        # get target device
        self.phaser_dev = self.get_device(self.phaser_dev_name)

        # set frequencies and bands
        self.freq_nco_threshold_hz =    300 * MHz
        self.freq_band1_center_hz = 175 * MHz
        self.freq_band2_center_hz = 425 * MHz

        self.freq_nco_band1_hz = self.freq_band1_center_hz - 302 * MHz
        self.freq_nco_band2_hz = self.freq_band2_center_hz - 302 * MHz

    def _prepare_argument_checks(self) -> TNone:
        """
        Check experiment arguments for validity.
        """
        if isinstance(self.phaser_att_db, float):
            if not (0. <= self.phaser_att_db <= 31.5):
                raise ValueError("Invalid value for phaser_att_db. Must be in [0, 31.5].")
        else:
            raise ValueError("Invalid type for phaser_att_db. Must be a float.")

        # todo: check that we're using correct phaser somehow


    """
    KERNEL FUNCTIONS
    """
    @kernel(flags={"fast-math"})
    def initialize_subsequence(self) -> TNone:
        """
        Prepare the subsequence immediately before run.
        """
        self.device_cleanup()

    @kernel(flags={"fast-math"})
    def cleanup_subsequence(self) -> TNone:
        """
        Clean up the subsequence immediately after run.
        """
        self.device_cleanup()

    @kernel(flags={"fast-math"})
    def run_config(self, config_list: TArray(TFloat, 2), repetitions: TInt32 = 1) -> TNone:
        """
        todo: document
        """
        # todo: verify inputs?
        # set attenuator
        self.phaser_dev.channel[self.phaser_dev_chan].set_att(self.phaser_att_db * dB)
        delay_mu(10000)

        # create holder variable
        freq_nco_center_hz = 0.

        # run pulses for given number of reps
        for i in range(repetitions):
            delay_mu(25000)

            for config_vals in config_list:
                # unpack config list
                freq_hz =   config_vals[0]
                ampl_pct =  config_vals[1]
                time_mu =   np.int64(config_vals[2])
                if freq_hz < self.freq_nco_threshold_hz:
                    freq_nco_center_hz = self.freq_nco_band1_hz
                    freq_center_hz = self.freq_band1_center_hz
                else:
                    freq_nco_center_hz = self.freq_nco_band2_hz
                    freq_center_hz = self.freq_band2_center_hz
                delay_mu(50000)

                # set NCO and DUC freqs
                at_mu(self.phaser_dev.get_next_frame_mu())
                self.phaser_dev.channel[self.phaser_dev_chan].set_nco_frequency(freq_nco_center_hz)
                at_mu(self.phaser_dev.get_next_frame_mu())
                self.phaser_dev.dac_sync()

                at_mu(self.phaser_dev.get_next_frame_mu())
                self.phaser_dev.channel[self.phaser_dev_chan].set_duc_frequency(freq_hz - freq_center_hz)
                at_mu(self.phaser_dev.get_next_frame_mu())
                self.phaser_dev.duc_stb()

                # prepare oscillators
                at_mu(self.phaser_dev.get_next_frame_mu())
                self.phaser_dev.channel[self.phaser_dev_chan].oscillator[0].set_frequency(0.01 * MHz)
                delay_mu(40)
                self.phaser_dev.channel[self.phaser_dev_chan].oscillator[0].set_amplitude_phase(ampl_pct, 0.)

                # run pulse
                delay_mu(time_mu)
                self.phaser_dev.channel[self.phaser_dev_chan].oscillator[0].set_amplitude_phase(0., 0., clr=1)

        # ensure outputs are disabled
        self.device_cleanup()

    @kernel(flags={"fast-math"})
    def device_cleanup(self) -> TNone:
        """
        Ensure all device outputs are disabled.
        """
        self.phaser_dev.channel[self.phaser_dev_chan].set_att_mu(0x00)
        delay_mu(10000)

        for i in range(5):
            self.phaser_dev.channel[self.phaser_dev_chan].oscillator[i].set_amplitude_phase(0., 0., clr=1)
            delay_mu(8000)

    @kernel(flags={"fast-math"})
    def run(self) -> TNone:
        pass

