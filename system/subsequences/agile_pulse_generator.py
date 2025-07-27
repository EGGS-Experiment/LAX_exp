import numpy as np
from artiq.experiment import *
from artiq.coredevice.ad9910 import PHASE_MODE_CONTINUOUS

from LAX_exp.language import *


class AgilePulseGenerator(LAXSubsequence):
    """
    Subsequence: Agile Pulse Generator
    todo: document
    """
    name = 'agile_pulse_generator'
    kernel_invariants = {
        "profile_agile", "pulse_config", "att_pulse_db",
        "att_pulse_mu", "num_pulses", "dds_configs", "pulse_times_mu"
    }

    def build_subsequence(self, profile_agile: TInt32 = 2, att_pulse_db: TFloat = 31.5,
                          pulse_config: TArray(TFloat, 2) = np.zeros(1, 3)):
        """
        Defines the main interface for the subsequence.
        Arguments:
            profile_agile: the AD9910 RAM profile to use.
            att_pulse_db: the DDS attenuation to set during the pulses.
            pulse_config: the pulse configuration - an array of [freq_mhz, ampl_pct, time_us].
        """
        # set subsequence parameters
        self.profile_agile =    profile_agile
        self.att_pulse_db =     att_pulse_db
        self.pulse_config =     pulse_config

        # get relevant devices
        self.setattr_device("qubit")

    def prepare_subsequence(self):
        """
        Prepare and precompute experiment values for speedy evaluation.
        """
        '''VALIDATE INPUTS'''
        self._prepare_argument_checks()

        '''CONVERT PULSE CONFIGURATION'''
        self.att_pulse_mu = att_to_mu(self.att_pulse_db)

        # note: create separate arrays for freq/ampl & time to avoid int32 conversion later on
        # (b/c all variables have to have same type)
        self.num_pulses = len(self.pulse_config)
        self.dds_configs = [
            [self.qubit.frequency_to_ftw(config[0] * MHz),
             self.qubit.amplitude_to_asf(config[1] / 100.)]
            for config in self.pulse_config
        ]
        self.pulse_times_mu = [
            self.core.seconds_to_mu(config[2] * us)
            for config in self.pulse_config
        ]

    def _prepare_argument_checks(self) -> TNone:
        """
        Check experiment arguments for validity.
        """
        # note: validate inputs here to get around bugs where args are passed from setattr_argument
        if self.profile_agile not in range(0, 7):
            raise ValueError("Invalid AD9910 profile for agile_pulse_generator: {:d}. Must be in [0, 7].".format(self.profile_agile))

        # parse pulse_config
        pulse_freqs_mhz = self.pulse_config[:, 0]
        pulse_ampls_pct = self.pulse_config[:, 1]
        pulse_times_us = self.pulse_config[:, 2]
        if not (all(pulse_freqs_mhz <= 400.) and all(pulse_freqs_mhz >= 30.)):
            raise ValueError("Invalid pulse configuration frequencies: {:d}. Must be in [30, 400] MHz.".format(pulse_freqs_mhz))
        elif not (all(pulse_ampls_pct <= 50.) and all(pulse_ampls_pct >= 0.01)):
            raise ValueError("Invalid pulse configuration amplitudes: {:d}. Must be in [0.01, 50] MHz.".format(pulse_ampls_pct))
        elif not (all(pulse_times_us <= 100000.) and all(pulse_times_us >= 0.5)):
            raise ValueError("Invalid pulse configuration times: {:d}. Must be in [0.5, 100000] MHz.".format(pulse_times_us))


    '''
    MAIN SEQUENCE
    '''
    @kernel(flags={"fast-math"})
    def run(self) -> TNone:
        """
        Run pulse configuration on qubit DDS.
        """
        # prepare DDS
        self.qubit.set_profile(self.profile_agile)
        self.qubit.cpld.io_update.pulse_mu(8)
        self.qubit.set_att_mu(self.att_pulse_mu)

        # run pulses
        for i in range(self.num_pulses):
            self.qubit.set_mu(self.dds_configs[i, 0], asf=self.dds_configs[i, 1], profile=self.profile_agile,
                              phase_mode=PHASE_MODE_CONTINUOUS)
            self.qubit.on()
            delay_mu(self.pulse_times_mu[i])
            self.qubit.off()

