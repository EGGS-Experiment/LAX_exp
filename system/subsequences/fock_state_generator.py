import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence
from math import factorial
from scipy.special import genlaguerre


class FockStateGenerator(LAXSubsequence):
    """
    Subsequence: Rabi Flop on Sidebands to create arbitrary n state
    Assume we start in |S1/2,mj=-1/2>|n=0> state and always end in |S1/2,mj=-1/2> electronic state
    """

    name = 'fock_state_generator'
    kernel_invariants = {
        "freq_carrier_rabiflop_ftw",
        "time_carrier_pi_pulse_mu",
        "dds_configs",
        "time_pi_pulses_mu"
    }

    def build_subsequence(self):
        self.setattr_device('qubit')

        # laser frequencies
        self.setattr_argument("freq_carrier_rabiflop_mhz", NumberValue(default=101.1187, step=1, precision=4, min=50., max=400.),
                              group='fock_state_generation')
        self.setattr_argument("freq_rsb_rabiflop_mhz", NumberValue(default=100.7868, step=1, precision=4, min=50., max=400.),
                              group='fock_state_generation')
        self.setattr_argument("freq_bsb_rabiflop_mhz", NumberValue(default=101.4516, step=1, precision=4, min=50., max=400.),
                              group='fock_state_generation')

        # laser parameters
        self.setattr_argument("ampl_qubit_pct", NumberValue(default=50, step=5, precision=1, min=1, max=50),
                              group='fock_state_generation')
        self.setattr_argument("att_flopping_db", NumberValue(default=8, step=0.5, precision=1, min=8, max=31.5),
                              group='fock_state_generation')

        # pulse timings
        self.setattr_argument('time_carrier_pi_pulse_us',
                                                           NumberValue(default=2, step=0.1, precision=3, min=0, max=1000),
                                                           group='fock_state_generation')
        self.setattr_argument('time_sideband_pi_pulse_us',
                                                                 NumberValue(default=30, step=0.1, precision=3, min=0, max=1000),
                                                                 group='fock_state_generation')

        # fock state
        self.setattr_argument("final_fock_state", NumberValue(default=0, step=1, precision=0, min=0, max=10),
                                                      group='fock_state_generation')

    def prepare_subsequence(self):
        # convert input arguments to machine units
        self.freq_rsb_rabiflop_ftw = self.qubit.frequency_to_ftw(self.freq_rsb_rabiflop_mhz * MHz)
        self.freq_bsb_rabiflop_ftw = self.qubit.frequency_to_ftw(self.freq_bsb_rabiflop_mhz * MHz)
        self.freq_carrier_rabiflop_ftw = self.qubit.frequency_to_ftw(self.freq_carrier_rabiflop_mhz * MHz)
        self.ampl_qubit_asf = self.qubit.amplitude_to_asf(self.ampl_qubit_pct / 100.)
        self.att_readout_mu = att_to_mu(self.att_flopping_db * dB)
        self.time_carrier_pi_pulse_mu = us_to_mu(self.time_carrier_pi_pulse_us)

        ### tmp remove
        lamb_dicke = 0.05
        ### tmp remove

        ### create array of fock state to create
        self.fock_states = np.arange(0, self.final_fock_state + 1)
        self.time_pi_pulses_us = np.zeros(len(self.fock_states))
        self.dds_configs = np.zeros((len(self.fock_states), 2))

        # configure arrays so we can alternate between rsb and bsb pi pulses
        for idx, fock_state in enumerate(self.fock_states):
            ### configure bsb pulses
            if idx % 2 == 0:
                self.time_pi_pulses_us[idx] = self.time_sideband_pi_pulse_us * np.sqrt(
                    factorial(fock_state) / factorial(fock_state + 1)) * genlaguerre(fock_state, 1)(lamb_dicke)

                self.dds_configs[idx, :] = [self.freq_bsb_rabiflop_ftw, self.ampl_qubit_asf]
            ### configure  rsb pulse
            else:
                self.time_pi_pulses_us[idx] = self.time_sideband_pi_pulse_us * np.sqrt(
                    factorial(fock_state - 1) / factorial(fock_state)) * genlaguerre(fock_state - 1, 1)(lamb_dicke)
                self.dds_configs[idx, :] = [self.freq_rsb_rabiflop_ftw, self.ampl_qubit_asf]

        self.time_pi_pulses_mu = np.array([us_to_mu(time_pi_pulses_us) for time_pi_pulses_us in self.time_pi_pulses_us])

        # determine if we need an additional pi pulse at end to reset to S1/2, mj = -1/2
        self.apply_carrier = self.final_fock_state % 2 == 0

    @kernel(flags={"fast-math"})
    def run(self) -> TNone:
        self.core.break_realtime()

        # continually apply bsb and rsb (alternating) to achieve the fock state we want
        for idx in range(self.final_fock_state):
            self.qubit.set_mu(np.int32(self.dds_configs[idx, 0]), asf=np.int32(self.dds_configs[idx, 1]), profile=0)
            self.qubit.on()
            delay_mu(self.time_pi_pulses_mu[idx])
            self.qubit.off()

        # reset to S1/2, mj=-1/2 if needed
        if self.apply_carrier:
            self.qubit.set_mu(np.int32(self.freq_carrier_rabiflop_ftw), asf=np.int32(self.ampl_qubit_asf), profile=0)
            self.qubit.on()
            delay_mu(self.time_carrier_pi_pulse_mu)
            self.qubit.off()
