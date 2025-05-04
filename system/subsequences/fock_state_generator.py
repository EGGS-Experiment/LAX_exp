import numpy as np
from EGGS_labrad.config.dds_config import dds_config
from artiq.build.lib.artiq.language.environment import EnumerationValue
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence
from math import factorial
from scipy.special import genlaguerre
from LAX_exp.extensions.physics_constants import *


class FockStateGenerator(LAXSubsequence):
    """
    Subsequence: Rabi Flop on Sidebands to create arbitrary n state
    Assume we start in |S1/2,mj=-1/2>|n=0> state and always end in |S1/2,mj=-1/2> electronic state
    """

    name = 'fock_state_generator'
    kernel_invariants = {
        "freq_carrier_rabiflop_ftw", "time_carrier_pi_pulse_mu", 'att_readout_mu', 'profile_fock',
        "dds_configs",
        "time_pi_pulses_mu",
        "apply_carrier"
    }

    def build_subsequence(self, profile_fock=2):
        self.profile_fock = profile_fock

        self.setattr_device('qubit')

        # laser frequencies
        self.setattr_argument("freq_carrier_rabiflop_mhz", NumberValue(default=101.1051, step=1, precision=4, min=50., max=400.),
                              group='fock_state_generation')
        self.setattr_argument("freq_rsb_rabiflop_mhz", NumberValue(default=100.7735, step=1, precision=4, min=50., max=400.),
                              group='fock_state_generation')
        self.setattr_argument("freq_bsb_rabiflop_mhz", NumberValue(default=101.4412, step=1, precision=4, min=50., max=400.),
                              group='fock_state_generation')

        # laser parameters
        self.setattr_argument("ampl_qubit_pct", NumberValue(default=50, step=5, precision=1, min=1, max=50),
                              group='fock_state_generation')
        self.setattr_argument("att_flopping_db", NumberValue(default=8, step=0.5, precision=1, min=8, max=31.5),
                              group='fock_state_generation')

        # pulse timings
        self.setattr_argument('time_carrier_pi_pulse_us', NumberValue(default=2.29, step=0.1, precision=3, min=0, max=1000),
                              group='fock_state_generation')
        self.setattr_argument('time_sideband_pi_pulse_us', NumberValue(default=27.1, step=0.1, precision=3, min=0, max=1000),
                              group='fock_state_generation')
        self.setattr_argument('motional_mode', EnumerationValue(["EGGS", "RF", "AXIAL"]))

        # fock state
        self.setattr_argument("final_fock_state", NumberValue(default=10, step=1, precision=0, min=0, max=10),
                                                      group='fock_state_generation')

    def prepare_subsequence(self):
        """
        Prepare and precompute experiment values.
        """
        # convert input arguments to machine units
        freq_rsb_rabiflop_ftw = self.qubit.frequency_to_ftw(self.freq_rsb_rabiflop_mhz * MHz)
        freq_bsb_rabiflop_ftw = self.qubit.frequency_to_ftw(self.freq_bsb_rabiflop_mhz * MHz)
        self.freq_carrier_rabiflop_ftw = self.qubit.frequency_to_ftw(self.freq_carrier_rabiflop_mhz * MHz)
        self.ampl_qubit_asf = self.qubit.amplitude_to_asf(self.ampl_qubit_pct / 100.)
        self.att_readout_mu = att_to_mu(self.att_flopping_db * dB)
        self.time_carrier_pi_pulse_mu = self.core.seconds_to_mu(self.time_carrier_pi_pulse_us*us)

        # calculate lamb dicke parameter to scale successive RSB/BSB pulse times
        if self.motional_mode == 'EGGS' or self.motional_mode == 'RF':
            lambe_dicke_projection = 1/2
        else:
            lambe_dicke_projection = 1/np.sqrt(2)

        omega = 4 * np.pi * (self.freq_carrier_rabiflop_mhz - self.freq_rsb_rabiflop_mhz) * MHz # extra factor of 2 for AOM units
        lamb_dicke = lambe_dicke_projection * (2 * np.pi / 729e-9) * np.sqrt(hbar/(2*mCa*omega))

        ### create config array for each RSB/BSB pulse
        self.time_pi_pulses_us = np.zeros(self.final_fock_state)
        self.dds_configs = np.zeros((self.final_fock_state, 2), dtype = np.int32)

        # configure arrays so we can alternate between rsb and bsb pi pulses
        for idx, fock_state in enumerate(range(self.final_fock_state)):
            ### configure pulses
            self.time_pi_pulses_us[idx] = self.time_sideband_pi_pulse_us / (
                    np.sqrt(1 / (fock_state + 1)) * genlaguerre(fock_state, 1)(lamb_dicke**2)
            )

            if idx % 2 == 0:
                self.dds_configs[idx, :] = [freq_bsb_rabiflop_ftw, self.ampl_qubit_asf]
            else:
                self.dds_configs[idx, :] = [freq_rsb_rabiflop_ftw, self.ampl_qubit_asf]

        self.time_pi_pulses_mu = np.array([self.core.seconds_to_mu(time_pi_pulses_us*us) for time_pi_pulses_us in self.time_pi_pulses_us])

        # determine if we need an additional pi pulse at end to reset to S1/2, mj = -1/2
        self.apply_carrier = (self.final_fock_state % 2 != 0 and self.final_fock_state != 0)

        if self.final_fock_state == 0 :
            self.dds_configs = np.array([[0,0]], dtype = np.int32)
            self.time_pi_pulses_mu = np.array([0], dtype = np.int64)

    @kernel(flags={"fast-math"})
    def run(self) -> TNone:
        # prepare target profile
        self.qubit.set_profile(self.profile_fock)
        self.qubit.cpld.io_update.pulse_mu(8)
        self.qubit.set_att_mu(self.att_readout_mu)
        delay_mu(8)

        # continually apply bsb and rsb (alternating) to achieve the fock state we want
        for idx in range(self.final_fock_state):
            # note: this 8ns delay is CRITICAL to proper operation
            delay_mu(8)
            self.qubit.set_mu(self.dds_configs[idx,0], asf=self.dds_configs[idx,1], profile=self.profile_fock)
            self.qubit.on()
            delay_mu(self.time_pi_pulses_mu[idx])
            self.qubit.off()


        # reset to S1/2, mj=-1/2 if needed via carrier pi-pulse
        if self.apply_carrier:
            self.qubit.set_mu(self.freq_carrier_rabiflop_ftw, asf=self.ampl_qubit_asf, profile=self.profile_fock)
            self.qubit.on()
            delay_mu(self.time_carrier_pi_pulse_mu)
            self.qubit.off()

