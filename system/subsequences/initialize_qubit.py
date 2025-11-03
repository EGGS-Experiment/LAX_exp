from artiq.experiment import *
from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence


class InitializeQubit(LAXSubsequence):
    """
    Subsequence: Initialize Qubit

    Initialize the ion in the S-1/2 mj=-1/2 state, then cool to the doppler limit.
    """
    name = 'initialize_qubit'
    kernel_invariants = {
        "time_spinpol_mu",
        "time_repump_qubit_mu",
        "time_doppler_cooling_mu"
    }

    def build_subsequence(self):
        self.setattr_device('probe')
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')
        self.setattr_device('repump_qubit')

    def prepare_subsequence(self):
        self.time_spinpol_mu =          self.get_parameter('time_spinpol_us', group='timing', override=True,
                                                           conversion_function=seconds_to_mu, units=us)
        self.time_repump_qubit_mu =     self.get_parameter('time_repump_qubit_us', group='timing', override=True,
                                                           conversion_function=seconds_to_mu, units=us)
        self.time_doppler_cooling_mu =  self.get_parameter('time_doppler_cooling_us', group='timing', override=True,
                                                           conversion_function=seconds_to_mu, units=us)

    @kernel(flags={"fast-math"})
    def run(self) -> TNone:
        """
        Quench ion from D-5/2, doppler, then run spinpol.
        Typical stateprep sequence sans SBC.
        """
        # ensure 397nm spinpol beam is off before starting
        self.probe.off()

        # set cooling waveform
        self.pump.cooling()

        # enable cooling repump (866nm)
        # 2025/03/27: evidently we don't turn 866nm off - should we?
        self.repump_cooling.on()

        # repump pulse
        self.repump_qubit.on()
        delay_mu(self.time_repump_qubit_mu)
        # 2025/03/27: why bother turning off BEFORE doppler?
        # should leave on instead b/c 397 has 393 component, which can scatter into D-5/2
        # self.repump_qubit.off()

        # doppler cooling
        self.pump.on()
        delay_mu(self.time_doppler_cooling_mu)
        self.pump.off()

        # spin polarization
        self.probe.on()
        delay_mu(self.time_spinpol_mu)
        self.probe.off()

        # 2025/03/27: ensure 854nm off
        self.repump_qubit.off()

    @kernel(flags={"fast-math"})
    def quench(self) -> TNone:
        """
        Quick quench function to simplify quenching for other exps.
        Note: 866nm is also on.
        """
        # set target profile
        self.pump.readout()

        # run quench
        # note: no 8ns delays here b/c done judiciously
        self.repump_qubit.on()
        self.repump_cooling.on()
        delay_mu(self.time_repump_qubit_mu)
        self.repump_cooling.off()
        self.repump_qubit.off()

    @kernel(flags={"fast-math"})
    def spin_polarize(self) -> TNone:
        """
        Quick spinpol function to simplify spinpol for other exps.
        """
        # set target profile
        self.pump.readout()

        # run spinpol
        # note: yes 8ns delays here b/c spinpol & 866 have no ext sw
        self.probe.on()
        delay_mu(8)
        self.repump_cooling.on()
        delay_mu(self.time_spinpol_mu)
        self.probe.off()
        delay_mu(8)
        self.repump_cooling.off()
        delay_mu(8)
