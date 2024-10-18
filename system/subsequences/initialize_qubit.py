from artiq.experiment import *
from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence


class InitializeQubit(LAXSubsequence):
    """
    Subsequence: Initialize Qubit

    Initialize the ion in the S-1/2 mj=-1/2 state and cool to the doppler limit.
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
        self.time_spinpol_mu =          self.get_parameter('time_spinpol_us', group='timing', override=True, conversion_function=seconds_to_mu, units=us)
        self.time_repump_qubit_mu =     self.get_parameter('time_repump_qubit_us', group='timing', override=True, conversion_function=seconds_to_mu, units=us)
        self.time_doppler_cooling_mu =  self.get_parameter('time_doppler_cooling_us', group='timing', override=True, conversion_function=seconds_to_mu, units=us)

    @kernel(flags={"fast-math"})
    def run(self) -> TNone:
        # set cooling waveform
        self.pump.cooling()

        # enable cooling repump
        self.repump_cooling.on()

        # repump pulse
        self.repump_qubit.on()
        delay_mu(self.time_repump_qubit_mu)
        self.repump_qubit.off()

        # doppler cooling
        self.pump.on()
        delay_mu(self.time_doppler_cooling_mu)
        self.pump.off()

        # spin polarization
        self.probe.on()
        delay_mu(self.time_spinpol_mu)
        self.probe.off()
