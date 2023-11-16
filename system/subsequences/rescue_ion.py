from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence


class RescueIon(LAXSubsequence):
    """
    Subsequence: Rescue Ion

    Rescue the ion by running the 397nm cooling beam at rescue parameters.
    """
    name = 'rescue_ion'

    def build_subsequence(self):
        # get devices
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')
        self.setattr_device('repump_qubit')

        # rescue cooling (happens at end of each repetition)
        self.setattr_argument("rescue_enable",                          BooleanValue(default=False), group=self.name)
        self.setattr_argument("repetitions_per_rescue",                 NumberValue(default=1, ndecimals=0, step=1, min=1, max=10000), group=self.name)
        self.setattr_argument("resuscitate_ion",                        BooleanValue(default=False), group=self.name)
        # self.setattr_argument("iterations_per_resuscitate",             NumberValue(default=1, ndecimals=0, step=1, min=1, max=10000), group='rescue_ion')
        # self.setattr_argument("additional_cooling_time_s",              NumberValue(default=1, ndecimals=5, step=0.1, min=0, max=10000), group='rescue_ion')

    def prepare_subsequence(self):
        self.time_rescue_mu =               self.get_parameter('time_rescue_us', group='timing', override=True, conversion_function=seconds_to_mu, units=us)
        self.time_resuscitate_mu =          self.get_parameter('time_resuscitate_us', group='timing', override=True, conversion_function=seconds_to_mu, units=us)

        # configure variable behavior for self.resuscitate
        self.resuscitate = self._no_op
        if self.resuscitate_ion is True:
            self.resuscitate = self._resuscitate

    @kernel(flags={"fast-math"})
    def run(self, i):
        # check whether rescuing is enabled
        if self.rescue_enable == True:

            # check whether it's time to rescue the ion
            if (i > 0) and (i % self.repetitions_per_rescue == 0):

                # set rescue waveform and ensure 866 is on
                self.pump.rescue()
                self.repump_cooling.on()
                self.repump_qubit.on()

                # doppler cooling
                self.pump.on()
                delay_mu(self.time_rescue_mu)
                self.pump.off()

    @kernel(flags={"fast-math"})
    def _resuscitate(self):
        # set rescue waveform and ensure 866 is on
        self.pump.rescue()
        self.repump_cooling.on()
        self.repump_qubit.on()

        # doppler cooling
        self.pump.on()
        delay_mu(self.time_resuscitate_mu)
        self.pump.off()

    @kernel(flags={"fast-math"})
    def _no_op(self):
        delay_mu(10)
