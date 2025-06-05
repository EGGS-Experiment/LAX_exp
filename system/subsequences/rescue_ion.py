import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence

_DEATHCOUNT_STATUS_MESSAGES = [
    "CLEAR",                #0
    "ERROR: DEATH",         #1
    "ERROR: TRANSITION"     #2
]


class RescueIon(LAXSubsequence):
    """
    Subsequence: Rescue Ion

    Gallantly attempt to stave off ion death.
    Processes ion symptom data in real-time and responds dynamically with rescue light.
    Rescue the ion by running the 397nm cooling beam at rescue parameters.
    """
    name = 'rescue_ion'
    kernel_invariants = {
        # timing
        "time_rescue_mu", "time_resuscitate_mu", "time_aperture_pulse_s",
        # death counting
        "deathcount_length", "deathcount_tolerance", "deathcount_threshold"
    }

    def build_subsequence(self):
        # get relevant devices
        self.setattr_device('pump')
        self.setattr_device('probe')
        self.setattr_device('repump_cooling')
        self.setattr_device('repump_qubit')
        self.setattr_device('aperture')

        # user configurable arguments to adjust behavior
        self.setattr_argument("enable_rescue",          BooleanValue(default=False), group=self.name,
                              tooltip="todo: document")
        self.setattr_argument("enable_resuscitate",     BooleanValue(default=False), group=self.name,
                              tooltip="todo: document")
        self.setattr_argument("enable_aperture",        BooleanValue(default=False), group=self.name,
                              tooltip="Allow aperture to be opened upon death conditions.")
        self.setattr_argument("add_397nm_spinpol",      BooleanValue(default=False), group=self.name,
                              tooltip="Use 397nm spin polarization light in addition to normal 397nm light.")
        self.setattr_argument("repetitions_per_rescue", NumberValue(default=1, precision=0, step=1, min=1, max=10000),
                              group=self.name, tooltip="Number of experiment repetitions before running a rescue.")

        # magic numbers
        self.time_aperture_pulse_s = 2.0

    def prepare_subsequence(self):
        """
        Prepare values for the subsequence.
        """
        '''PARAMETERS'''
        # parameters - sequence timing
        self.time_rescue_mu =       self.get_parameter('time_rescue_us', group='timing',
                                                       override=True, conversion_function=seconds_to_mu, units=us)
        self.time_resuscitate_mu =  self.get_parameter('time_resuscitate_us', group='timing',
                                                       override=True, conversion_function=seconds_to_mu, units=us)

        # parameters - death detection
        self.deathcount_length =    self.get_parameter('deathcount_length', group='management.death', override=False)
        self.deathcount_tolerance = self.get_parameter('deathcount_tolerance', group='management.death', override=False)
        self.deathcount_threshold = self.get_parameter('deathcount_threshold', group='management.death', override=False)

        '''KERNEL DATA STRUCTURES'''
        # holder arrays
        self._deathcount_arr =              np.zeros(self.deathcount_length, dtype=np.int32)
        self._deathcount_iter =             0
        self._deathcount_sum_counts =       0
        # status flags
        self._deathcount_status =           0
        self._deathcount_status_latched =   0

    @kernel(flags={"fast-math"})
    def initialize_subsequence(self) -> TNone:
        # clear ion status here to prevent experiments in the pipeline from overriding it
        self.set_dataset('management.dynamic.ion_status', 'CLEAR', broadcast=True)

    @kernel(flags={"fast-math"})
    def run(self, i: TInt32) -> TNone:
        # check whether it's time to rescue the ion
        if (self.enable_rescue is True) and (i > 0) and ((i % self.repetitions_per_rescue) == 0):
            # set rescue waveform and ensure repumps are on
            self.pump.rescue()
            self.repump_cooling.on()
            self.repump_qubit.on()

            # rescue cooling
            self.pump.on()
            if self.add_397nm_spinpol is True:
                self.probe.on()
            delay_mu(self.time_rescue_mu)
            self.pump.off()
            self.probe.off()

    @kernel(flags={"fast-math"})
    def resuscitate(self) -> TNone:
        """
        Resuscitate ion: a fast pre-exposure prophylaxis against ion death.
        Intended for regular execution.
        """
        if self.enable_resuscitate:
            # set rescue waveform and ensure repumps are on
            self.pump.rescue()
            self.repump_cooling.on()
            self.repump_qubit.on()

            # resuscitate cooling
            self.pump.on()
            if self.add_397nm_spinpol is True:
                self.probe.on()
            delay_mu(self.time_resuscitate_mu)
            self.pump.off()
            self.probe.off()

    @kernel(flags={"fast-math"})
    def detect_death(self, counts: TInt32) -> TNone:
        """
        Process counts in real-time for low-overhead, dynamic response to certain error conditions.
        :param counts: the most recent ion counts from a standard readout pulse.
        """
        '''UPDATE FILTER'''
        # threshold incoming counts and store data in array
        if counts > self.deathcount_threshold:
            self._deathcount_arr[self._deathcount_iter % self.deathcount_length] = 1
            self._deathcount_sum_counts += 1
        else:
            self._deathcount_arr[self._deathcount_iter % self.deathcount_length] = 0
        delay_mu(15000) # 15us


        '''PROCESS FILTER RESULTS'''
        # start removing values from filter once we have stored enough counts
        if self._deathcount_iter >= self.deathcount_length:
            # subtract history from running average (i.e. circular buffer)
            self._deathcount_sum_counts -= self._deathcount_arr[(self._deathcount_iter + 1) % self.deathcount_length]

            # process syndromes
            if self._deathcount_sum_counts < self.deathcount_tolerance:
                self._deathcount_status = 1     # syndrome: ion death (no bright counts)
            elif self._deathcount_sum_counts > (self.deathcount_length - self.deathcount_tolerance):
                self._deathcount_status = 2     # syndrome: bad transition (no dark counts)
            else:
                self._deathcount_status = 0     # syndrome: clear (no error)
            delay_mu(25000) # 25us

            # process status only if it's changed
            # note: this is important to prevent e.g. multiple updates if symptoms are persistent
            if self._deathcount_status != self._deathcount_status_latched:
                # set syndrome message
                self.set_dataset('management.dynamic.ion_status',
                                 _DEATHCOUNT_STATUS_MESSAGES[self._deathcount_status],
                                 broadcast=True)

                # clear holder variables upon error condition
                if self._deathcount_status != 0:
                    self._deathcount_sum_counts =   0
                    self._deathcount_iter =         0
                    for i in range(self.deathcount_length):
                        self._deathcount_arr[i] =   0

                    # pulse aperture only upon transition into death (not for continuous death)
                    if self.enable_aperture and (self._deathcount_status == 1):
                        self.core.wait_until_mu(now_mu())
                        print('Ion death detected - opening aperture in response.')
                        self.aperture.pulse_aperture_open(self.time_aperture_pulse_s)
                        self.core.break_realtime()
                        delay_mu(500000) # 500us

                # update status flag
                self._deathcount_status_latched = self._deathcount_status
            delay_mu(50000) # 50us

        # update count array iterator
        self._deathcount_iter += 1

