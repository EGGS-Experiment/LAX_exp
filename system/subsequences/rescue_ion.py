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
        "deathcount_length", "deathcount_tolerance",
        # count processing
        "input_type", "count_threshold", "_enable_thresholding"
    }

    def build_subsequence(self, input_type: TStr='counts'):
        """
        Defines the main interface for the subsequence.
        :param input_type: input type for detect_death. Can be one of ['counts', 'probability'].
            If 'counts', then input is thresholded according to the dataset pmt.count_rate_bright_3ms
            and pmt.count_rate_dark_3ms.
            If 'probability', then input is directly processed (no thresholding).
        """
        # build arguments
        self.input_type = input_type

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

        # calculate state discrimination threshold
        time_readout_us =   self.get_parameter('time_readout_us', group='timing', override=False)
        count_rate_bright = self.get_parameter('count_rate_bright_3ms', group='pmt', override=False) * (time_readout_us / 3000.)
        count_rate_dark =   self.get_parameter('count_rate_dark_3ms', group='pmt', override=False) * (time_readout_us / 3000.)
        self.count_threshold = count_rate_bright / np.log(1 + count_rate_bright / count_rate_dark)

        # process input type
        if not isinstance(self.input_type, str) or (self.input_type not in ('counts', 'probability')):
            raise ValueError("Invalid input type. Must be one of ['counts', 'probability'].")
        elif self.input_type == "counts":
            self._enable_thresholding = True
        elif self.input_type == "probability":
            self._enable_thresholding = False

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
        # threshold incoming counts
        counts_tmp = counts
        if self._enable_thresholding:
            if counts_tmp > self.count_threshold:
                counts_tmp = 0  # bright => P = 0
            else:
                counts_tmp = 1  # dark => P = 1
        # store result and update average
        self._deathcount_arr[self._deathcount_iter % self.deathcount_length] = counts_tmp
        self._deathcount_sum_counts += counts_tmp
        self._deathcount_iter += 1
        delay_mu(15000) # 15us


        '''PROCESS FILTER RESULTS'''
        # only start processing filter once we have filled the array once
        if self._deathcount_iter > self.deathcount_length:
            # subtract history from running average (i.e. circular buffer)
            self._deathcount_sum_counts -= self._deathcount_arr[self._deathcount_iter % self.deathcount_length]

            # process syndromes
            if self._deathcount_sum_counts > (self.deathcount_length - self.deathcount_tolerance):
                self._deathcount_status = 1     # syndrome: ion death (no bright counts)
            elif self._deathcount_sum_counts < self.deathcount_tolerance:
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

