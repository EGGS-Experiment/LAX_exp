import numpy as np
import typing

import artiq.coredevice.ad9910
import artiq.coredevice.ttl
import artiq.coredevice.edge_counter

from dax.experiment import *
from dax.interfaces.detection import DetectionInterface


class DetectionModule(DaxModule, DetectionInterface):
    # System dataset keys
    DETECTION_TIME_KEY = 'detection_time'
    DETECTION_FREQ_KEY = 'detection_freq'
    STATE_DETECTION_THRESHOLD_KEY = 'state_detection_threshold'
    ACTIVE_CHANNELS_KEY = 'active_channels'

    # PMT channels
    MAX_NUM_CHANNELS = 4
    _PMT_TTL_KEYS = [f'ttl{i}' for i in range(MAX_NUM_CHANNELS)]
    _PMT_EC_KEYS = [f'ttl_counter{i}' for i in range(MAX_NUM_CHANNELS)]

    def build(self):
        # Detection DDS and switch
        self._detection_dds = self.get_device('urukul0_ch0', artiq.coredevice.ad9910.AD9910)
        self._detection_sw = self.get_device('ttl_urukul0_sw0', artiq.coredevice.ttl.TTLOut)
        self.update_kernel_invariants('_detection_dds', '_detection_sw')

        # PMT channels
        self._pmt_ttl = [self.get_device(k, artiq.coredevice.ttl.TTLInOut) for k in self._PMT_TTL_KEYS]
        self._pmt_counter = [self.get_device(k, artiq.coredevice.edge_counter.EdgeCounter) for k in self._PMT_EC_KEYS]
        self.update_kernel_invariants('_pmt_ttl', '_pmt_counter')

    def init(self, *, force: bool = False) -> None:
        """Initialize this module.

        :param force: Force full initialization
        """
        # Obtain system datasets
        self._detection_time: float = self.get_dataset_sys(self.DETECTION_TIME_KEY, 200 * us)
        self._detection_freq: float = self.get_dataset_sys(self.DETECTION_FREQ_KEY, 205 * MHz)
        self._state_detection_threshold: int = self.get_dataset_sys(self.STATE_DETECTION_THRESHOLD_KEY, 2)
        self._active_channels: typing.List[np.int32] = self.get_dataset_sys(self.ACTIVE_CHANNELS_KEY,
                                                                            [np.int32(i) for i in range(2)])
        self.update_kernel_invariants('_detection_time', '_detection_freq',
                                      '_state_detection_threshold', '_active_channels')

        if force:
            # Initialize devices
            self.init_kernel()

    @kernel
    def init_kernel(self):
        # For initialization, always reset core first
        self.core.reset()

        # Set direction of PMT TTL pins
        for pmt_ttl in self._pmt_ttl:
            pmt_ttl.input()
            delay_mu(np.int64(self.core.ref_multiplier))  # Added minimal delay to make the events sequential

        # Configure detection DDS
        delay(200 * us)
        self._detection_sw.off()
        self._detection_dds.set(self._detection_freq)
        delay(5 * us)
        self._detection_dds.set_att(17 * dB)

        # Wait until all events have been submitted, always required for initialization
        self.core.wait_until_mu(now_mu())

    def post_init(self):
        pass

    """Module functionality"""

    @kernel
    def detect_channels_mu(self, channels: TList(TInt32), duration: TInt64 = np.int64(0),
                           detection_beam: TBool = True) -> TInt64:
        """Detect ions in parallel using the PMT array (symmetric operation).

        This method returns a timestamp, but this timestamp is not required for obtaining the count.
        Counts can be obtained using the :func:`count` and :func:`measure` functions.

        :param channels: List of PMT channels
        :param duration: Duration of detection in machine units, default value if none given
        :param detection_beam: True to control the detection beam during detection
        :return: The timestamp of the detection window end
        """
        if len(channels) == 0:
            # Check that channel list is not empty
            raise ValueError('Channel list can not be empty')

        if duration <= 0:
            # Use default duration
            duration = self.core.seconds_to_mu(self._detection_time)

        try:
            if detection_beam:
                # Switch detection laser on
                self._detection_sw.on()

            # Perform parallel detection
            for c in channels:
                self._pmt_counter[c].set_config(True, False, False, True)
            delay_mu(duration)
            for c in channels:
                self._pmt_counter[c].set_config(False, False, True, False)

            if detection_beam:
                # Switch detection laser off
                self._detection_sw.off()

        except RTIOUnderflow:
            self.core.reset()
            self._detection_sw.off()
            self.core.wait_until_mu(now_mu())
            raise

        # Return the timestamp (not required to retrieve count from EdgeCounter)
        return now_mu()

    @kernel
    def detect_channels(self, channels: TList(TInt32), duration: TFloat = 0.0,
                        detection_beam: TBool = True) -> TInt64:
        """Detect ions in parallel using the PMT array (symmetric operation).

        This method returns a timestamp, but this timestamp is not required for obtaining the count.
        Counts can be obtained using the :func:`count` and :func:`measure` functions.

        :param channels: List of PMT channels
        :param duration: Duration of detection, default value if none given
        :param detection_beam: True to control the detection beam during detection
        :return: The timestamp of the detection window end
        """
        return self.detect_channels_mu(channels, self.core.seconds_to_mu(duration),
                                       detection_beam=detection_beam)

    @kernel
    def detect_active(self, duration: TFloat = 0.0,
                      detection_beam: TBool = True) -> TInt64:
        """Detect ions using active PMT channels (symmetric operation).

        This method is a convenience function for calling :func:`detect_channels` with only active channels.
        Note that the active channels parameter has to be set earlier.

        :param duration: Duration of detection, default value if none given
        :param detection_beam: True to control the detection beam during detection
        :return: The timestamp of the detection window end
        """
        return self.detect_channels_mu(self._active_channels, self.core.seconds_to_mu(duration),
                                       detection_beam=detection_beam)

    @kernel
    def count(self, channel: TInt32) -> TInt32:
        """Return the PMT count of a specific channel.

        This function can be used in a list comprehension to obtain the counts of a list of channels:

        ``results = [self.detect.count(c) for c in channels]``

        This function can not directly return the array with results due to
        limitations in the compiler (dynamic memory management).

        :param channel: The PMT channel
        :return: The count of the given PMT channel
        """
        return self._pmt_counter[channel].fetch_count()

    @kernel
    def measure(self, channel: TInt32) -> TBool:
        """Read the PMT count of a specific channel and discriminate against the state detection threshold.

        This function can be used in a list comprehension to obtain the measurements of a list of channels:

        ``results = [self.detect.measure(c) for c in channels]``

        This function can not directly return the array with results due to
        limitations in the compiler (dynamic memory management).

        :param channel: The PMT channel
        :return: True if the number of detected events was above the threshold
        """
        return self.count(channel) > self._state_detection_threshold

    """Channel helpers"""

    @portable
    def active_channels(self) -> TList(TInt32):
        """Get the list of active channels.

        :return: List of active channels
        """
        return self._active_channels

    @host_only
    def set_active_channels(self, active_channels: typing.List[np.int32]):
        """Update the list of active channels after ion loading.

        :param active_channels: The list of active channels
        """
        # Store value locally and in the dataset
        assert len(active_channels) <= self.MAX_NUM_CHANNELS
        assert all(0 <= c < self.MAX_NUM_CHANNELS for c in active_channels)
        assert len(set(active_channels)) == len(active_channels)
        self._active_channels = [np.int32(c) for c in active_channels]
        self.set_dataset_sys(self.ACTIVE_CHANNELS_KEY, self._active_channels)

    """Interface functions"""

    @host_only
    def get_pmt_array(self) -> typing.List[artiq.coredevice.edge_counter.EdgeCounter]:
        return self._pmt_counter

    @host_only
    def get_state_detection_threshold(self) -> int:
        return self._state_detection_threshold

    @host_only
    def get_default_detection_time(self) -> float:
        return self._detection_time
