import numpy as np
import typing

from dax.experiment import *
from dax.modules.time_resolved_context import TimeResolvedContext

from dax_example.modules.detection import DetectionModule


class TimeResolvedService(DaxService):
    SERVICE_NAME = 'time_resolved'

    # System dataset keys
    PARTITION_SIZE_KEY: str = 'partition_size'

    def build(self):
        # Obtain required modules
        self._detect = self.registry.find_module(DetectionModule)
        self.update_kernel_invariants('_detect')

        # Create time resolved context
        self._time_resolved_context = TimeResolvedContext(self, 'time_resolved_context',
                                                          plot_base_key='{scheduler.rid}', plot_group_base_key="")
        self.update_kernel_invariants('context')

    def init(self):
        # Get system datasets
        self._max_partition_size: int = self.get_dataset_sys(self.PARTITION_SIZE_KEY, 64)

    def post_init(self):
        pass

    """Service functionality"""

    @property
    def context(self) -> TimeResolvedContext:
        """Return the time resolved context object.

        This context can be used with a `with` statement.
        Only inside this context it is possible to use the
        :func:`detect` and :func:`count` type functions.

        This context can be used inside or outside kernel context
        and relies on async RPC calls for enter and exit.

        The time resolved context can be further configured by calling its functions.

        :return: The time resolved context object
        """
        return self._time_resolved_context

    @kernel
    def detect_channels(self, channels: TList(TInt32),
                        num_bins: TInt32, bin_width: TFloat, bin_spacing: TFloat, offset: TFloat = 0.0,
                        detection_beam: TBool = False):
        """Time resolved detection on specific channels.

        Note that the offset delay is applied in this function!

        :param channels: List of PMT channels
        :param num_bins: The number of bins
        :param bin_width: The width of each individual bin
        :param bin_spacing: The spacing between bins
        :param offset: The known fixed offset of this trace in seconds, used for partitioning (defaults to no offset)
        :param detection_beam: True to control the detection beam during detection
        """
        # Convert bin width and spacing to machine units
        bin_width_mu = self.core.seconds_to_mu(bin_width)
        bin_spacing_mu = self.core.seconds_to_mu(bin_spacing)

        if bin_spacing_mu < self.core.ref_multiplier:
            # Guarantee a minimum bin spacing to separate detection windows
            # Subtract the difference from the bin width to not change the bin start times
            bin_width_mu -= self.core.ref_multiplier - bin_spacing_mu
            bin_spacing_mu = np.int64(self.core.ref_multiplier)

        # Apply offset
        delay(offset)

        for _ in range(num_bins - 1):
            # Detect and add spacing for the bins
            self._detect.detect_channels_mu(channels, bin_width_mu, detection_beam=detection_beam)
            delay_mu(bin_spacing_mu)

        # No spacing after the last bin
        self._detect.detect_channels_mu(channels, bin_width_mu, detection_beam=detection_beam)

        # Store metadata in context object
        self.context.append_meta(bin_width, bin_spacing, offset=offset)

    @kernel
    def detect(self, channel: TInt32, num_bins: TInt32, bin_width: TFloat, bin_spacing: TFloat, offset: TFloat = 0.0,
               detection_beam: TBool = False):
        """Time resolved detection on one channel.

        Note that the offset delay is applied in this function!

        :param channel: The target channel
        :param num_bins: The number of bins
        :param bin_width: The width of each individual bin
        :param bin_spacing: The spacing between bins
        :param offset: The known fixed offset of this trace in seconds, used for partitioning (defaults to no offset)
        :param detection_beam: True to control the detection beam during detection
        """
        self.detect_channels([channel], num_bins, bin_width, bin_spacing, offset=offset, detection_beam=detection_beam)

    @kernel
    def remove(self):
        """Removes the metadata of the last call to a ``detect*`` function.

        When a ``detect*`` function is called, metadata of the detection is stored.
        In case a matching ``count*`` function can not be called because data is inconsistent,
        the metadata can be removed by calling this function.
        Queued data must be consistent, and every call to ``detect*`` must therefore have a
        matching call to a ``count*`` function or the :func:`remove` function.

        The implementation of the :func:`remove` function leaves it up to the user to decide
        if the last detection is repeated to still yield a data point, or if the data point
        is dropped. Note that if a data point is dropped, results from the time resolved context
        might be inconsistent with other data series (e.g. scannables).
        """
        self.context.remove_meta()

    @kernel
    def count_channels(self, channels: TList(TInt32), num_bins: TInt32, offset_mu: TInt64 = np.int64(0)):
        """Record the PMT counts of a list of channels.

        The count values are requested from the detection module and
        the results are stored in the time resolved context buffer.

        Note that corrections for delayed events should result in **negative offset**.
        The negative offset represents the fact that detection started before the event happened.

        :param channels: The channels to record the counts of
        :param num_bins: The number of bins
        :param offset_mu: An offset to correct any shifts of events in machine units (defaults to no offset)
        """

        # Append the list of detection counts to the time resolved context buffer
        self.context.append_data([[self._detect.count(c) for _ in range(num_bins)] for c in channels],
                                 offset_mu=offset_mu)

    @kernel
    def count(self, channel: TInt32, num_bins: TInt32, offset_mu: TInt64 = np.int64(0)):
        """Record the PMT counts of one channel.

        The count values are requested from the detection module and
        the results are stored in the time resolved context buffer.

        Note that corrections for delayed events should result in **negative offset**.
        The negative offset represents the fact that detection started before the event happened.

        :param channel: The target channel
        :param num_bins: The number of bins
        :param offset_mu: An offset to correct any shifts of events in machine units (defaults to no offset)
        """

        # Append the list of detection counts to the time resolved context buffer
        self.count_channels([channel], num_bins, offset_mu=offset_mu)

    """Partitioning functions"""

    @host_only
    def partition_bins(self, num_bins: int,
                       bin_width: float, bin_spacing: float) -> typing.List[typing.Tuple[np.int32, float]]:
        """Partition a number of bins.

        This function returns a list of tuples that can be used at runtime for partitioning in a loop.
        The format of each element is ``(current_num_bins, current_offset)`` which can be used accordingly.

        This function returns the partition table as a list. Hence, it can only be called from the host.
        This module must be initialized for partitioning to be available. Hence, this function will
        most often be called immediately after ``dax_init()``.

        :param num_bins: The total number of bins desired
        :param bin_width: The width of each bin
        :param bin_spacing: The spacing between bins
        :return: A list with tuples that can be used for automatic partitioning at runtime
        """
        if not self.hasattr('_max_partition_size'):
            raise RuntimeError('dax_init() needs to be called before requesting a partition table')
        return self.context.partition_bins(num_bins, self._max_partition_size, bin_width, bin_spacing)

    @host_only
    def partition_window(self, window_size: float,
                         bin_width: float, bin_spacing: float) -> typing.List[typing.Tuple[np.int32, float]]:
        """Partition a time window.

        This function returns a list of tuples that can be used at runtime for partitioning in a loop.
        The format of each element is ``(current_num_bins, current_offset)`` which can be used accordingly.

        This function returns the partition table as a list. Hence, it can only be called from the host.
        This module must be initialized for partitioning to be available. Hence, this function will
        most often be called immediately after ``dax_init()``.

        :param window_size: The total window size
        :param bin_width: The width of each bin
        :param bin_spacing: The spacing between bins
        :return: A list with tuples that can be used for automatic partitioning at runtime
        """
        if not self.hasattr('_max_partition_size'):
            raise RuntimeError('dax_init() needs to be called before requesting a partition table')
        return self.context.partition_window(window_size, self._max_partition_size, bin_width, bin_spacing)
