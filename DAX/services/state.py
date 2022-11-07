from dax.experiment import *
from dax.modules.hist_context import HistogramContext

from dax_example.modules.detection import DetectionModule


class StateService(DaxService):
    SERVICE_NAME = 'state'

    # System dataset keys
    INIT_TIME_KEY = 'init_time'

    def build(self) -> None:
        # Obtain required modules
        self._detect = self.registry.find_module(DetectionModule)
        self.update_kernel_invariants('_detect')

        # Create histogram context
        self._histogram_context = HistogramContext(self, 'hist_context',
                                                   plot_base_key='{scheduler.rid}', plot_group_base_key="")
        self.update_kernel_invariants('context')

    def init(self) -> None:
        pass

    def post_init(self) -> None:
        pass

    """Service functionality"""

    @property
    def context(self) -> HistogramContext:
        """Return the histogram context object.

        This context can be used with a `with` statement.
        Only inside this context it is possible to use the
        :func:`count` and :func:`measure` type functions.

        This context can be used inside or outside kernel context
        and relies on async RPC calls for enter and exit.

        The histogram context can be further configured by calling its functions.

        :return: The histogram context object
        """
        return self._histogram_context

    @kernel
    def detect_channels(self, channels: TList(TInt32), duration: TFloat = 0.0):
        """Detect state on specific channels.

        Call is forwarded to detection module.

        :param channels: List of PMT channels
        :param duration: Duration of detection, default value if none given
        """
        self._detect.detect_channels(channels, duration)

    @kernel
    def detect_active(self, duration: TFloat = 0.0):
        """Detect state on active channels.

        Call is forwarded to detection module.

        :param duration: Duration of detection, default value if none given
        """
        self._detect.detect_active(duration)

    @kernel
    def count_channels(self, channels: TList(TInt32)):
        """Record the PMT counts of a list of channels.

        The count values are requested from the detection module and
        the results are stored in the histogram buffer.

        :param channels: The channels to record the counts of
        """

        # Append the list of detection counts to the histogram buffer
        self.context.append([self._detect.count(c) for c in channels])

    @kernel
    def count_active(self):
        """Record the PMT counts of active channels.

        The count values are requested from the detection module and
        the results are stored in the histogram buffer.
        """

        # Append the list of detection counts to the histogram buffer
        self.context.append([self._detect.count(c) for c in self._detect.active_channels()])

    @kernel
    def measure_channels(self, channels: TList(TInt32)):
        """Record the PMT counts of a list of channels discriminated against the state detection threshold.

        The count values are requested from the detection module and
        the results are discriminated and stored in the histogram buffer.

        :param channels: The channels to record the measurements of
        """
        self.context.append([self._detect.measure(c) for c in channels])

    @kernel
    def measure_active(self):
        """Record the PMT counts of active channels discriminated against the state detection threshold.

        The count values are requested from the detection module and
        the results are discriminated and stored in the histogram buffer.
        """
        self.context.append([self._detect.measure(c) for c in self._detect.active_channels()])
