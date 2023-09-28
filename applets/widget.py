"""
Forked from CampbellGroup @ UCLA: https://github.com/CampbellGroup/Yaax
"""
import typing
import logging

import PyQt5  # noqa: F401
import pyqtgraph

__all__ = ['PlotWidget']

# Initialize logging infrastructure
logging.basicConfig(level=logging.WARNING)


class NoDefault:
    """Class used to indicate that no default was set."""
    pass


class SkipUpdateException(KeyError):
    """Raised when this update will be skipped."""
    pass


class PlotWidget(pyqtgraph.PlotWidget):
    """Minor extension over the regular PlotWidget with a few extra conveniences.

    Documentation about plotting functions can be found at:
    https://pyqtgraph.readthedocs.io/en/latest/graphicsItems/plotitem.html

    Function calls to `self` are dynamically forwarded to the underlying plot item.
    For direct access to the plot item, use `self.plotItem` or call `self.getPlotItem()`.
    """

    def __init__(self, args: typing.Any, **kwargs: typing.Any):
        """The init function as it is called by ARTIQ.

        :param args: The arguments from argparse
        :param kwargs: Keyword arguments forwarded to the superclass
        """

        # Call super
        pyqtgraph.PlotWidget.__init__(self, **kwargs)
        # Store arguments
        self.__args = args  # type: typing.Any
        # Title, acts as a buffer
        self.__title = None  # type: typing.Optional[str]
        # Flag to indicate the title was set
        self.__title_flag = False  # type: bool
        # Dict with data, acts as a buffer
        self.__data_buffer = {}  # type: typing.Dict[str, typing.Any]
        # Logger
        self.__logger = logging.getLogger(self.__class__.__name__)

    @property
    def logger(self) -> logging.Logger:
        """Logger object."""
        return self.__logger

    def update_applet(self, args):
        """This function replaces the :func:`data_changed` function and will be called whenever data changes.

        Originally the :func:`data_changed` function would be implemented for custom applets.
        Now the :func:`data_changed` function provides generic functionality and instead
        custom applets should overwrite this method.

        The signature of this method changed and data can now be
        accessed using the :func:`get_data` function and the title is already set.

        :param args: The arguments object returned by argparse
        """
        raise NotImplementedError  # Not using abc to prevent metaclass conflict

    def get_dataset(self, key: str, default: typing.Any = NoDefault) -> typing.Any:
        """Get data from the latest buffer.

        If the data is not available and no default was set, the update function will gracefully return.

        :param key: The key
        :param default: A default value if desired
        :return: The requested value or the default if given
        """
        try:
            # Try to return the data
            return self.__data_buffer[key][1]  # Extra index required
        except KeyError:
            if default is NoDefault:
                # Raise if no default was set
                raise SkipUpdateException
            else:
                return default

    def data_changed(self, data: typing.Dict[str, typing.Any], mods: typing.List[typing.Any],
                     title: typing.Optional[str] = None) -> None:
        """This function is called when a subscribed dataset changes.

        It now provides some standard functionality and custom applets should override
        the :func:`update_applet` function.

        :param data: Raw data in the form of a dict
        :param mods: A list of unknown objects
        :param title: The title, if this is a TitleApplet and the title was set
        """
        assert isinstance(data, dict)
        assert isinstance(mods, list)
        assert isinstance(title, str) or title is None

        # Store data in the buffer
        self.__data_buffer = data
        # Store the formatted title
        self.__title = title
        # Set the title flag to false
        self.__title_flag = False

        if self.plotItem is not None:
            try:
                # Call the update function
                self.update_applet(self.__args)
            except SkipUpdateException:
                pass
            else:
                if not self.__title_flag:
                    # Title was not set, set the formatted title
                    self.setTitle(title)

    def extend_title(self, text: str) -> None:
        """Extend the title.

        :param text: The text to extend the title with
        """
        if self.__title is not None:
            text = '{}:&nbsp;&nbsp;&nbsp;&nbsp;{}'.format(self.__title, text)
        self.setTitle(text)

    # noinspection PyPep8Naming
    def setTitle(self, title: typing.Optional[str] = None, **kwargs: typing.Any) -> None:
        """Override the default set title function (call is forwarded to the underlying plot item)."""
        if title is not None:
            # Add formatting
            title = '<h1>{}</h1>'.format(title)

        # Set title on underlying plot item
        self.plotItem.setTitle(title, **kwargs)
        # Set title flag
        self.__title_flag = True