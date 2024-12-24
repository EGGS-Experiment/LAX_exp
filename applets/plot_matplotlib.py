import numpy
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from PyQt5.QtCore import Qt
from sipyco.sync_struct import Subscriber
import asyncio
from asyncio import set_event_loop, Event
from threading import Thread

matplotlib.use('Qt5Agg')
from sipyco import pyon

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from artiq.applets.simple import TitleApplet
from LAX_exp.applets.widget import QMainWindow


# class DatasetCollector:
#     def __init__(self, server='::1', port_notify=3250):
#         self.server = server
#         self.port_notify = port_notify
#         self.data = dict()
#         self.loop = asyncio.new_event_loop()
#         asyncio.set_event_loop(self.loop)
#         self.stop_event = Event()
#
#     def sub_init(self, data):
#         self.data = data
#         return self.data
#
#     def sub_mod(self, mod):
#         self.data = mod['struct']
#         self.stop_event.set()
#
#     def subscribe(self):
#         self.subscriber = Subscriber("datasets",
#                                      self.sub_init, self.sub_mod)
#         self.loop.run_until_complete(self.subscriber.connect(
#             self.server, self.port_notify))
#
#     def get_data(self):
#         self.loop.run_until_complete(self.stop_event.wait())

class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None, nsubplots=1):
        self.fig, self.axes = plt.subplots(nsubplots)
        plt.ion()
        super().__init__(self.fig)


class MatplotlibPlot(QMainWindow):
    def __init__(self, args, req, **kwargs):
        # Call super
        super().__init__(args, req, **kwargs)
        self.sc = MplCanvas(self, nsubplots=args.num_subplots)
        self.points = []
        self.sc.fig.canvas.mpl_connect('button_press_event', self.mouse_event)
        self.boxes = []
        self.box_clicked = None
        # if args.x_label is not None:
        #     self.sc.fig.supxlabel(args.x_label, fontsize=14)
        # if args.y_label is not None:
        #     self.sc.fig.supylabel(args.y_label, fontsize=14)
        if args.title is not None:
            self.sc.fig.suptitle(args.title, fontsize=18)
        self.default_keys = ['x', 'y', 'errors', 'fit_x',  'fit_y',  'subplot_x_labels', 'subplot_y_labels',
                             'subplot_titles', 'rid']

    def update_applet(self, args):

        results = pyon.decode(self.get_dataset(args.results))

        ys = np.array(self.get_from_dict(results, 'y', None))
        xs = np.array(self.get_from_dict(results, 'x', None))
        errors = np.array(self.get_from_dict(results, 'errors', None))
        fit_xs = np.array(self.get_from_dict(results, 'fit_x', None))
        fit_ys = np.array(self.get_from_dict(results, 'fit_y', None))
        x_labels = numpy.array(self.get_from_dict(results, 'x_label', None))
        y_labels = numpy.array(self.get_from_dict(results, 'y_label', None))
        titles = numpy.array(self.get_from_dict(results, 'title', None))
        rid = numpy.array(self.get_from_dict(results, 'rid', None))
        # data = np.array(self.get_dataset(args.data))
        # if len(data) == 2:
        #     xs, ys = data
        #     errors = np.array(None)
        # elif len(data) == 3:
        #     xs, ys, errors = data
        # else:
        #     raise ValueError("data should be (2,N,M) or (3,N,M) array contain x,y, and errors (optional)")

        # ys = np.array(self.get_dataset(args.y))
        # xs = np.array(self.get_dataset(args.x, None))
        # errors = np.array(self.get_dataset(args.error, None))
        # fit_xs = np.array(self.get_dataset(args.fit_x, None))
        # fit_ys = np.array(self.get_dataset(args.fit_y, None))
        # x_labels = numpy.array(self.get_dataset(args.subplot_x_labels, None))
        # y_labels = numpy.array(self.get_dataset(args.subplot_y_labels, None))
        # titles = np.array(self.get_dataset(args.subplot_titles, None))
        # rid = np.array(self.get_dataset(args.rid, None))

        unused_keys = [key for key in results.keys() if key not in self.default_keys]
        for key in unused_keys:
            self.logger.warning(f'Key {key} is listed in results dictionary for plotting, but is unused')
            self.logger.warning(f'The keys used are: {self.default_keys}')

        # determine number of datasets there are
        if len(ys.shape) == 1:
            num_datasets = 1
        else:
            num_datasets = ys.shape[0]

        # ensure number of datasets is equal to number of subplots
        if num_datasets != args.num_subplots:
            raise ValueError("Number of datasets does not match number of subplots")

        if num_datasets == 1:

            # ensure all variables are the same size as x or are NoneType
            if xs is not None and len(xs.shape) != 0:
                if len(xs) != len(ys):
                    raise ValueError("x must have the same length as y")

            if errors is not None and len(errors.shape) != 0:
                if len(errors) != len(ys):
                    raise ValueError("error must have the same length as y")

            if fit_ys is not None and len(fit_ys.shape) != 0:
                if len(fit_ys) != len(ys):
                    raise ValueError("fit_y must have the same length as y")
            else:
                fit_xs = None

            if fit_xs is not None and len(fit_xs.shape) != 0:
                if len(fit_xs) != len(fit_ys):
                    raise ValueError("fit_x must have the same length as fit_y")

            if x_labels is not None and len(x_labels.shape) != 0:
                if np.max(x_labels.shape) > 1:
                    raise ValueError("There can only be a single x label for a single subplot")
            if y_labels is not None and len(y_labels.shape) != 0:
                if np.max(y_labels.shape) > 1:
                    raise ValueError("There can only be a single y label for a single subplot")

            if titles is not None and len(titles.shape) != 0:
                if np.max(titles.shape) > 1:
                    raise ValueError("There can only be a single title for a single subplot")

            self.plot(xs, ys, errors, fit_xs, fit_ys, x_labels.item(), y_labels.item(), titles.item(), rid=rid)

        # check if the variables are multi_dimensional that all variables have the same shape
        elif num_datasets >= 2:

            # enumerate through the datasets in y
            for ind, y in enumerate(ys):
                # if a multidimensional array grab the necessary dataset or if one-dimensional just return the variable
                x = self.get_plot_data(xs, ind)
                error = self.get_plot_data(errors, ind)
                fit_x = self.get_plot_data(fit_xs, ind)
                fit_y = self.get_plot_data(fit_ys, ind)
                x_label = self.get_plot_label(x_labels, ind)
                y_label = self.get_plot_label(y_labels, ind)
                title = self.get_plot_label(titles, ind)

                self.plot(x, y, error, fit_x, fit_y, x_label, y_label, title, ind, rid=rid)

        # set the matplotlib plot in the Qt widget and show the widget
        self.setCentralWidget(self.sc)
        self.show()

    """PLOTTING FUNCTIONS"""

    def plot(self, x, y, error, fit_x, fit_y, x_label="", y_label="", title="", ind=0, rid=False):

        """
        Plot the data in a matplotlib subplots

        Args:
            x : x-data
            y : y-data
            error : error bars
            fit_x : domain for the fitted data
            fit_y : the fit for the data
            x_label : label for the x-axis of the subplot
            y_label : label for the y-axis of the subplot
            title : title for the subplot
            ind: index (location) of the subplot in the matplotlib figure
        """

        # verify the variables are of the correct type
        x = self.verify_aux_arg_type(y, x, "X")
        error = self.verify_aux_arg_type(y, error, "Error")
        fit_y = self.verify_aux_arg_type(y, fit_y, "Fit Y")
        fit_x = self.verify_aux_arg_type(y, fit_x, "Fit X")

        if fit_y is not None and fit_x is not None:
            if len(fit_y) > len(fit_x):
                fit_y = None
            elif len(fit_x) > len(fit_y):
                # trim fit x data
                fit_x = fit_x[:len(fit_y)]
        elif fit_y is not None and x is not None:
            fit_x = x
        else:
            fit_x = None

        if not isinstance(self.sc.axes, numpy.ndarray):
            self.sc.axes = np.array([self.sc.axes])

        # only plot if a new rid is update (prevents legend from becoming flooded)
        handles, labels = self.sc.axes[ind].get_legend_handles_labels()
        if labels == [] or (rid not in np.int32(np.array(labels))):
            line = self.sc.axes[ind].errorbar(x, y, error, marker="o", linestyle="-", label=rid)
            self.points.append(line)
        if fit_y is not None and fit_x is not None:
            self.sc.axes[ind].plot(fit_x, fit_y)
        if title is not None:
            self.sc.axes[ind].set_title(title)
        if x_label is not None:
            self.sc.axes[ind].set_xlabel(x_label)
        if y_label is not None:
            self.sc.axes[ind].set_ylabel(y_label)

        self.sc.axes[ind].ticklabel_format(axis='x', style='plain', useOffset=False)
        self.sc.axes[ind].ticklabel_format(axis='y', style='plain', useOffset=False)

        self.sc.axes[ind].legend()

    """VERIFICATION AND HELPER FUNCTIONS"""

    def get_from_dict(self, d, key, default=None):
        if key in list(d.keys()):
            return d[key]
        else:
            return default

    def get_plot_data(self, plot_data, ind):
        """
        Get a dataset from a mutltidimensional list or return the variable

        Args:
            plot_data : multidimensional list of datasets
            ind: index of the dataset in the multidimensional list
        """
        if len(plot_data.shape) > 1:
            plot_data = plot_data[ind]
        return plot_data

    def get_plot_label(self, plot_label, ind):
        """
        Get a label (title, xlabel, etc.) from a list or returns the label

        Args:
            plot_label : list of labels
            ind: index of the label in the list
        """

        # see if plot label comes in array
        if isinstance(plot_label, np.ndarray):
            # check if array is filled with any Nones
            if (plot_label == None).any():
                plot_label = None
            elif len(plot_label.shape) < 1:
                plot_label = plot_label
            # otherwise just grab the label
            else:
                plot_label = plot_label[ind]
        if isinstance(plot_label, list):
            plot_label = plot_label[ind]
        return plot_label

    def verify_aux_arg_type(self, y, aux_arg, aux_arg_name=""):
        """
        Verify the type of the not required datasets

        y: the required dataset
        aux_arg: the not required dataset
        aux_arg_name: the name of the not required dataset
        """
        if aux_arg is not None:
            if (aux_arg == None).any():
                aux_arg = None

            elif not len(aux_arg):
                aux_arg = None
            else:
                if aux_arg.shape != y.shape:
                    self.logger.warning(f'{aux_arg_name} data does not match shape of Y data')
                    aux_arg = None

        return aux_arg

    """EVENT HANDLING"""

    def mouse_event(self, event):
        axes = self.sc.axes.tolist()
        if not event.inaxes:
            return

        box_clicked = [box_clicked for box_clicked in self.boxes if box_clicked.contains(event)[0]]
        if len(box_clicked) > 0:
            self.box_clicked = box_clicked[0]
        else:
            self.box_clicked = None
            artist_clicked = [picked.lines[0] for picked in self.points if picked.lines[0].contains(event)[0]]
            if len(artist_clicked) > 0:
                artist = artist_clicked[0]
                axis = artist.axes
                self.boxes.append(axis.text(event.xdata, event.ydata, f"x: {event.xdata:.2f} \ny: {event.ydata:.2f}",
                                            fontsize=12,
                                            bbox={'facecolor': 'white', 'pad': 4, 'edgecolor': 'black'}))

    def keyPressEvent(self, event):
        if event.key() == Qt.Key_Backspace or event.key() == Qt.Key_Delete:
            self.boxes = [box for box in self.boxes if box != self.box_clicked]
            if self.box_clicked is not None:
                self.box_clicked.remove()
                self.box_clicked = None


def main():
    # Create applet object

    applet = TitleApplet(MatplotlibPlot)

    # applet.argparser.add_argument("--x", type=str, default=None, help="X values", required=False)
    applet.argparser.add_argument("--x-label", type=str, default="", required=False)
    applet.argparser.add_argument("--y-label", type=str, default="", required=False)
    applet.argparser.add_argument("--num-subplots", type=int, default=1, required=False)
    # applet.argparser.add_argument("--subplot-x-labels", type=str, default=None, help="x labels for subplots",
    #                               required=False)
    # applet.argparser.add_argument("--subplot-y-labels", type=str, default=None, help="y labels for subplots",
    #                               required=False)
    # applet.argparser.add_argument("--subplot-titles", type=str, default=None, help="title data for subplots",
    #                               required=False)
    # applet.argparser.add_argument("--rid", type=str, default=None, help="rid for experiment", required=False)
    # applet.argparser.add_argument("--error", type=str, default=None, help="Error data (multiple graphs)",
    #                               required=False)
    # applet.argparser.add_argument("--fit-x", type=str, default=None, help="X values for fit data", required=False)
    # applet.argparser.add_argument("--fit-y", type=str, default=None, help="Fit for Y data", required=False)

    applet.add_dataset("results", "experimental results")
    # applet.add_dataset("y", "Y values")
    # applet.add_dataset("x", "X values")
    # applet.add_dataset('data', 'data values')
    # applet.add_dataset("subplot-x-labels", "x labels for subplots",  required=False)
    # applet.add_dataset("subplot-y-labels", "y labels for subplots", required=False)
    # applet.add_dataset("subplot-titles", "title data for subplots", required=False)
    # # applet.add_dataset("error", "Error data (multiple graphs)", required=False)
    # applet.add_dataset("fit-x", "X values for fit data", required=False)
    # applet.add_dataset("fit-y", "Fit for Y data", required=False)
    # applet.add_dataset("rid", "Experiment RID", required=False)

    applet.run()


if __name__ == "__main__":
    main()
