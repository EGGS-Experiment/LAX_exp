import numpy
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from PyQt5.QtCore import Qt

from sipyco import pyon
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from artiq.applets.simple import TitleApplet
from LAX_exp.applets.widget import QMainWindow

matplotlib.use('Qt5Agg')


class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None, nsubplots=1):
        # self.fig, self.axes = plt.subplots(nsubplots)
        self.axes = []
        self.fig = plt.figure()
        rows = nsubplots // 2 + 1
        if rows - 1:
            cols = 2
        else:
            cols = 1

        for axis in range(nsubplots):
            if nsubplots % 2 != 0 and axis == nsubplots - 1:
                self.axes.append(plt.subplot2grid(shape=(rows, cols), loc=(axis // 2, axis % 2), colspan=2))
            else:
                self.axes.append(plt.subplot2grid(shape=(rows, cols), loc=(axis//2, axis%2), colspan=1))


        self.axes = np.array(self.axes)

        plt.ion()
        plt.tight_layout()
        self.fig.tight_layout()
        super().__init__(self.fig)


class MatplotlibPlot(QMainWindow):
    """
    Class for creating an interactive Matplotlib plot in an applet

    Attributes:
        sc (MatplotlibCanvas): a matplotlib canvas to hold the plots
        data_points_list (list): list of data points clicked
        boxes (list): list of text boxes on figure
        box_clicked (ax.text): latest text box clicked
        default_keys (list): keys that will be used in the results dictionary passed to this class
    """

    def __init__(self, args, req, **kwargs):
        # Call super
        super().__init__(args, req, **kwargs)
        # create matplotlib canvas to insert plots into
        self.sc = MplCanvas(self, nsubplots=args.num_subplots)

        # instantiate lists of boxes and points for later events
        self.data_points_list = []
        self.fit_lines = []
        self.boxes = []
        self.box_clicked = None

        # connect canvas to button press events
        self.sc.fig.canvas.mpl_connect('button_press_event', self.mouse_event)
        # if args.x_label is not None:
        #     self.sc.fig.supxlabel(args.x_label, fontsize=14)
        # if args.y_label is not None:
        #     self.sc.fig.supylabel(args.y_label, fontsize=14)
        if args.title is not None:
            self.sc.fig.suptitle(args.title, fontsize=18)

        # default keys for dictionary passed to applet
        self.default_keys = ['x', 'y', 'errors', 'fit_x', 'fit_y', 'subplot_x_labels', 'subplot_y_labels',
                             'subplot_titles', 'ylims', 'rid']

    def update_applet(self, args):
        """
        Process run everytime dataset is modified and produces and updates the plot

        Args:
            args (dict): Dictionary of arguments
        """

        # extract data from dictionary
        results = pyon.decode(self.get_dataset(args.results))

        # parse dictionary
        ys = np.array(self.get_from_dict(results, 'y', None), dtype=object)
        xs = np.array(self.get_from_dict(results, 'x', None))
        errors = np.array(self.get_from_dict(results, 'errors', None))
        fit_xs = np.array(self.get_from_dict(results, 'fit_x', None))
        fit_ys = np.array(self.get_from_dict(results, 'fit_y', None))
        x_labels = numpy.array(self.get_from_dict(results, 'subplot_x_labels', None))
        y_labels = numpy.array(self.get_from_dict(results, 'subplot_y_labels', None))
        titles = numpy.array(self.get_from_dict(results, 'subplot_titles', None))
        ylims = numpy.array(self.get_from_dict(results, 'ylims', None))
        rid = numpy.array(self.get_from_dict(results, 'rid', None))

        # inform user of unused keys and data
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

            if ylims is not None and len(ylims.shape) != 0:
                if np.max(ylims.shape) > 1:
                    raise ValueError("There can only be a single y limit for a single subplot")

            self.plot(xs, ys, errors, fit_xs, fit_ys, x_labels.item(), y_labels.item(),
                      titles.item(), ylim = ylims, rid=rid)

        # check if the variables are multi_dimensional that all variables have the same shape
        elif num_datasets >= 2:

            # enumerate through the datasets in y
            for ind, y in enumerate(ys):
                # if a multidimensional array grab the necessary dataset or if one-dimensional just return the variable
                x = self.get_plot_element(xs, ind)
                error = self.get_plot_element(errors, ind)
                fit_x = self.get_plot_element(fit_xs, ind)
                fit_y = self.get_plot_element(fit_ys, ind)
                x_label = self.get_plot_element(x_labels, ind)
                y_label = self.get_plot_element(y_labels, ind)
                title = self.get_plot_element(titles, ind)
                ylim = self.get_plot_element(ylims, ind)

                self.plot(x, y, error, fit_x, fit_y, x_label, y_label, title, ylim, ind, rid=rid)

        # set the matplotlib plot in the Qt widget and show the widget
        self.setCentralWidget(self.sc)
        self.show()

    """PLOTTING FUNCTIONS"""

    def plot(self, x, y, error, fit_x, fit_y, x_label="", y_label="", title="", ylim = None, ind=0, rid=None):
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
            rid: run number of experiment
            ylim: y limits
        """

        # verify the variables are of the correct type
        x = self.verify_aux_arg_type(y, x, "X")
        error = self.verify_aux_arg_type(y, error, "Error")
        fit_y = self.verify_aux_arg_type(fit_y, fit_y, "Fit Y")

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

        # get features of the plot legend
        handles, labels = self.sc.axes[ind].get_legend_handles_labels()

        # only plot if a new experiment is run
        if labels == [] or (rid not in np.int32(np.array(labels))):
            data_points = self.sc.axes[ind].errorbar(x, y, error, marker="o", linestyle="-", label=rid,
                                              markersize = 2)
            self.data_points_list.append(data_points)

        # plot fit to data and set titles and labels
        if fit_y is not None and fit_x is not None:
            fit_lines = self.sc.axes[ind].plot(fit_x, fit_y)
            self.fit_lines.append(fit_lines)
        if title is not None:
            self.sc.axes[ind].set_title(title)
        if x_label is not None:
            self.sc.axes[ind].set_xlabel(x_label)
        if y_label is not None:
            self.sc.axes[ind].set_ylabel(y_label)
        if ylim is not None:
            self.sc.axes[ind].set_ylim(np.min(ylim), np.max(ylim))

        # style plot
        self.sc.axes[ind].ticklabel_format(axis='x', style='plain', useOffset=False)
        self.sc.axes[ind].ticklabel_format(axis='y', style='plain', useOffset=False)
        # set xticks by removing Nones from list
        x_filtered = [x_ele for x_ele in x if x_ele is not None]
        self.sc.axes[ind].set_xticks(np.round(np.linspace(np.min(x_filtered), np.max(x_filtered), 5),3))
        self.sc.axes[ind].legend()
        self.sc.axes[ind].grid(True)

    """VERIFICATION AND HELPER FUNCTIONS"""

    def get_from_dict(self, d, key, default=None):
        if key in list(d.keys()):
            return d[key]
        else:
            return default

    def get_plot_element(self, plot_element, ind):
        """
        Get an element (title, xlabel, data, etc.) from a list or returns the element

        Args:
            plot_element : list of element (data or label)
            ind: index of the label in the list
        """

        # see if plot label comes in array
        if isinstance(plot_element, np.ndarray):
            # check if array is filled with all Nones
            if (plot_element == None).all():
                plot_element = None
            elif len(plot_element.shape) < 1:
                plot_element = plot_element
            # otherwise just grab the label
            elif ind > len(plot_element)-1:
                plot_element = None
            else:
                plot_element = plot_element[ind]
        elif isinstance(plot_element, list) and ind < len(plot_element)-1:
            plot_element = plot_element[ind]

        return plot_element

    def verify_aux_arg_type(self, y, aux_arg, aux_arg_name=""):
        """
        Verify the type of the not required datasets

        Args:
            y: the required dataset
            aux_arg: the not required dataset
            aux_arg_name: the name of the not required dataset
        """
        if aux_arg is not None:
            if (aux_arg == None).all():
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
        """
        Record a mouse click on data point or text box

        Args:
            event (Event): event recording a mouse click
        """
        # if click was on a text box already created then record and don't create new text box
        box_clicked = [box_clicked for box_clicked in self.boxes if box_clicked.contains(event)[0]]
        if len(box_clicked) > 0:
            self.box_clicked = box_clicked[0]

        # if data point was clicked create text box and record
        else:
            self.box_clicked = None
            artist_clicked = [picked.lines[0] for picked in self.data_points_list if picked.lines[0].contains(event)[0]]
            if len(artist_clicked) > 0:
                self.boxes.append(artist_clicked[0].axes.text(event.xdata, event.ydata,
                                                              f"x: {event.xdata:.4f} \ny: {event.ydata:.4f}",
                                                              fontsize=12,
                                                              bbox={'facecolor': 'white', 'pad': 4,
                                                                    'edgecolor': 'black'}))

            # check fit lines
            else:
                self.box_clicked = None
                artist_clicked = [line for lines in self.fit_lines for line in lines if line.contains(event)[0]]

                if len(artist_clicked) > 0:
                    self.boxes.append(artist_clicked[0].axes.text(event.xdata, event.ydata,
                                                                  f"x: {event.xdata:.4f} \ny: {event.ydata:.4f}",
                                                                  fontsize=12,
                                                                  bbox={'facecolor': 'white', 'pad': 4,
                                                                        'edgecolor': 'black'}))

    def keyPressEvent(self, event):
        """
        If a text box was clicked and then the delete or backspace key is pressed, remove the text box

        Args:
            event (Event): An event object recording a key press
        """
        if event.key() == Qt.Key_Backspace or event.key() == Qt.Key_Delete:
            # remove clicked box from list of boxes after key event
            self.boxes = [box for box in self.boxes if box != self.box_clicked]
            # if box was clicked prior to key event remove it
            if self.box_clicked is not None:
                self.box_clicked.remove()
                self.box_clicked = None


def main():
    # Create applet object
    applet = TitleApplet(MatplotlibPlot)

    # get arguments
    applet.argparser.add_argument("--x-label", type=str, default="", required=False)
    applet.argparser.add_argument("--y-label", type=str, default="", required=False)
    applet.argparser.add_argument("--num-subplots", type=int, default=1, required=False)

    # get datasets
    applet.add_dataset("results", "dictionary of experimental results")

    # run applet
    applet.run()


if __name__ == "__main__":
    main()
