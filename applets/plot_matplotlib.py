import typing

import numpy
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use('Qt5Agg')

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from artiq.applets.simple import TitleApplet
from LAX_exp.applets.widget import QMainWindow


class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None, nsubplots=1):
        self.fig, self.axes = plt.subplots(nsubplots)
        super().__init__(self.fig)


class MatplotlibPlot(QMainWindow):
    def __init__(self, args, req, **kwargs):
        # Call super
        super().__init__(args, req, **kwargs)
        self.sc = MplCanvas(self, nsubplots=args.num_subplots)
        if args.x_label is not None:
            self.sc.fig.supxlabel(args.x_label)
    def update_applet(self, args):

        # grab dataset values
        xs = np.array(self.get_dataset(args.x), None)
        ys = np.array(self.get_dataset(args.y))
        errors = np.array(self.get_dataset(args.error, None))
        fit_xs = np.array(self.get_dataset(args.fit_x, None))
        fit_ys = np.array(self.get_dataset(args.fit_y, None))
        x_labels = numpy.array(self.get_dataset(args.subplot_x_labels, None))
        y_labels = numpy.array(self.get_dataset(args.subplot_y_labels, None))
        titles = np.array(self.get_dataset(args.subplot_titles, None))

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
            if xs is not None:
                if len(xs) != len(ys):
                    raise ValueError("x must have the same length as y")

            if errors is not None:
                if len(errors) != len(ys):
                    raise ValueError("error must have the same length as y")

            if fit_ys is not None:
                if len(fit_ys) != len(ys):
                    raise ValueError("fit_y must have the same length as y")
            else:
                fit_xs = None

            if fit_xs is not None:
                if len(fit_xs) != len(fit_ys):
                    raise ValueError("fit_x must have the same length as fit_y")
            self.plot(xs, ys, errors, fit_xs, fit_ys)

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

                self.plot(x, y, error, fit_x, fit_y, x_label, y_label, title, ind)

        # set the matplotlib plot in the Qt widget and show the widget
        self.setCentralWidget(self.sc)
        self.show()

    def plot(self, x, y, error, fit_x, fit_y, x_label="", y_label="", title="", ind=0):

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

        self.sc.axes[ind].errorbar(x, y, error)
        if fit_y is not None and fit_x is not None:
            self.sc.axes[ind].plot(fit_x, fit_y)

        if title is not None:
            self.sc.axes[ind].set_title(title)
        if x_label is not None:
            self.sc.axes[ind].set_xlabel(x_label)
        if y_label is not None:
            self.sc.axes[ind].set_ylabel(y_label)

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


def main():
    # Create applet object

    applet = TitleApplet(MatplotlibPlot)
    applet.argparser.add_argument("--x-label", type=str, default="", required=False)
    applet.argparser.add_argument("--y-label", type=str, default="", required=False)
    applet.argparser.add_argument("--num-subplots", type=int, default=1, required=False)

    applet.add_dataset("x", "X values")
    applet.add_dataset("y", "Y values")
    applet.add_dataset("error", "Error data (multiple graphs)", required=False)
    applet.add_dataset("fit-x", "X values for fit data", required=False)
    applet.add_dataset("fit-y", "Fit for Y data", required=False)
    applet.add_dataset("subplot-x-labels", "x labels for subplots", required=False)
    applet.add_dataset("subplot-y-labels", "y labels for subplots ", required=False)
    applet.add_dataset("subplot-titles", "title data for subplots", required=False)

    applet.run()

if __name__ == "__main__":
    main()
