#!/usr/bin/env python3
"""
Forked from CampbellGroup @ UCLA: https://github.com/CampbellGroup/Yaax
"""
import pyqtgraph
import numpy as np

from artiq.tools import scale_from_metadata
from artiq.applets.simple import TitleApplet
from LAX_exp.applets.widget import PlotWidget


class MultiXYPlot(PlotWidget):
    """Plot XY applet for multiple graphs."""

    def __init__(self, args, req, **kwargs):
        # Call super
        super().__init__(args, req, **kwargs)

        # Set plot name format
        self._name_format = args.plot_names
        if type(self._name_format) is not list:
            self._name_format = [self._name_format]
        # if '{index}' not in args.plot_names:
        #     # Append index
        #     self._name_format += ' {index}'

        # Set labels
        self.plotItem.setLabel('bottom', args.x_label)
        self.plotItem.setLabel('left', args.y_label)

        # Disable mouse tracking (mitigates qt hover bug)
        self.setMouseTracking(False)
        # Enable legend
        self.plotItem.addLegend()

    def update_applet(self, args):
        # Obtain input data
        y = np.asarray(self.get_dataset(args.y))
        x = np.asarray(self.get_dataset(args.x, np.arange(len(y))))
        error = self.get_dataset(args.error, None)
        fit = self.get_dataset(args.fit, None)
        fit_x = self.get_dataset(args.fit_x, None)
        h_lines = self.get_dataset(args.h_lines, [])
        v_lines = self.get_dataset(args.v_lines, [])

        # Verify input data
        if not len(y) or len(y) > len(x):
            return
        if y.ndim != 2:
            raise ValueError('Y data must have two dimensions')
        if len(x) > len(y):
            # Trim x data
            x = x[:len(y)]

        if error is not None:
            if not len(error):
                error = None
            else:
                error = np.asarray(error)
                if error.shape != y.shape:
                    self.logger.warning('Error data does not match shape of Y data')
                    error = None
        if fit is not None:
            if not len(fit):
                fit = None
            else:
                fit = np.asarray(fit)
                if fit_x is None:
                    if fit.shape != y.shape:
                        self.logger.warning('Fit data does not match shape of Y data')
                        fit = None
                    else:
                        fit_x = x
                else:
                    fit_x = np.asarray(fit_x)
                    if len(fit) > len(fit_x):
                        self.logger.warning('Fit data does not match shape of fit X data')
                        fit = None
                        fit_x = None
                    elif len(fit_x) > len(fit):
                        # trim fit x data
                        fit_x = fit_x[:len(fit)]

        # Handle sliding window
        if args.sliding_window > 0:
            # Get window size
            window_size = args.sliding_window

            # Truncate input data based on the window size
            v_lines = [v for v in v_lines if len(y) - window_size <= v < len(y)]
            y = y[-window_size:]
            x = x[-window_size:]
            if error is not None:
                error = error[-window_size:]
            if fit is not None:
                win = (fit_x >= np.min(x)) & (fit_x <= np.max(x))
                fit = fit[win]
                fit_x = fit_x[win]

        if args.index:
            # Indices are provided, only plot selected sequences
            graph_index = args.index
            y = y[:, graph_index]
            if error is not None:
                error = error[:, graph_index]
            if fit is not None:
                fit = fit[:, graph_index]
        else:
            graph_index = np.arange(y.shape[1])

        if args.subsample > 1:
            # Get subsampling factor
            factor = args.subsample

            # Subsample data
            y = y[::factor]
            x = x[::factor]
            if error is not None:
                error = error[::factor]
            if fit is not None:
                fit = fit[::factor]
                fit_x = fit_x[::factor]

        if args.multiplier is not None:
            # Apply multiplier on data
            y *= args.multiplier

        # Clear plot
        self.plotItem.clear()

        # Plot horizontal and vertical lines
        for h in h_lines:
            self.plotItem.addLine(y=h)
        for v in v_lines:
            self.plotItem.addLine(x=v)

        if error is not None:
            for y_values, error_values in zip(y.transpose(), error.transpose()):
                # Plot error bars (note: https://github.com/pyqtgraph/pyqtgraph/issues/211)
                self.plotItem.addItem(pyqtgraph.ErrorBarItem(x=x, y=y_values, height=error_values))

        if fit is not None:
            # Sort based on x data
            idx = fit_x.argsort()
            for fit_values, i in zip(fit.transpose(), graph_index):
                # Pick a color and plot graph
                color = pyqtgraph.intColor(i)
                self.plotItem.plot(x=fit_x[idx], y=fit_values[idx], pen=color)


        for y_values, i in zip(y.transpose(), graph_index):
            # Assemble name of the plot
            try:
                name = self._name_format[i]
            except IndexError:
                name = self._name_format[0]
                # ensure multiply-listed named are differentiated via index
                if '{index}' not in name:
                    name += ' {index}'
                name = name.format(index=i)

            # Pick a color and plot graph
            color = pyqtgraph.intColor(i)
            pen = pyqtgraph.mkPen(color, width=1)
            self.plotItem.plot(x=x, y=y_values, name=name, pen=pen, connect='all')


def main():
    # Create applet object
    applet = TitleApplet(MultiXYPlot, default_update_delay=0.1)

    # Add custom arguments
    applet.argparser.add_argument("--x-label", default=None, type=str, help="The X label")
    applet.argparser.add_argument("--y-label", default=None, type=str, help="The Y label")
    applet.argparser.add_argument("--sliding-window", default=0, type=int,
                                  help="Only show the latest data points")
    applet.argparser.add_argument("--subsample", default=0, type=int,
                                  help="Subsample data by the given factor before plotting")
    applet.argparser.add_argument("--multiplier", default=None, type=float,
                                  help="Multiply data before plotting")
    applet.argparser.add_argument("--index", nargs='+', default=[], type=int,
                                  help="Only plot the graphs at the given indices (plots all by default)")
    applet.argparser.add_argument("--plot-names", nargs="+", default='Plot', type=str,
                                  help="Prefix of numbered plot names (formatting with `{index}` possible)")

    # Add datasets
    applet.add_dataset("y", "Y data (multiple graphs)")
    applet.add_dataset("x", "X values", required=False)
    applet.add_dataset("error", "Error bars for Y data (multiple graphs)", required=False)
    applet.add_dataset("fit", "Fit for Y data (multiple graphs)", required=False)
    applet.add_dataset("fit-x", "X values for fit data", required=False)
    applet.add_dataset("v-lines", "Vertical lines", required=False)
    applet.add_dataset("h-lines", "Horizontal lines", required=False)
    applet.run()


if __name__ == "__main__":
    main()
