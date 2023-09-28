#!/usr/bin/env python3

"""
Forked from CampbellGroup @ UCLA: https://github.com/CampbellGroup/Yaax
"""

import numpy as np
import PyQt5  # make sure pyqtgraph imports Qt5
from PyQt5.QtCore import QTimer
import pyqtgraph

from artiq.applets.simple import TitleApplet


class DifferentialModePlot(pyqtgraph.PlotWidget):
    def __init__(self, args):
        pyqtgraph.PlotWidget.__init__(self)
        self.args = args
        self.timer = QTimer()
        self.timer.setSingleShot(True)
        self.timer.timeout.connect(self.length_warning)
        self.mismatch = {'X values': False,
                         'Error bars': False,
                         'Fit values': False}

        self.darkcounts = self.addItem(self.ScatterPlotItem)
        self.brightcounts = self.addItem(self.ScatterPlotItem)
        self.diffcounts = self.addItem(self.ScatterPlotItem)

    def data_changed(self, data, mods, title):
        try:
            y = data[self.args.y][1]
        except KeyError:
            return

        x = data.get(self.args.x, (False, None))[1]
        if x is None:
            x = np.arange(len(y[:, 0]))
        error = data.get(self.args.error, (False, None))[1]
        fit = data.get(self.args.fit, (False, None))[1]

        if not len(y[:, 0]) or len(y[:, 0]) != len(x):
            self.mismatch['X values'] = True
        else:
            self.mismatch['X values'] = False
        if error is not None and hasattr(error, "__len__"):
            if not len(error):
                error = None
            elif len(error) != len(y):
                self.mismatch['Error bars'] = True
            else:
                self.mismatch['Error bars'] = False
        if fit is not None:
            if not len(fit):
                fit = None
            elif len(fit) != len(y):
                self.mismatch['Fit values'] = True
            else:
                self.mismatch['Fit values'] = False
        if not any(self.mismatch.values()):
            self.timer.stop()
        else:
            if not self.timer.isActive():
                self.timer.start(1000)
            return

        self.darkcounts.clear()
        self.darkcounts.setData(x, y[:, 0], pen=None, symbol='o', symbolBrush='y')

        self.brightcounts.clear()
        self.brightcounts.setData(x, y[:, 1], pen=None, symbol='o', symbolBrush='r')

        self.diffcounts.clear()
        self.diffcounts.setData(x, y[:, 2], pen=None, symbol='x', symbolBrush='b')

        # self.plot(x, y[:,0], pen=None, symbol="x")
        # self.plot(x, y[:,1], pen=None, symbol="x")
        # self.plot(x, y[:,2], pen=None, symbol="o")

        self.setTitle(title)
        if error is not None:
            # See https://github.com/pyqtgraph/pyqtgraph/issues/211
            if hasattr(error, "__len__") and not isinstance(error, np.ndarray):
                error = np.array(error)
            errbars = pyqtgraph.ErrorBarItem(
                x=np.array(x), y=np.array(y), height=error)
            self.addItem(errbars)
        if fit is not None:
            xi = np.argsort(x)
            self.plot(x[xi], fit[xi])

    def length_warning(self):
        self.clear()
        text = "⚠️ dataset lengths mismatch:\n"
        errors = ', '.join([k for k, v in self.mismatch.items() if v])
        text = ' '.join([errors, "should have the same length as Y values"])
        self.addItem(pyqtgraph.TextItem(text))


def main():
    applet = TitleApplet(DifferentialModePlot)
    applet.add_dataset("y", "Y values")
    applet.add_dataset("x", "X values", required=False)
    applet.add_dataset("error", "Error bars for each X value", required=False)
    applet.add_dataset("fit", "Fit values for each X value", required=False)
    applet.run()


if __name__ == "__main__":
    main()