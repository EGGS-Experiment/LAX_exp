#!/usr/bin/env python3

from PyQt5 import QtWidgets, QtCore, QtGui

from artiq.gui.tools import LayoutWidget
from artiq.applets.simple import SimpleApplet
from artiq.applets.progress_bar import ProgressWidget

from time import time


class CompletionMonitor(LayoutWidget):
    """
    artiq.applets.progress_bar but with expected completion time.
    """

    def __init__(self, args, req):
        LayoutWidget.__init__(self)
        self.dataset_value = args.value
        self.subsample = args.subsample

        # create widgets
        self.progress_bar = ProgressWidget(args, req)
        self.progress_bar.setTextVisible(True)
        self.progress_bar.setAlignment(QtCore.Qt.AlignCenter)
        self.addWidget(self.progress_bar, 0, 0)

        self.eta_display = QtWidgets.QLabel("Time Remaining:\t00:00.00")
        self.eta_display.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignCenter)
        font = QtGui.QFont()
        font.setPointSize(19)
        self.eta_display.setFont(font)
        self.addWidget(self.eta_display, 1, 0)

        # set up timer
        self.prev_pct = 0
        self.prev_time = 0
        self.mbar = 0
        self._idx_update = 0

    def data_changed(self, value, metadata, persist, mods):
        # update progress bar
        self.progress_bar.data_changed(value, metadata, persist, mods)

        # get value from dataset_db and update completion
        try:
            val = value[self.dataset_value]
        except (KeyError, ValueError, TypeError):
            val = 0

        # process statistics
        if val >= self.prev_pct:
            time_now = time()

            # update predicted completion time
            if (self._idx_update > 0) and (self._idx_update % self.subsample == 0):
                # calculate statistics
                d_mbar = (val - self.prev_pct) / (time_now - self.prev_time)
                self.mbar = self.mbar * (1. - 1. / (self._idx_update//self.subsample)) + d_mbar / (self._idx_update//self.subsample)
                # print("{}\t{}\t{:.5g}\t{:.5g}".format(self._idx_update, self._idx_update//self.subsample, d_mbar, self.mbar))

                # calculate remaining time
                time_remaining_s = (100. - val) / self.mbar
                # print("\n\tpct: {:2.2f}\t time_remaining_s: {:2.2f} => {:02d}:{:04.2f}\n".format(val, time_remaining_s, int(time_remaining_s//60), time_remaining_s % 60))
                self.eta_display.setText("Time Remaining:\t{:02d}:{:04.2f}".format(int(time_remaining_s // 60), time_remaining_s % 60))

            # update values
            self.prev_pct = val
            self.prev_time = time_now
            self._idx_update += 1

        # reset and start anew
        elif val < self.prev_pct:
            self.prev_pct = 0
            self.prev_time = 0
            self.mbar = 0
            self._idx_update = 0


def main():
    applet = SimpleApplet(CompletionMonitor)
    applet.add_dataset("value", "counter")
    applet.argparser.add_argument("--min", type=int, default=0,
                                  help="minimum (left) value of the bar")
    applet.argparser.add_argument("--max", type=int, default=100,
                                  help="maximum (right) value of the bar")
    applet.argparser.add_argument("--subsample", type=int, default=3,
                                  help="number of updates to buffer for completion time statistics")
    applet.run()


if __name__ == "__main__":
    main()
