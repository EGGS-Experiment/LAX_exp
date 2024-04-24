#!/usr/bin/env python3

from PyQt5 import QtWidgets
from artiq.applets.simple import SimpleApplet

from playsound import playsound


class TextWidget(QtWidgets.QLabel):
    def __init__(self, args):
        QtWidgets.QLabel.__init__(self)
        self.dataset_name = args.dataset
        self.setStyleSheet('background-color: green; color: white')

    def data_changed(self, data, mods):
        try:
            # tmp remove
            print(data[self.dataset_name])

            # get status text
            text = data[self.dataset_name][1]

            # configure error display
            if 'ERROR' in text.split(':')[0].upper():
                self.setStyleSheet('background-color: red; color: white')

                # todo: play sound

            # configure display for normal operation
            else:
                self.setStyleSheet('background-color: green; color: white')

        except (KeyError, ValueError, TypeError):
            n = "UNKNOWN"

        # set display text
        self.setText(text)


def main():
    applet = SimpleApplet(TextWidget)
    applet.add_dataset("dataset", "dataset to show")
    applet.run()

if __name__ == "__main__":
    main()
