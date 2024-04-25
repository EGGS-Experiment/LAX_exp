#!/usr/bin/env python3

from PyQt5 import QtWidgets
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QFont

from artiq.applets.simple import SimpleApplet

from playsound import playsound

# _PLAYSOUND_ENABLE = False
# try:
#     # todo: try to make resource imports programmatic
#     #from inspect import getfile, getmodule
#     from playsound import playsound
#     from twisted.internet.threads import deferToThread
#     _UNLOCKED_SOUND_PATH = os.path.join(os.environ['EGGS_LABRAD_ROOT'], 'EGGS_labrad\\clients\\wavemeter_client\\channel_unlocked_short.mp3')
#     _PLAYSOUND_ENABLE = True
# except Exception as e:
#     print('Error: playsound package not installed.')


class TextWidget(QtWidgets.QLabel):
    def __init__(self, args):
        QtWidgets.QLabel.__init__(self)
        self.dataset_name = args.dataset

        # configure display
        self.setStyleSheet('background-color: green; color: white')
        self.setFont(QFont('MS Shell Dlg 2', pointSize=19))
        self.setAlignment(Qt.AlignCenter | Qt.AlignHCenter)


    def data_changed(self, data, mods):
        try:
            # get status text
            text = data[self.dataset_name][1]

            # configure error display
            if 'ERR' in text.split(':')[0].upper():
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
