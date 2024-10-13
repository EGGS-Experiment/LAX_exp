#!/usr/bin/env python3

import os
from PyQt5 import QtWidgets
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QFont

from artiq.applets.simple import SimpleApplet
# todo: clean up playsound structure stuff


'''SET UP AUDIO WARNING'''
_SOUND_PATH_LOCAL = "EGGS_labrad\\clients\\wavemeter_client\\channel_unlocked_short.mp3"
_PLAYSOUND_SUPPORTED = False
def playsound_warning(reps):
    return

try:
    from playsound import playsound
    _SOUND_PATH_WARNING = os.path.join(os.environ['EGGS_LABRAD_ROOT'], _SOUND_PATH_LOCAL)
    _PLAYSOUND_SUPPORTED = True
    def playsound_warning(reps):
        for i in range(reps):
            playsound(_SOUND_PATH_WARNING)
except Exception as e:
    print('Warning: playsound package not installed.')


class TextWidget(QtWidgets.QLabel):
    def __init__(self, args, req):
        QtWidgets.QLabel.__init__(self)
        self.dataset_name = args.dataset
        self.req = req

        # configure display
        self.setStyleSheet('background-color: green; color: white')
        self.setFont(QFont('MS Shell Dlg 2', pointSize=19))
        self.setAlignment(Qt.AlignCenter | Qt.AlignHCenter)

        # configure audible warnings
        self.play_sound = args.playsound
        self.sound_reps = 3

        # play sound to test
        if self.play_sound:
            playsound_warning(self.sound_reps)

    def data_changed(self, value, metadata, persist, mods):
        try:
            # get status text
            text = value[self.dataset_name]

            # configure display - error status
            if 'ERR' in text.split(':')[0].upper():
                self.setStyleSheet('background-color: red; color: white')
                # play sound to signal error
                if self.play_sound:
                    playsound_warning(self.sound_reps)

            # configure display - normal operation
            else:
                self.setStyleSheet('background-color: green; color: white')

        except (KeyError, ValueError, TypeError):
            text = "UNKNOWN"
            self.setStyleSheet('background-color: black; color: white')

        # set display text
        self.setText(text)


def main():
    applet = SimpleApplet(TextWidget)
    applet.argparser.add_argument("--playsound", default=False, type=bool,
                                  help="Play an audible warning upon a error condition.")
    applet.add_dataset("dataset", "dataset to show")
    applet.run()


if __name__ == "__main__":
    main()
