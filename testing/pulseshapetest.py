import labrad
import numpy as np
from os import environ
from artiq.experiment import *
from datetime import datetime
from artiq.coredevice.urukul import urukul_sta_rf_sw, SPI_CONFIG
from artiq.coredevice.rtio import (rtio_output, rtio_input_timestamp,
                                   rtio_input_data)

from LAX_exp.base import LAXExperiment
from LAX_exp.system.objects.PhaserPulseShape import PhaserPulseShape
from LAX_exp.system.objects.SpinEchoWizard import SpinEchoWizard


class pulseshapingtest(LAXExperiment, Experiment):
    """
    pulseshapingtest
    """

    def build(self):
        # core devices
        self.setattr_device("core")
        self.setattr_device("core_dma")
        self.setattr_device("ttl8")
        self.setattr_device("ttl9")

        self.setattr_device("ttl0_counter")
        self.setattr_device('urukul2_ch2')
        self.setattr_device('urukul2_cpld')

        self.setattr_device('phaser0')
        self.setattr_device('phaser_eggs')

        # objects
        self.spinecho_wizard = SpinEchoWizard(self)
        self.pulse_shaper = PhaserPulseShape(self)

    def prepare(self):
        pass

    def run(self):
        self.pulse_shaper.prepare()
        self.pulse_shaper.calculate_pulseshape()
        self.pulse_shaper.compile_waveform()

        print(self.pulse_shaper._ampl_tmp_arr)
        print('\tsize: {}\n'.format(len(self.pulse_shaper._ampl_tmp_arr)))
        print(self.pulse_shaper._phas_tmp_arr)
        print('\tsize: {}\n'.format(len(self.pulse_shaper._phas_tmp_arr)))
        print(self.pulse_shaper._time_tmp_arr)
        print('\tsize: {}\n'.format(len(self.pulse_shaper._time_tmp_arr)))

