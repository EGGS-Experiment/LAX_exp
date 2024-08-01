import labrad
import numpy as np
from os import environ
from artiq.experiment import *
from datetime import datetime
from artiq.coredevice.urukul import urukul_sta_rf_sw, SPI_CONFIG
from artiq.coredevice.rtio import (rtio_output, rtio_input_timestamp,
                                   rtio_input_data)

from LAX_exp.base import LAXExperiment, LAXEnvironment
from LAX_exp.system.objects.PhaserPulseShaper import PhaserPulseShaper
from LAX_exp.system.objects.SpinEchoWizard import SpinEchoWizard

import matplotlib.pyplot as plt


class SpinEchoWizardTest(LAXEnvironment, Experiment):
    """
    SpinEchoWizardTest
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
        self.pulse_shaper = PhaserPulseShaper(self)

    def prepare(self):
        pass

    def run(self):
        self.spinecho_wizard.prepare()
        self.spinecho_wizard.calculate_pulseshape()
        self.spinecho_wizard.compile_waveform()

        # print(self.spinecho_wizard._ampl_tmp_arr)
        # print('\tampl size: {}\n'.format(np.shape(self.spinecho_wizard._ampl_tmp_arr)))
        # # print(self.spinecho_wizard._phas_tmp_arr)
        # print('\tphas size: {}\n'.format(np.shape(self.spinecho_wizard._phas_tmp_arr)))
        # # print(self.spinecho_wizard._time_tmp_arr)
        # print('\ttime size: {}\n'.format(np.shape(self.spinecho_wizard._time_tmp_arr)))

        # plt.plot(self.spinecho_wizard._ampl_tmp_arr[0])
        # plt.show()
        self.spinecho_wizard.display_waveform()

