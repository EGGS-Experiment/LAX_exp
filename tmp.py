import numpy as np
import labrad
from os import environ
from artiq.experiment import *
from artiq.coredevice.urukul import urukul_sta_rf_sw, SPI_CONFIG
from artiq.coredevice.rtio import (rtio_output, rtio_input_timestamp,
                                   rtio_input_data)

class testarg12(EnvExperiment):
    """
    testarg34
    Testing
    """

    def build(self):
        self.setattr_device("core")
        self.setattr_device("ttl8")

        self.setattr_argument("freq_eggs_heating_secular_mhz",          NumberValue(default=1.6, ndecimals=5, step=0.1, min=0.001, max=1000000))
        self.setattr_argument("freq_eggs_heating_mhz_list",             Scannable(
                                                                            default=CenterScan(85, 5, 0.2, randomize=True),
                                                                            global_min=30, global_max=400, global_step=1,
                                                                            unit="MHz", scale=1, ndecimals=5
                                                                        ))

        # self.setattr_device("core_dma")
        # self.setattr_device('urukul1_ch0')

        # self.set_dataset('ampl_qubit_pct', 50.0, broadcast=True, persist=True)
        # self.set_dataset('ampl_repump_cooling_pct', 10.0, broadcast=True, persist=True)

    def prepare(self):
        self.set_dataset('freq_probe_shutter_mhz',              210.0, broadcast=True,  persist=True)
        self.set_dataset('freq_pump_shutter_mhz',               220.0, broadcast=True,  persist=True)
        self.set_dataset('freq_repump_cooling_shutter_mhz',     230.0, broadcast=True,  persist=True)
        self.set_dataset('freq_repump_qubit_shutter_mhz',       240.0, broadcast=True,  persist=True)
        self.set_dataset('freq_qubit_shutter_mhz',              250.0, broadcast=True,  persist=True)

        self.set_dataset('ampl_shutter_pct',                    5.0, broadcast=True,    persist=True)


        # f_tmp = self.freq_eggs_heating_secular_mhz
        # freq_eggs_sidebands_mhz =                                       np.array([[freq_mhz - f_tmp, freq_mhz + f_tmp] for freq_mhz in self.freq_eggs_heating_mhz_list])
        #
        # from scipy.interpolate import Akima1DInterpolator
        # yz1 = self.get_dataset('calibration.eggs.resonance_ratio_curve_mhz')
        # spl1 = Akima1DInterpolator(yz1[:, 0], yz1[:, 1])
        #
        # # look up values
        # th0de = np.zeros((len(freq_eggs_sidebands_mhz), 2), dtype=float)
        # for i, val in enumerate(freq_eggs_sidebands_mhz):
        #     a1, a2 = (spl1(val[0]), spl1(val[1]))
        #     th0de[i] = np.array([a2/(a1+a2), a1/(a1+a2)])

        # self.call_child_method('prepare')
        # self.hasenv = testhasenv(self, {'aa':11,'bbb':22}, 9031)
        # self.cxn = labrad.connect(environ['LABRADHOST'], port=7682, tls_mode='off', username='', password='lab')
        # self.os = self.cxn.oscilloscope_server
        # print(self.os.list_devices())

    #@kernel
    def run(self):
        pass
        #self.set_dataset('pmt.input_channe', 'rising', broadcast=True, persist=True)

    def analyze(self):
        print('test done')
