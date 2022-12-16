import numpy as np
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
        self.setattr_device("core_dma")
        self.setattr_argument("test1", PYONValue({'a':1,'b':2}))
        self.setattr_device('urukul1_ch0')

        self.set_dataset('timing.time_rabiflop_us', 15, broadcast=True, persist=True)
        # self.set_dataset('ampl_qubit_pct', 50.0, broadcast=True, persist=True)
        # self.set_dataset('ampl_repump_cooling_pct', 10.0, broadcast=True, persist=True)

    # def prepare(self):
    #     val = self.get_dataset('timing.time_profileswitch_delay_us')
    #     print('val: {}'.format(self.core.seconds_to_mu(us*val)))


    #@kernel
    def run(self):
        #pass
        # self.core.reset()
        # self.core.break_realtime()
        # self.urukul1_ch0.init()
        # delay(10*ms)
        management_parameters = {
            # experiments
            "management.completition_pct":                              0.00,

            # datasets
            "management.dataset_save_locations":                        ["Z:\\Motion\\Data"],
        }
        for parameter_name, parameter_value in management_parameters.items():
            self.set_dataset(parameter_name, parameter_value, broadcast=True, persist=True)

        #self.set_dataset('dds.dds_board_tickle_num', 0, broadcast=True, persist=True)
        #self.set_dataset('dds.dds_tickle_channel', 3, broadcast=True, persist=True)
        #self.set_dataset('pmt_gating_edge', 'rising', broadcast=True, persist=True)
        #print('type: {}'.format(type(self.test1)))
        #print(self.test1)

