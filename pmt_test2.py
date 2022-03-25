from artiq.experiment import *


class pmt_test2(EnvExperiment):
    """
    testing pmt2
    """

    def build(self):
        # get core devices
        self.setattr_device('core')
        self.setattr_device('core_dma')

    @kernel
    def run(self):
        self.core.reset()
        handle = self.core_dma.get_handle('PMT_exp')
        self.core.reset()
        self.core_dma.playback_handle(handle)

    def analyze(self):
        print(self.get_dataset('pmt_test_dataset'))
