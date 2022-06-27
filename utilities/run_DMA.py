from artiq.experiment import *


class DMARun(EnvExperiment):
    """
    DMA Run
    Run a DMA Exp by handle.
    """

    def build(self):
        self.setattr_device("core")
        self.setattr_device("core_dma")

    def prepare(self):
        pass

    @kernel
    def run(self):
        self.core.break_realtime()
        handle_tmp = self.core_dma.get_handle('PMT_exp')
        self.core.break_realtime()
        self.core_dma.playback_handle(handle_tmp)

    def analyze(self):
        th1 = self.get_dataset('pmt_test_dataset')
        print(th1)
