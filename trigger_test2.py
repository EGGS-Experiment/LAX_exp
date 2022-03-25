from artiq.experiment import *


class trigger_test2(EnvExperiment):
    """
    testing trigger2
    """

    def build(self):
        # get core devices
        self.setattr_device('core')
        self.setattr_device('core_dma')

    def run(self):
        pass
    #
    # @kernel
    # def run(self):
    #     self.core.reset()
    #     handle = self.core_dma.get_handle('PMT_exp')
    #     self.core.reset()
    #     self.core_dma.playback_handle(handle)

    def analyze(self):
        th1 = self.get_dataset('ttt')
        th1 = list(filter(None, th1))
        print(th1)
