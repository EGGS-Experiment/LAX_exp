from artiq.experiment import *



class dma_test1(EnvExperiment):
    """
    testing dma1
    """

    def build(self):
        # get core devices
        self.setattr_device('core')
        self.setattr_device('core_dma')
        self.setattr_device('ttl5')
        self.setattr_device('ttl6')

    @kernel
    def run(self):
        self.core.reset()
        with self.core_dma.record('th1'):
            self.ttl5.on()

