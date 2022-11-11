import numpy as np
from artiq.experiment import *
_DMA_HANDLE_READOUT = "thkim"


class TTLSchmitt(EnvExperiment):
    """
    TTL Schmitt Trigger
    Correlate photon counts with the phase of an RF signal.
    """

    def build(self):
        """
        Set devices and arguments for the experiment.
        """
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # tmp remove
        self.pmt_input_channel = 0
        self.rf_input_channel = 3

        # devices
        self.pmt_counter =                                      self.get_device("ttl_counter{:d}".format(self.pmt_input_channel))
        self.rf_trigger =                                       self.get_device("ttl_counter{:d}".format(self.rf_input_channel))

        # values
        self.setattr_argument("repetitions",                    NumberValue(default=100, ndecimals=0, step=1, min=1, max=10000))
        self.setattr_argument("time_count_us",                  NumberValue(default=1000, ndecimals=5, step=1, min=1, max=10000000))


    def prepare(self):
        """
        Set up the dataset and prepare things such that the kernel functions have minimal overhead.
        """
        # timing
        self.time_count_mu =                                    self.core.seconds_to_mu(self.time_count_us * us)

        # dataset
        self.set_dataset("ttl_trigger_test", [])
        self.setattr_dataset("ttl_trigger_test")


    @kernel
    def run(self):
        self.core.break_realtime()

        # readout sequence
        with self.core_dma.record(_DMA_HANDLE_READOUT):
            self.rf_trigger.gate_falling_mu(self.time_count_mu)

        # get handles
        handle_readout = self.core_dma.get_handle(_DMA_HANDLE_READOUT)
        self.core.break_realtime()

        # main loop
        for i in range(self.repetitions):

            # read out
            self.core_dma.playback_handle(handle_readout)

            # record data
            self.append_to_dataset("ttl_trigger_test", self.rf_trigger.fetch_count())
            self.core.break_realtime()


    def analyze(self):
        """
        Analyze the results from the experiment.
        """
        # turn dataset into numpy array for ease of use
        self.ttl_trigger_test = np.array(self.ttl_trigger_test)
        print(self.ttl_trigger_test)

        # print results
        print("counts: {:.4f} +/- {:.4f}".format(np.mean(self.ttl_trigger_test), np.std(self.ttl_trigger_test)))
