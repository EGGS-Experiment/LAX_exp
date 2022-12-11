import numpy as np
from artiq.experiment import *
from LAX_exp.LAX.devices.beam_397_pump import Beam397Pump
from LAX_exp.LAX.devices.pmt import PMTCounter


class yzdeTest(EnvExperiment):
    """
    yzde test
    """

    def build(self):
        """
        Set devices and arguments for the experiment.
        """
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # experiment runs
        self.setattr_argument("repetitions",                    NumberValue(default=10, ndecimals=0, step=1, min=1, max=10000))

        # frequency scan
        self.setattr_argument("freq_qubit_scan_mhz",            Scannable(default=RangeScan(104.24, 104.96, 801),
                                                                    global_min=60, global_max=200, global_step=1,
                                                                    unit="MHz", scale=1, ndecimals=5), group = 'thkim.yzde')

        self.time_readout_mu = self.core.seconds_to_mu(10 * ms)


    def prepare(self):
        """
        Set up the dataset and prepare things such that
        the kernel functions have minimal overhead.
        """
        # devices
        self.pmt = PMTCounter(self)
        self.pump_397 = Beam397Pump(self)

        self.set_dataset("storage_tmp", [])
        self.setattr_dataset("storage_tmp")

        self.pmt._prepare_parameters()
        self.pmt.prepare_class()

        self.pump_397._prepare_parameters()
        self.pump_397.prepare_hardware()


    #@kernel(flags={"fast-math"})
    def run(self):
        """
        Run the experimental sequence.
        """
        #self.core.reset()

        yz = list(range(self.repetitions))
        for i in yz:
            print('\ta1: {}'.format(i))
            val = self.yzdetmp()
            self.append_to_dataset("storage_tmp", val)


    #@kernel
    @kernel(flags='fast-math')
    def yzdetmp(self):
        self.core.break_realtime()
        self.pump_397.cfg_sw(True)
        self.pmt.count(self.time_readout_mu)
        self.pump_397.cfg_sw(False)
        self.core.break_realtime()
        return self.pmt.fetch_count()

    def analyze(self):
        print('avg: {}'.format(np.mean(self.storage_tmp)))
