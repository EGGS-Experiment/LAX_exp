from artiq.experiment import *
from artiq.coredevice.ad9910 import _AD9910_REG_FTW

from numpy import zeros, ones
from datetime import datetime


class pmt_test(EnvExperiment):
    """
    testing pmt
    """

    def build(self):
        # get arguments
        self.setattr_argument('ttl_number', NumberValue(ndecimals=0, step=1, min=0, max=23))
        self.setattr_argument('trigger_ttl_number', NumberValue(default=-1, ndecimals=0, step=1, min=-1, max=23))
        self.setattr_argument('trigger_active', BooleanValue(default=False))
        self.setattr_argument('bin_time_us', NumberValue(default=1, ndecimals=0, step=1, min=1, max=100))
        self.setattr_argument('reset_time_us', NumberValue(default=1, ndecimals=0, step=1, min=1, max=100))
        self.setattr_argument('length_us', NumberValue(default=100, ndecimals=0, step=1, min=1, max=100000))
        self.setattr_argument('edge_method', StringValue(default='rising'))
        # get core devices
        self.setattr_device('core')
        self.setattr_device('core_dma')

    def prepare(self):
        # get ttl input
        pmt_name = 'ttl{:d}'.format(self.ttl_number)
        self.ttl_input = self.get_device(pmt_name)
        # get trigger input
        trigger_name = 'ttl{:d}'.format(self.trigger_ttl_number)
        if self.trigger_active:
            self.trigger_input = self.get_device(trigger_name)
        else:
            self.trigger_input = self.ttl_input
        # check TTLs are input
        if (self.ttl_input.__class__.__name__ != 'TTLInOut') or (self.trigger_input.__class__.__name__ != 'TTLInOut'):
            raise Exception('Error: TTLs are output-only.')
        # check edge gating input valid
        if self.edge_method not in ('rising', 'falling', 'both'):
            raise Exception('Error: invalid edge method.')
        # get edge gating method
        if self.edge_method == 'rising':
            self.gate_edge = getattr(self.ttl_input, 'gate_rising_mu')
        elif self.edge_method == 'falling':
            self.gate_edge = getattr(self.ttl_input, 'gate_falling_mu')
        elif self.edge_method == 'both':
            self.gate_edge = getattr(self.ttl_input, 'gate_both_mu')
        # process input values into mu
        self.bin_time_mu = self.core.seconds_to_mu(self.bin_time_us * us)
        self.reset_time_mu = self.core.seconds_to_mu(self.reset_time_us * us)
        self.length_mu = self.core.seconds_to_mu(self.length_us * us)
        # composite values
        self.total_edges = 0
        self.num_bins = int(self.length_us / (self.bin_time_us + self.reset_time_us))
        # create dataset
        date = datetime.now()
        self.dataset_name = 'pmt_{:s}_{:02d}_{:02d}_{:02d}{:02d}'.format(str(date.year), date.month, date.day, date.hour, date.minute)
        self.dataset_name = 'pmt_test_dataset' #todo: undo
        self.set_dataset(self.dataset_name, zeros(self.num_bins), broadcast=True, persist=True)
        self.setattr_dataset(self.dataset_name)
        #todo: undo print values
        print(self.get_dataset('pmt_test_dataset'))
        print('bin time (us):', self.bin_time_us)
        print('reset time (us):', self.reset_time_us)
        print('length (us):', self.length_us)
        print('bins:', self.num_bins)

    @kernel
    def run(self):
        self.core.reset()
        self.ttl_input.input()
        self.trigger_input.input()
        self.core.break_realtime()
        print('scde')
        self.core.reset()
        with self.core_dma.record('pt'):
            for i in range(self.num_bins):
                # edges = self.ttl_input.count(self.gate_edge(self.bin_time_mu))
                # self.total_edges += edges
                # self.mutate_dataset(self.dataset_name, i, edges)
                self.mutate_dataset(self.dataset_name, i, self.ttl_input.count(self.gate_edge(self.bin_time_mu)))
                delay_mu(self.reset_time_mu)
        self.core.break_realtime()

    def analyze(self):
        # print('total edges:', self.total_edges)
        pass
