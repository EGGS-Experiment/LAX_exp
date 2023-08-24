import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment

# todo: sort out labrad imports


class PMTAlignment(LAXExperiment, Experiment):
    """
    Utility: PMT Alignment

    Read PMT counts over time with cooling repump on/off to compare signal/background.
    """
    name = 'PMT Alignment'


    def build_experiment(self):
        """
        Set devices and arguments for the experiment.
        """
        # timing
        self.setattr_argument('time_total_s',                       NumberValue(default=10000, ndecimals=6, step=1, min=0, max=100000), group='timing')
        self.setattr_argument("samples_per_point",                  NumberValue(default=50, ndecimals=0, step=10, min=1, max=500), group='timing')

        # relevant devices
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')
        self.setattr_device('repump_qubit')
        self.setattr_device('pmt')

        # todo: idk
        self.setattr_device('urukul2_cpld')

    def prepare_experiment(self):
        # get relevant timings and calculate the number of repetitions
        self.time_readout_mu =                                      self.get_parameter('time_readout_us', group='timing', override=True, conversion_function=seconds_to_mu, units=us)
        self.repetitions =                                          round(self.core.seconds_to_mu(self.time_total_s) / (2 * self.time_readout_mu))

        # todo: create holder variable that can store averaged counts
        # todo: store absolute time variables
        

    @property
    def results_shape(self):
        return (self.repetitions, 3)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self):
        self.core.break_realtime()

        # ensure beams are prepared for DMA
        self.pump.readout()

        # record custom sequence
        with self.core_dma.record('_PMT_ALIGNMENT'):
            # set readout beams and
            self.urukul2_cpld.cfg_switches(0b1110)
            self.pmt.count(self.time_readout_mu)
            self.pump.off()
            # self.

        # tmp remove
        # set profile
        # record DMA
        # for loop
        # 0b1110
        # delay readout time
        # 0b0010
        # delay readout time

        # note: bulk of this code is copy

        # set readout waveform
        self.pump.readout()

        # readout pulse
        self.pump.on()
        self.repump_cooling.on()
        self.pmt.count(self.time_readout_mu)
        self.pump.off()
        
        
    @kernel(flags={"fast-math"})
    def run_main(self):
        self.core.reset()

        # retrieve pulse sequence handle
        handle = self.core_dma.get_handle(_DMA_HANDLE)
        self.core.break_realtime()


        # MAIN LOOP
        for i in self.repetitions:

            # run pulse sequence from core DMA
            self.core_dma.playback_handle(handle)

            # record pmt counts to dataset
            with parallel:
                self.update_dataset(i, self.pmt_counter.fetch_count(), self.pmt_counter.fetch_count())
                delay_mu(self.time_reset_mu)


    @rpc(flags={"async"})
    def updateLabrad(self, index, counts_signal, counts_background):
        # create update array todo: document better
        update_arr = np.array([counts_signal, counts_background, counts_signal - counts_background])

        # update dataset
        self.mutate_dataset('pmt_dataset', index, update_arr)

        # todo: if index % update_iter, update labrad


    # ANALYZE
    def analyze(self):
        """
        Analyze the results from the experiment.
        """
        # print results
        print('\tCounts: {:.3f} +/- {:.3f}'.format(mean(self.pmt_dataset), std(self.pmt_dataset)))
