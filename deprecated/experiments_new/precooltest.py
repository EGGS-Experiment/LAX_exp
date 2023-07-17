import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import DopplerCool, Readout, RescueIon


class precooltest(LAXExperiment, Experiment):
    """
    Experiment: precooltest

    idk
    """
    name = 'precooltest'


    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions",                        NumberValue(default=30, ndecimals=0, step=1, min=1, max=10000))

        # precooling configuration
        self.setattr_argument("ampl_readout_list_pct",              Scannable(
                                                                        default=RangeScan(10, 50, 41, randomize=True),
                                                                        global_min=1, global_max=50, global_step=1,
                                                                        unit="pct", scale=1, ndecimals=2
                                                                    ), group='precool')
        self.setattr_argument("freq_readout_mhz",                   NumberValue(default=120, ndecimals=5, step=1, min=70, max=400), group='precool')
        self.setattr_argument("time_readout_us",                    NumberValue(default=100, ndecimals=3, step=10, min=1, max=10000000), group='precool')

        # get devices
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')
        self.setattr_device('repump_qubit')
        self.setattr_device('pmt')

        # get subsequences
        self.dopplercool_subsequence =                              DopplerCool(self)

    def prepare_experiment(self):
        # get previous readout time to be used for precooling
        self.time_precool_mu =                                      self.get_parameter('time_readout_us', group='timing', override=True, conversion_function=seconds_to_mu, units=us)

        # convert readout config to machine units
        self.ampl_readout_list_asf =                                np.array([pct_to_asf(ampl_pct) for ampl_pct in list(self.ampl_readout_list_pct)])
        self.freq_readout_ftw =                                     hz_to_ftw(self.freq_readout_mhz * MHz)
        self.time_readout_mu =                                      self.core.seconds_to_mu(self.time_readout_us * us)

    @property
    def results_shape(self):
        return (self.repetitions * len(self.ampl_readout_list_pct),
                2)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self):
        self.core.break_realtime()

        # record custom precooling test sequence
        with self.core_dma.record('_PRECOOLTEST'):
            # ensure repump beams are turned on correctly
            self.repump_cooling.on()
            self.repump_qubit.on()

            # precool @ readout params
            self.pump.readout()
            self.pump.on()
            delay_mu(self.time_precool_mu)
            self.pump.off()

            # doppler cooling
            self.dopplercool_subsequence.run()

            # new readout
            self.pump.set_profile(3)
            self.pump.on()
            self.pmt.count(self.time_readout_mu)
            self.pump.off()

        # self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def run_main(self):
        self.core.reset()

        # get custom readout handle
        _handle_precool_test = self.core_dma.get_handle('_PRECOOLTEST')
        self.core.break_realtime()

        for trial_num in range(self.repetitions):
            self.core.break_realtime()

            # sweep frequency
            for ampl_asf in self.ampl_readout_list_asf:

                # set new readout params in profile 3
                self.pump.set_mu(self.freq_readout_ftw, asf=ampl_asf, profile=3)

                # run precool test
                self.core_dma.playback_handle(_handle_precool_test)

                # update dataset
                with parallel:
                    self.update_results(self.pump.asf_to_amplitude(ampl_asf) * 100, self.pmt.fetch_count())
                    self.core.break_realtime()
