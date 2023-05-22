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

        # todo: set readout freq and ampl

        # scan parameters
        self.setattr_argument("freq_qubit_scan_mhz",                Scannable(
                                                                        default=CenterScan(103.7385, 0.02, 0.0005, randomize=True),
                                                                        global_min=60, global_max=200, global_step=1,
                                                                        unit="MHz", scale=1, ndecimals=5
                                                                    ), group=self.name)
        self.setattr_argument("time_qubit_us",                      NumberValue(default=3000, ndecimals=5, step=1, min=1, max=10000000), group=self.name)

        # get devices
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')
        self.setattr_device('repump_qubit')
        self.setattr_device('pmt')

        # get subsequences
        self.dopplercool_subsequence =                              DopplerCool(self)

    def prepare_experiment(self):
        # convert frequencies to machine units
        self.freq_qubit_scan_ftw =                                  np.array([hz_to_ftw(freq_mhz * MHz) for freq_mhz in self.freq_qubit_scan_mhz])

        # convert attenuation to machine units
        self.att_qubit_mu =                                         att_to_mu(self.att_qubit_db * dB)

        # get previous readout time to be used for precooling
        self.time_precool_mu =                                      self.get_parameter('time_readout_us', group='timing', override=True, conversion_function=seconds_to_mu, units=us)

        # get new actual readout time
        self.time_readout_mu =                                      self.core.seconds_to_mu(100 * us)

    @property
    def results_shape(self):
        return (self.repetitions * len(self.freq_qubit_scan_mhz),
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

        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def run_main(self):
        self.core.reset()

        # get custom readout handle
        _handle_precool_test = self.core_dma.get_handle('_PRECOOLTEST')
        self.core.break_realtime()

        for trial_num in range(self.repetitions):
            self.core.break_realtime()

            # sweep frequency
            for freq_ftw in self.freq_qubit_scan_ftw:

                # set new readout params in profile 3
                self.pump.set_mu(freq_ftw, asf=0x1FFF, profile=3)
                self.core.break_realtime()

                # run precool test
                self.core_dma.playback_handle('_PRECOOLTEST')

                # update dataset
                with parallel:
                    self.update_results(freq_ftw, self.pmt.fetch_count())
                    self.core.break_realtime()
