import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import InitializeQubit, RabiFlop, Readout, RescueIon
from LAX_exp.system.subsequences import NoOperation, SidebandCoolContinuous, SidebandCoolPulsed


class RabiFlopping(LAXExperiment, Experiment):
    """
    Experiment: Rabi Flopping

    Measures ion fluorescence vs 729nm pulse time and frequency.
    """
    name = 'Rabi Flopping'


    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions",                        NumberValue(default=30, ndecimals=0, step=1, min=1, max=10000))

        # rabi flopping arguments
        self.setattr_argument("cooling_type",                       EnumerationValue(["Doppler", "SBC - Continuous", "SBC - Pulsed"], default="Doppler"))
        self.setattr_argument("time_rabi_us_list",                  Scannable(
                                                                        default=RangeScan(0, 100, 101, randomize=True),
                                                                        global_min=1, global_max=100000, global_step=1,
                                                                        unit="us", scale=1, ndecimals=5
                                                                    ), group=self.name)
        self.setattr_argument("freq_rabiflop_mhz",                  NumberValue(default=104.06, ndecimals=5, step=1, min=1, max=10000), group=self.name)
        self.setattr_argument("att_readout_db",                     NumberValue(default=8, ndecimals=1, step=0.5, min=8, max=31.5), group=self.name)


        # get devices
        self.setattr_device('qubit')

        # prepare sequences
        self.initialize_subsequence =                               InitializeQubit(self)
        self.doppler_subsequence =                                  NoOperation(self)
        self.sidebandcool_pulsed_subsequence =                      SidebandCoolPulsed(self)
        self.sidebandcool_continuous_subsequence =                  SidebandCoolContinuous(self)
        self.readout_subsequence =                                  Readout(self)
        self.rescue_subsequence =                                   RescueIon(self)

    def prepare_experiment(self):
        # choose correct cooling subsequence
        if self.cooling_type == "Doppler":
            self.cooling_subsequence =                              self.doppler_subsequence
        elif self.cooling_type == "SBC - Continuous":
            self.cooling_subsequence =                              self.sidebandcool_continuous_subsequence
        elif self.cooling_type == "SBC - Pulsed":
            self.cooling_subsequence =                              self.sidebandcool_pulsed_subsequence


        # convert rabi flopping arguments
        max_time_us =                                               np.max(list(self.time_rabi_us_list))
        self.time_rabiflop_mu_list =                                np.array([
                                                                        [seconds_to_mu((max_time_us - time_us) * us), seconds_to_mu(time_us * us)]
                                                                        for time_us in self.time_rabi_us_list
                                                                    ])
        self.freq_rabiflop_ftw =                                    hz_to_ftw(self.freq_rabiflop_mhz * MHz)
        self.att_readout_mu =                                       att_to_mu(self.att_readout_db * dB)

    @property
    def results_shape(self):
        return (self.repetitions * len(self.time_rabiflop_mu_list),
                2)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self):
        self.core.break_realtime()

        # record subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.cooling_subsequence.record_dma()
        self.readout_subsequence.record_dma()

        # set qubit readout waveform
        self.qubit.set_mu(self.freq_rabiflop_ftw, asf=self.qubit.ampl_qubit_asf, profile=0)

    @kernel(flags={"fast-math"})
    def run_main(self):
        self.core.reset()

        for trial_num in range(self.repetitions):

            # sweep time
            for time_rabi_pair_mu in self.time_rabiflop_mu_list:
                self.core.break_realtime()

                # initialize ion in S-1/2 state
                self.initialize_subsequence.run_dma()

                # run sideband cooling
                self.cooling_subsequence.run_dma()

                # prepare qubit beam for readout
                self.qubit.set_profile(0)
                self.qubit.set_att_mu(self.att_readout_mu)
                delay_mu(time_rabi_pair_mu[0])

                # rabi flop
                self.qubit.on()
                delay_mu(time_rabi_pair_mu[1])
                self.qubit.off()

                # do readout
                self.readout_subsequence.run_dma()

                # update dataset
                with parallel:
                    self.update_results(time_rabi_pair_mu[1], self.readout_subsequence.fetch_count())
                    self.core.break_realtime()
