import numpy as np
from random import shuffle
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import InitializeQubit, SidebandCoolContinuous, SidebandCoolPulsed, RabiFlop, Readout, RescueIon


class SidebandCooling(LAXExperiment, Experiment):
    """
    Experiment: Sideband Cooling

    Measures temperature after a given number of RSB pulses.
    """
    name = 'Sideband Cooling'


    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions",                            NumberValue(default=1, ndecimals=0, step=1, min=1, max=10000))

        # sideband cooling type
        self.setattr_argument("cooling_type",                           EnumerationValue(["Continuous", "Pulsed"], default="Continuous"))


        # sideband cooling readout
        self.setattr_argument("freq_rsb_scan_mhz",                      Scannable(
                                                                            # default=CenterScan(103.078, 0.02, 0.0005),
                                                                            default=(ExplicitScan([103.078])),
                                                                            global_min=30, global_max=200, global_step=1,
                                                                            unit="MHz", scale=1, ndecimals=5
                                                                        ), group='sideband_readout')
        self.setattr_argument("freq_bsb_scan_mhz",                      Scannable(
                                                                            # default=CenterScan(103.783, 0.02, 0.0005),
                                                                            default=(ExplicitScan([103.783])),
                                                                            global_min=30, global_max=200, global_step=1,
                                                                            unit="MHz", scale=1, ndecimals=5
                                                                        ), group='sideband_readout')
        self.setattr_argument("time_readout_pipulse_us",                NumberValue(default=90, ndecimals=5, step=1, min=1, max=10000), group='sideband_readout')
        self.setattr_argument("ampl_readout_pipulse_pct",               NumberValue(default=50, ndecimals=5, step=1, min=1, max=100), group='sideband_readout')
        self.setattr_argument("att_readout_db",                         NumberValue(default=8, ndecimals=1, step=0.5, min=8, max=31.5), group='sideband_readout')

        # get relevant devices
        self.setattr_device('qubit')

        # get subsequences
        self.initialize_subsequence =                                   InitializeQubit(self)
        self.sidebandcool_pulsed_subsequence =                          SidebandCoolPulsed(self)
        self.sidebandcool_continuous_subsequence =                      SidebandCoolContinuous(self)
        self.rabiflop_subsequence =                                     RabiFlop(self, time_rabiflop_us=self.time_readout_pipulse_us)
        self.readout_subsequence =                                      Readout(self)
        self.rescue_subsequence =                                       RescueIon(self)

    def prepare_experiment(self):
        # choose correct cooling subsequence
        if self.cooling_type == "Continuous":
            self.sidebandcool_subsequence =                             self.sidebandcool_continuous_subsequence
        elif self.cooling_type == "Pulsed":
            self.sidebandcool_subsequence =                             self.sidebandcool_pulsed_subsequence

        # convert readout frequencies to machine units
        self.freq_readout_ftw_list =                                    np.array([self.qubit.frequency_to_ftw(freq_mhz * MHz)
                                                                        for freq_mhz in (list(self.freq_rsb_scan_mhz) + list(self.freq_bsb_scan_mhz))])
        # combine & shuffle readout frequencies
        shuffle(self.freq_readout_ftw_list)

        # convert readout parameters
        self.time_readout_pipulse_mu =                                  self.core.seconds_to_mu(self.time_readout_pipulse_us * us)
        self.ampl_readout_pipulse_asf =                                 self.qubit.amplitude_to_asf(self.ampl_readout_pipulse_pct / 100)

        # convert attenuation to machine units
        self.att_readout_mu =                                           att_to_mu(self.att_readout_db * dB)

    @property
    def results_shape(self):
        return (self.repetitions * len(self.freq_readout_ftw_list),
                2)


    # MAIN SEQUENCE
    @kernel
    def initialize_experiment(self):
        self.core.break_realtime()

        # record subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.sidebandcool_subsequence.record_dma()

        # record custom readout sequence
        # note: this is necessary since DMA sequences will preserve urukul attenuation register
        with self.core_dma.record('_SBC_READOUT'):
            # set readout waveform for qubit
            self.qubit.set_profile(0)
            self.qubit.set_att_mu(self.att_readout_mu)

            # transfer population to D-5/2 state
            self.rabiflop_subsequence.run()

            # read out fluorescence
            self.readout_subsequence.run()

    @kernel
    def run_main(self):
        self.core.reset()

        # get custom readout handle
        _handle_sbc_readout = self.core_dma.get_handle('_SBC_READOUT')
        self.core.break_realtime()

        for trial_num in range(self.repetitions):

            # sweep frequency
            for freq_ftw in self.freq_readout_ftw_list:

                # set frequency
                self.qubit.set_mu(freq_ftw, asf=self.ampl_readout_pipulse_asf, profile=0)
                self.core.break_realtime()

                # initialize ion in S-1/2 state
                self.initialize_subsequence.run_dma()

                # sideband cool
                self.sidebandcool_subsequence.run_dma()

                # custom SBC readout
                self.core_dma.playback_handle(_handle_sbc_readout)

                # update dataset
                with parallel:
                    self.update_results(freq_ftw, self.readout_subsequence.fetch_count())
                    self.core.break_realtime()

            # rescue ion as needed
            self.rescue_subsequence.run(trial_num)
