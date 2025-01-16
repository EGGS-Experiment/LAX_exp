import numpy as np
from artiq.experiment import *
from artiq.coredevice import ad9910

from LAX_exp.analysis import *
from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import InitializeQubit, RabiFlop, QubitPulseShape, Readout, RescueIon

from itertools import product


class BichromaticCalibration(LAXExperiment, Experiment):
    """
    Experiment: Bichromatic Calibration

    todo: document
    """
    name = 'Bichromatic Calibration'
    kernel_invariants = {
        'initialize_subsequence', 'pulseshape_subsequence', 'readout_subsequence', 'rescue_subsequence',
        'profile_target',
        'qubit_carrier', 'freq_qubit_carrier_default_ftw', 'ampl_qubit_carrier_default_asf',
        'att_qubit_carrier_default_mu',
        'ampl_729_carrier_asf', 'ampl_qubit_asf', 'att_729_carrier_mu', 'att_qubit_mu',
        'config_experiment_list'
    }

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions", NumberValue(default=20, precision=0, step=1, min=1, max=100000))

        # scan parameters - frequency
        self.setattr_argument("freq_qubit_mhz", NumberValue(default=101.3341, precision=6, step=1, min=50., max=400.),
                                                group="scan.frequency")
        self.setattr_argument("freq_729_carrier_center_mhz", NumberValue(default=80., precision=6, step=1, min=50., max=400.),
                                                                group="scan.frequency")
        self.setattr_argument("freq_729_carrier_sweep_khz_list",  Scannable(
                                                            default=[
                                                                CenterScan(0., 2000., 10, randomize=True),
                                                                ExplicitScan([0.]),
                                                                ExplicitScan([1303.29, -1303.29]),
                                                                RangeScan(-1303.29, 1303.29, 200, randomize=True),
                                                            ],
                                                            global_min=-20000., global_max=200000., global_step=10,
                                                            unit="kHz", scale=1, precision=3
                                                        ), group="scan.frequency")

        # scan parameters - time
        self.setattr_argument("equalize_delays",        BooleanValue(default=False), group="scan.time")
        self.setattr_argument("time_rabi_us_list",      Scannable(
                                                            default=[
                                                                ExplicitScan([6.05]),
                                                                RangeScan(1, 50, 200, randomize=True),
                                                                CenterScan(3.05, 5., 0.1, randomize=True),
                                                            ],
                                                            global_min=1, global_max=100000, global_step=1,
                                                            unit="us", scale=1, precision=5
                                                        ), group="scan.time")

        # carrier beam parameters
        self.setattr_argument("ampl_729_carrier_pct",   NumberValue(default=20, precision=3, step=5, min=0.01, max=50), group="beam.carrier")
        self.setattr_argument("att_729_carrier_db",     NumberValue(default=31.5, precision=1, step=0.5, min=13., max=31.5), group="beam.carrier")

        # beam parameters
        self.setattr_argument("enable_pulseshaping", BooleanValue(default=True), group="beam.qubit")
        self.setattr_argument("ampl_qubit_pct", NumberValue(default=20, precision=3, step=5, min=1, max=50), group="beam.qubit")
        self.setattr_argument("att_qubit_db",   NumberValue(default=31.5, precision=1, step=0.5, min=8, max=31.5), group="beam.qubit")

        # relevant devices
        self.setattr_device('qubit')

        # subsequences
        self.profile_target = 6
        self.initialize_subsequence =   InitializeQubit(self)
        self.pulseshape_subsequence =   QubitPulseShape(self, ram_profile=self.profile_target,
                                                        ampl_max_pct=self.ampl_qubit_pct, num_samples=1000,
                                                        pulse_shape="blackman")
        self.readout_subsequence =      Readout(self)
        self.rescue_subsequence =       RescueIon(self)

    def prepare_experiment(self):
        """
        Prepare & precompute experimental values.
        """
        '''
        CONVERT VALUES TO MACHINE UNITS
        '''
        # tmp remove
        self.qubit_carrier = self.get_device("urukul0_ch1")
        self.freq_qubit_carrier_default_ftw =   self.qubit_carrier.frequency_to_ftw(80. * MHz)
        self.ampl_qubit_carrier_default_asf =   self.qubit_carrier.amplitude_to_asf(50. / 100.)
        self.att_qubit_carrier_default_mu =     att_to_mu(13. * dB)
        # tmp remove

        # beam parameters
        self.ampl_729_carrier_asf = self.qubit.amplitude_to_asf(self.ampl_729_carrier_pct / 100.)
        self.att_729_carrier_mu =   att_to_mu(self.att_729_carrier_db * dB)

        self.ampl_qubit_asf =       self.qubit.amplitude_to_asf(self.ampl_qubit_pct / 100.)
        self.att_qubit_mu =         att_to_mu(self.att_qubit_db * dB)

        # create frequency config to equalize output at ion
        freq_729_config_ftw_list = np.array([
            [
                # note: carrier DDS controls SLS single-pass => 1.0x
                self.qubit.frequency_to_ftw(self.freq_729_carrier_center_mhz * MHz + freq_khz * kHz),
                # note: qubit DDS controls double-pass => 0.5x
                self.qubit.frequency_to_ftw(self.freq_qubit_mhz * MHz + 0.5 * (freq_khz * kHz))
            ]
            for freq_khz in self.freq_729_carrier_sweep_khz_list
        ])

        # convert time to machine units
        max_time_us = np.max(list(self.time_rabi_us_list))
        # create timing list such that all shots have same length
        time_rabiflop_mu_list = np.array([
            [self.core.seconds_to_mu((max_time_us - time_us) * us), self.core.seconds_to_mu(time_us * us)]
            for time_us in self.time_rabi_us_list
        ])
        # turn off delay equalization based on user selection
        if not self.equalize_delays: time_rabiflop_mu_list[:, 0] = np.int64(8)

        '''
        CREATE EXPERIMENT CONFIG
        '''
        # create an array of values for the experiment to sweep
        # note: need to use product since the constituent config lists are 2D
        self.config_experiment_list = np.array([
            np.concatenate((vals))
            for vals in product(freq_729_config_ftw_list, time_rabiflop_mu_list)
        ], dtype=np.int64)
        np.random.shuffle(self.config_experiment_list)

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_experiment_list),
                4)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self) -> TNone:
        self.core.break_realtime()

        # ensure DMA sequences use correct profile
        self.qubit.set_profile(0)
        # reduce attenuation/power of qubit beam to resolve lines
        self.qubit.set_att_mu(self.att_qubit_mu)
        self.qubit_carrier.set_att_mu(self.att_729_carrier_mu)
        self.core.break_realtime()

        # ensure qubit carrier is set up correctly on ALL profiles
        for i in range(8):
            self.qubit_carrier.set_mu(self.freq_qubit_carrier_default_ftw,
                                      asf=self.ampl_qubit_carrier_default_asf,
                                      profile=i)
            self.qubit_carrier.cpld.io_update.pulse_mu(8)
            delay_mu(5000)
        # enable RAM mode and clear DDS phase accumulator
        # self.qubit_carrier.write32(ad9910._AD9910_REG_CFR1,
        #                  # (1 << 16) |  # select_sine_output
        #                  (1 << 13) | # phase_autoclear
        #                  2          # default serial I/O configs
        #                  )
        self.qubit_carrier.cpld.io_update.pulse_mu(8)
        self.qubit_carrier.sw.on()
        self.core.break_realtime()

        # record subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.readout_subsequence.record_dma()
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        self.core.break_realtime()

        for trial_num in range(self.repetitions):
            self.core.break_realtime()

            # sweep exp config
            for config_vals in self.config_experiment_list:

                '''PREPARE & CONFIGURE'''
                # extract values from config list
                freq_729_carrier_ftw =  np.int32(config_vals[0])
                freq_qubit_ftw =        np.int32(config_vals[1])
                time_equalize_mu =      config_vals[2]
                time_pulse_mu =         config_vals[3]
                self.core.break_realtime()

                # set pulse time
                if self.enable_pulseshaping:
                    self.pulseshape_subsequence.configure(time_pulse_mu)

                # set qubit carrier frequency
                self.qubit_carrier.set_mu(freq_729_carrier_ftw, asf=self.ampl_729_carrier_asf,
                                          profile=self.profile_target)
                # set qubit frequency
                if self.enable_pulseshaping:
                    self.qubit.set_ftw(freq_qubit_ftw)
                else:
                    self.qubit.set_mu(freq_qubit_ftw, asf=self.ampl_qubit_asf, profile=self.profile_target)
                self.core.break_realtime()

                # # tmp remove
                # self.qubit_carrier.set_att_mu(self.att_729_carrier_mu)
                # # tmp remove

                '''INITIALIZE'''
                # initialize ion in S-1/2 state
                self.initialize_subsequence.run_dma()

                '''MAIN PULSE'''
                # add delay to ensure each shot takes same time
                delay_mu(time_equalize_mu)

                # rabi flop & read out
                if self.enable_pulseshaping:
                    self.pulseshape_subsequence.run()
                else:
                    self.qubit.cpld.set_profile(self.profile_target)
                    self.qubit.cpld.io_update.pulse_mu(8)
                    self.qubit.on()
                    delay_mu(time_pulse_mu)
                    self.qubit.off()

                '''READOUT & STORE RESULTS'''
                # read out
                self.readout_subsequence.run_dma()

                # update dataset
                self.update_results(freq_qubit_ftw,
                                    self.readout_subsequence.fetch_count(),
                                    freq_729_carrier_ftw,
                                    time_pulse_mu)
                self.core.break_realtime()

                # resuscitate ion
                self.rescue_subsequence.resuscitate()

            # rescue ion as needed
            self.rescue_subsequence.run(trial_num)

            # support graceful termination
            self.check_termination()
            self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def cleanup_experiment(self) -> TNone:
        """
        Clean up the experiment.
        """
        self.core.break_realtime()

        # set qubit carrier to default value (b/c AOM thermal drift) on ALL profiles
        for i in range(8):
            self.qubit_carrier.set_mu(self.freq_qubit_carrier_default_ftw,
                                      asf=self.ampl_qubit_carrier_default_asf,
                                      profile=i)
            self.qubit_carrier.cpld.io_update.pulse_mu(8)
            delay_mu(5000)
        self.qubit_carrier.sw.on()
        self.qubit_carrier.set_att_mu(self.att_qubit_carrier_default_mu)
        self.core.break_realtime()

    # ANALYSIS
    def analyze_experiment(self):
        pass
