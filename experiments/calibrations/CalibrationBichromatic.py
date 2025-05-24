import numpy as np
from artiq.experiment import *
from artiq.coredevice import ad9910

from LAX_exp.experiment import *
from LAX_exp.system.subsequences import (
    InitializeQubit, SidebandCoolContinuousRAM, QubitPulseShape,
    Readout, RescueIon, NoOperation
)


class CalibrationBichromatic(LAXExperiment, Experiment):
    """
    Calibration: Bichromatic Amplitude

    Calibrate relative amplitudes for bi/multi-chromatic operations via the 729nm singlepass.
    """
    name = 'Calibration Bichromatic'
    kernel_invariants = {
        # subsequences
        'initialize_subsequence', 'pulseshape_subsequence', 'readout_subsequence', 'rescue_subsequence',

        # hardware values
        'singlepass0', 'freq_singlepass0_default_ftw', 'ampl_singlepass0_default_asf',
        'att_singlepass0_default_mu',
        'att_singlepass0_mu', 'ampl_qubit_asf', 'att_qubit_mu',

        # config
        'profile_729_target', 'profile_729_SBC', 'config_experiment_list'
    }

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions",    NumberValue(default=88, precision=0, step=1, min=1, max=100000))
        self.setattr_argument("cooling_type",   EnumerationValue(["Doppler", "SBC - Continuous"],
                                                                 default="SBC - Continuous"))

        # subsequences - with arguments
        self.profile_729_SBC =      5
        self.profile_729_target =   6
        self.sidebandcool_subsequence =  SidebandCoolContinuousRAM(
            self, profile_729=self.profile_729_SBC, profile_854=3,
            ram_addr_start_729=0, ram_addr_start_854=0, num_samples=250
        )

        # beam - qubit (729nm doublepass)
        self.setattr_argument("freq_qubit_mhz", NumberValue(default=101.0965, precision=6, step=1, min=50., max=400., unit="MHz", scale=1.),
                                                group="beam.qubit")
        self.setattr_argument("ampl_qubit_pct", NumberValue(default=50, precision=3, step=5, min=1, max=50, unit="%", scale=1.), group="beam.qubit")
        self.setattr_argument("att_qubit_db",   NumberValue(default=8., precision=1, step=0.5, min=8, max=31.5, unit="dB", scale=1.), group="beam.qubit")
        self.setattr_argument("enable_pulseshaping", BooleanValue(default=False), group="beam.qubit")

        # singlepass beam parameters
        self.setattr_argument("freq_singlepass_center_mhz", NumberValue(default=120.339, precision=6, step=1, min=50., max=400., unit="MHz", scale=1.),
                                                                group="beam.singlepass")
        self.setattr_argument("att_singlepass_db",     NumberValue(default=7., precision=1, step=0.5, min=2., max=31.5, unit="dB", scale=1.), group="beam.singlepass")
        self.setattr_argument("freq_singlepass_sweep_khz_list",  Scannable(
                                                            default=[
                                                                CenterScan(0., 2000., 10, randomize=True),
                                                                ExplicitScan([0.]),
                                                                ExplicitScan([-702.89, 702.89]),
                                                                RangeScan(-702.89, 702.29, 200, randomize=True),
                                                            ],
                                                            global_min=-20000., global_max=200000., global_step=10,
                                                            unit="kHz", scale=1, precision=3
                                                        ), group="beam.singlepass")
        self.setattr_argument("ampl_singlepass_pct_list",   Scannable(
                                                            default=[
                                                                RangeScan(25, 50, 10, randomize=True),
                                                                ExplicitScan([50.]),
                                                                CenterScan(30, 20., 4, randomize=True),
                                                            ],
                                                            global_min=1, global_max=70, global_step=5,
                                                            unit="%", scale=1, precision=3
                                                        ), group="beam.singlepass")

        # scan parameters - time
        self.setattr_argument("equalize_delays",        BooleanValue(default=False), group="scan.time")
        self.setattr_argument("time_rabi_us_list",      Scannable(
                                                            default=[
                                                                RangeScan(1, 25, 10, randomize=True),
                                                                ExplicitScan([6.05]),
                                                                CenterScan(3.05, 5., 0.1, randomize=True),
                                                            ],
                                                            global_min=1, global_max=100000, global_step=1,
                                                            unit="us", scale=1, precision=5
                                                        ), group="scan.time")

        # relevant devices
        self.setattr_device('qubit')

        # subsequences - other
        self.pulseshape_subsequence =   QubitPulseShape(
            self, ram_profile=self.profile_729_target, ram_addr_start=300, num_samples=250,
            ampl_max_pct=self.ampl_qubit_pct, pulse_shape="blackman"
        )
        self.initialize_subsequence =   InitializeQubit(self)
        self.doppler_subsequence =      NoOperation(self)
        self.readout_subsequence =      Readout(self)
        self.rescue_subsequence =       RescueIon(self)

    def prepare_experiment(self):
        """
        Prepare & precompute experimental values.
        """
        # choose target cooling subsequence
        if self.cooling_type == "Doppler":
            self.cooling_subsequence = self.doppler_subsequence
        elif self.cooling_type == "SBC - Continuous":
            self.cooling_subsequence = self.sidebandcool_subsequence

        '''
        CONVERT VALUES TO MACHINE UNITS
        '''
        # get default values for singlepass (so we can set them back afterwards)
        self.singlepass0 = self.get_device("urukul0_ch1")
        self.freq_singlepass0_default_ftw = self.singlepass0.frequency_to_ftw(120.339 * MHz)
        self.ampl_singlepass0_default_asf = self.singlepass0.amplitude_to_asf(0.5)
        self.att_singlepass0_default_mu =   att_to_mu(7. * dB)

        # beam parameters
        self.ampl_qubit_asf =       self.qubit.amplitude_to_asf(self.ampl_qubit_pct / 100.)
        self.att_singlepass0_mu =   att_to_mu(self.att_singlepass_db * dB)
        self.att_qubit_mu =         att_to_mu(self.att_qubit_db * dB)

        # create ampl sweep
        ampl_singlepass_asf_list =  [self.singlepass0.amplitude_to_asf(ampl_pct / 100.)
                                     for ampl_pct in self.ampl_singlepass_pct_list]

        # create frequency config to equalize output at ion
        freq_singlepass_config_ftw_list = np.array([
            [
                # note: singlepass DDS controls SLS single-pass => 1.0x
                self.qubit.frequency_to_ftw(self.freq_singlepass_center_mhz * MHz + freq_khz * kHz),
                # note: qubit DDS controls double-pass => 0.5x
                self.qubit.frequency_to_ftw(self.freq_qubit_mhz * MHz - 0.5 * (freq_khz * kHz))
            ]
            for freq_khz in self.freq_singlepass_sweep_khz_list
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
        self.config_experiment_list = create_experiment_config(
            freq_singlepass_config_ftw_list,
            ampl_singlepass_asf_list,
            time_rabiflop_mu_list,
            shuffle_config=True
        )

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_experiment_list),
                5)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self) -> TNone:
        # ensure singlepass is set up correctly on ALL profiles
        self.singlepass0.set_att_mu(self.att_singlepass0_mu)
        for i in range(8):
            self.singlepass0.set_mu(self.freq_singlepass0_default_ftw,
                                      asf=self.ampl_singlepass0_default_asf,
                                      profile=i)
            delay_mu(8000)
        self.singlepass0.sw.on()

        # record subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.cooling_subsequence.record_dma()
        self.readout_subsequence.record_dma()
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        for trial_num in range(self.repetitions):
            self.core.break_realtime()

            # sweep exp config
            for config_vals in self.config_experiment_list:

                '''PREPARE & CONFIGURE'''
                # extract values from config list
                freq_singlepass_ftw =   np.int32(config_vals[0])
                freq_qubit_ftw =        np.int32(config_vals[1])
                ampl_singlepass_asf =   np.int32(config_vals[2])
                time_equalize_mu =      config_vals[3]
                time_pulse_mu =         config_vals[4]
                self.core.break_realtime()

                # configure pulse times
                if self.enable_pulseshaping:
                    self.pulseshape_subsequence.configure(time_pulse_mu)
                    delay_mu(50000)

                # set 729nm frequencies
                self.singlepass0.set_mu(freq_singlepass_ftw, asf=ampl_singlepass_asf, profile=self.profile_729_target,
                                        phase_mode=ad9910.PHASE_MODE_CONTINUOUS)
                if self.enable_pulseshaping:
                    self.qubit.set_ftw(freq_qubit_ftw)
                else:
                    self.qubit.set_mu(freq_qubit_ftw, asf=self.ampl_qubit_asf, profile=self.profile_729_target,
                                      phase_mode=ad9910.PHASE_MODE_CONTINUOUS)
                delay_mu(10000)

                '''INITIALIZE'''
                # initialize ion in S-1/2 state and cool
                self.initialize_subsequence.run_dma()
                self.cooling_subsequence.run_dma()

                '''MAIN PULSE'''
                # add delay to ensure each shot takes same time
                delay_mu(time_equalize_mu)

                # rabi flop & read out
                if self.enable_pulseshaping:
                    self.pulseshape_subsequence.run()
                else:
                    self.qubit.cpld.set_profile(self.profile_729_target)
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
                                    freq_singlepass_ftw,
                                    time_pulse_mu,
                                    ampl_singlepass_asf)
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
        # set singlepass to default value (b/c AOM thermal drift) on ALL profiles
        for i in range(8):
            self.singlepass0.set_mu(self.freq_singlepass0_default_ftw,
                                      asf=self.ampl_singlepass0_default_asf,
                                      profile=i)
            delay_mu(5000)

        self.singlepass0.sw.on()
        self.singlepass0.set_att_mu(self.att_singlepass0_default_mu)
        delay_mu(25000)

