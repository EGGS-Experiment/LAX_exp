import numpy as np
from artiq.experiment import *

from LAX_exp.analysis import *
from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import (
    InitializeQubit, SidebandCoolContinuousRAM, QubitPulseShape, Readout, RescueIon
)

from itertools import product
from artiq.coredevice import ad9910


class CatStateCharacterize(LAXExperiment, Experiment):
    """
    Experiment: Cat State Characterize

    todo: document
    """
    name = 'Cat State Characterize'
    kernel_invariants = {
        # subsequences
        'initialize_subsequence', 'sidebandcool_subsequence', 'readout_subsequence', 'rescue_subsequence',

        # hardware values - default
        'singlepass0', 'singlepass1',
        'freq_singlepass_default_ftw_list', 'ampl_singlepass_default_asf_list', 'att_singlepass_default_mu_list',
        'ampl_doublepass_default_asf', 'att_doublepass_default_mu',
        'freq_sigmax_ftw', 'ampl_sigmax_asf', 'att_sigmax_mu', 'time_sigmax_mu',

        # hardware values - cat
        'ampls_cat_asf', 'atts_cat_mu', 'time_pulse1_cat_mu', 'phases_pulse1_cat_pow', 'phase_pulse3_sigmax_pow',
        'phases_pulse4_cat_pow',
        'phases_pulse4_cat_update_dir', 'ampl_pulse5_readout_asf', 'att_pulse5_readout_mu',

        # configs
        'profile_729_SBC', 'profile_729_target', 'config_experiment_list'
    }

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions", NumberValue(default=50, precision=0, step=1, min=1, max=100000))

        # allocate relevant beam profiles
        self.profile_729_SBC =      1
        self.profile_729_target =   6

        # get subsequences
        self.sidebandcool_subsequence = SidebandCoolContinuousRAM(
            self, profile_729=self.profile_729_SBC, profile_854=3,
            ram_addr_start_729=0, ram_addr_start_854=0, num_samples=500
        )
        self.initialize_subsequence =   InitializeQubit(self)
        self.readout_subsequence =      Readout(self)
        self.rescue_subsequence =       RescueIon(self)

        '''DEFAULT CONFIG ARGUMENTS'''
        # defaults - beam values
        self.max_ampl_singlepass_pct, self.min_att_singlepass_db = (58., 3.)
        self.setattr_argument("freq_singlepass_default_mhz_list",   PYONValue([80., 80.]), group='defaults.beams', tooltip="[rsb_mhz, bsb_mhz]")
        self.setattr_argument("ampl_singlepass_default_pct_list",   PYONValue([58., 58.]), group='defaults.beams', tooltip="[rsb_pct, bsb_pct]")
        self.setattr_argument("att_singlepass_default_db_list",     PYONValue([3., 3.]), group='defaults.beams', tooltip="[rsb_db, bsb_db]")

        self.setattr_argument("ampl_doublepass_default_pct",    NumberValue(default=50., precision=3, step=5, min=0.01, max=50), group="defaults.beams")
        self.setattr_argument("att_doublepass_default_db",      NumberValue(default=8., precision=1, step=0.5, min=8., max=31.5), group="defaults.beams")

        # defaults - sigma_x
        self.setattr_argument("freq_sigmax_mhz",    NumberValue(default=101.2401, precision=6, step=1, min=50., max=400.), group="defaults.sigmax")
        self.setattr_argument("ampl_sigmax_pct",    NumberValue(default=50., precision=3, step=5, min=0.01, max=50), group="defaults.sigmax")
        self.setattr_argument("att_sigmax_db",      NumberValue(default=8., precision=1, step=0.5, min=8., max=31.5), group="defaults.sigmax")
        self.setattr_argument("time_sigmax_us",     NumberValue(default=3.21, precision=2, step=5, min=0.1, max=10000), group="defaults.sigmax")

        # defaults - cat
        self.setattr_argument("freq_cat_center_mhz_list",   Scannable(
                                                                default=[
                                                                    ExplicitScan([101.2401]),
                                                                    CenterScan(101.2401, 0.01, 0.0001, randomize=True),
                                                                    RangeScan(101.2301, 101.2501, 50, randomize=True),
                                                                ],
                                                                global_min=60., global_max=400, global_step=1,
                                                                unit="MHz", scale=1, precision=6
                                                            ), group='defaults.cat')
        self.setattr_argument("freq_cat_secular_khz_list",  Scannable(
                                                                default=[
                                                                    ExplicitScan([702.01]),
                                                                    CenterScan(702.01, 4, 0.1, randomize=True),
                                                                    RangeScan(699.0, 704.2, 50, randomize=True),
                                                                ],
                                                                global_min=0, global_max=10000, global_step=1,
                                                                unit="kHz", scale=1, precision=3
                                                            ), group='defaults.cat')
        self.setattr_argument("ampls_cat_pct",  PYONValue([58., 58.]), group='defaults.cat', tooltip="[rsb_pct, bsb_pct]")
        self.setattr_argument("atts_cat_db",    PYONValue([5., 5.]), group='defaults.cat', tooltip="[rsb_db, bsb_db]")


        '''PULSE ARGUMENTS - CAT 1'''
        # pulse 0 - sigma_x0
        self.setattr_argument("enable_pulse0_sigmax", BooleanValue(default=False), group='pulse0.sigmax')

        # pulse 1 - cat 0
        self.setattr_argument("enable_pulse1_cat", BooleanValue(default=False), group='pulse1.cat')
        self.setattr_argument("time_pulse1_cat_us", NumberValue(default=100, precision=2, step=5, min=0.1, max=10000), group="pulse1.cat")
        self.setattr_argument("phases_pulse1_cat_turns",  PYONValue([0., 0.]), group='pulse1.cat', tooltip="[rsb_turns, bsb_turns]")

        # pulse 2 - quench
        self.setattr_argument("enable_pulse2_herald",   BooleanValue(default=False), group='pulse2.quench')
        self.setattr_argument("enable_pulse2_quench",   BooleanValue(default=True), group='pulse2.quench')

        '''PULSE ARGUMENTS - CAT 2'''
        # pulse 3 - sigma_x1
        self.setattr_argument("enable_pulse3_sigmax",       BooleanValue(default=True), group='pulse3.sigmax')
        self.setattr_argument("phase_pulse3_sigmax_turns",  NumberValue(default=0., precision=3, step=0.1, min=-1.0, max=1.0), group='pulse3.sigmax')

        # pulse 4 - cat 1
        self.setattr_argument("enable_pulse4_cat",  BooleanValue(default=True), group='pulse4.cat')
        self.setattr_argument("enable_pulse4_herald",   BooleanValue(default=True), group='pulse4.cat')
        self.setattr_argument("enable_pulse4_quench",   BooleanValue(default=True), group='pulse4.cat')
        self.setattr_argument("time_pulse4_cat_us_list",    Scannable(
                                                        default=[
                                                            ExplicitScan([100]),
                                                            RangeScan(0, 500, 50, randomize=True),
                                                        ],
                                                        global_min=1, global_max=100000, global_step=1,
                                                        unit="us", scale=1, precision=5
                                                    ), group="pulse4.cat")
        self.setattr_argument("phases_pulse4_cat_turns",        PYONValue([0., 0.]), group='pulse4.cat', tooltip="[rsb_turns, bsb_turns]")
        self.setattr_argument("target_pulse4_cat_phase",        EnumerationValue(['RSB', 'BSB', 'RSB-BSB', 'RSB+BSB'], default='RSB-BSB'), group="pulse4.cat")
        self.setattr_argument("phase_pulse4_cat_turns_list",    Scannable(
                                                                    default=[
                                                                        ExplicitScan([0.]),
                                                                        RangeScan(0, 1.0, 11, randomize=True),
                                                                    ],
                                                                    global_min=-1.0, global_max=1.0, global_step=0.1,
                                                                    unit="turns", scale=1, precision=3
                                                                ), group="pulse4.cat")

        '''PULSE ARGUMENTS - READOUT'''
        # pulse 5 - readout
        self.setattr_argument("enable_pulse5_readout",      BooleanValue(default=True), group="pulse5.readout")
        self.setattr_argument("ampl_pulse5_readout_pct",    NumberValue(default=50., precision=3, step=5, min=0.01, max=50), group="pulse5.readout")
        self.setattr_argument("att_pulse5_readout_db",      NumberValue(default=8., precision=1, step=0.5, min=8., max=31.5), group="pulse5.readout")
        self.setattr_argument("freq_pulse5_readout_mhz_list",   Scannable(
                                                                default=[
                                                                    ExplicitScan([101.9851]),
                                                                    CenterScan(101.9851, 0.01, 0.0002, randomize=True),
                                                                    RangeScan(101.9801, 101.9901, 50, randomize=True),
                                                                ],
                                                                global_min=60., global_max=400, global_step=1,
                                                                unit="MHz", scale=1, precision=6
                                                            ), group="pulse5.readout")
        self.setattr_argument("time_pulse5_readout_us_list",    Scannable(
                                                                default=[
                                                                    ExplicitScan([122.9]),
                                                                    RangeScan(0, 1000, 100, randomize=True),
                                                                ],
                                                                global_min=1, global_max=100000, global_step=1,
                                                                unit="us", scale=1, precision=5
                                                            ), group="pulse5.readout")

        # relevant devices
        self.setattr_device('qubit')
        self.setattr_device('pump')
        self.setattr_device('repump_qubit')

    def prepare_experiment(self):
        """
        Prepare & precompute experimental values.
        """
        '''
        CHECK INPUT FOR ERRORS
        '''
        # ensure single pass values are safe and valid
        if any((ampl_pct > self.max_ampl_singlepass_pct or ampl_pct < 0.
                for ampl_pct in self.ampl_singlepass_default_pct_list)):
            raise ValueError("Singlepass amplitude outside valid range - [0., {:f}].".format(self.max_ampl_singlepass_pct))
        
        if any((att_db > 31.5 or att_db < self.min_att_singlepass_db
                for att_db in self.att_singlepass_default_db_list)):
            raise ValueError("Singlepass attenuation outside valid range - [{:.1f}, 31.5].".format(self.min_att_singlepass_db))
        
        # ensure we only herald once
        if self.enable_pulse2_herald and self.enable_pulse4_herald:
            raise ValueError("Cannot herald twice. Must choose either pulse2_herald OR pulse4_herald.")

        '''
        CONVERT VALUES TO MACHINE UNITS - DEFAULTS
        '''
        # defaults - singlepass AOM
        self.singlepass0 = self.get_device("urukul0_ch1")
        self.singlepass1 = self.get_device("urukul0_ch2")

        self.freq_singlepass_default_ftw_list = [self.singlepass0.frequency_to_ftw(freq_mhz * MHz)
                                                 for freq_mhz in self.freq_singlepass_default_mhz_list]
        self.ampl_singlepass_default_asf_list = [self.singlepass0.amplitude_to_asf(ampl_asf / 100.)
                                                 for ampl_asf in self.ampl_singlepass_default_pct_list]
        self.att_singlepass_default_mu_list =   [att_to_mu(att_db * dB)
                                                 for att_db in self.att_singlepass_default_db_list]

        # defaults - doublepass AOM
        self.ampl_doublepass_default_asf =     self.qubit.amplitude_to_asf(self.ampl_doublepass_default_pct / 100.)
        self.att_doublepass_default_mu =       att_to_mu(self.att_doublepass_default_db * dB)

        # defaults - sigma_x pulses
        self.freq_sigmax_ftw =  self.qubit.frequency_to_ftw(self.freq_sigmax_mhz * MHz)
        self.ampl_sigmax_asf =  self.qubit.amplitude_to_asf(self.ampl_sigmax_pct / 100.)
        self.att_sigmax_mu =    att_to_mu(self.att_sigmax_db * dB)
        self.time_sigmax_mu =   self.core.seconds_to_mu(self.time_sigmax_us * us)

        # defaults - cat
        freq_cat_center_ftw_list =  np.array([self.qubit.frequency_to_ftw(freq_mhz * MHz)
                                              for freq_mhz in self.freq_cat_center_mhz_list])
        freq_cat_secular_ftw_list = np.array([self.singlepass0.frequency_to_ftw(freq_khz * kHz)
                                              for freq_khz in self.freq_cat_secular_khz_list])
        self.ampls_cat_asf =    np.array([self.singlepass0.amplitude_to_asf(ampl_pct / 100.)
                                          for ampl_pct in self.ampls_cat_pct], dtype=np.int32)
        self.atts_cat_mu =      np.array([att_to_mu(att_db * dB)
                                          for att_db in self.atts_cat_db], dtype=np.int32)

        '''
        CONVERT VALUES TO MACHINE UNITS - PULSES
        '''
        # pulse 1 - cat 0
        self.time_pulse1_cat_mu =       self.core.seconds_to_mu(self.time_pulse1_cat_us * us)
        self.phases_pulse1_cat_pow =    [self.singlepass0.turns_to_pow(phas_pow)
                                         for phas_pow in self.phases_pulse1_cat_turns]

        # pulse 3 - sigma_x1
        self.phase_pulse3_sigmax_pow =  self.qubit.turns_to_pow(self.phase_pulse3_sigmax_turns)

        # pulse 4 - cat 1
        self.phases_pulse4_cat_pow =    [self.singlepass0.turns_to_pow(phas_pow)
                                         for phas_pow in self.phases_pulse4_cat_turns]
        if self.enable_pulse4_cat:
            time_pulse4_cat_mu_list =   np.array([self.core.seconds_to_mu(time_us * us)
                                                 for time_us in self.time_pulse4_cat_us_list], dtype=np.int64)
            phase_pulse4_cat_pow_list = np.array([self.singlepass0.turns_to_pow(phas_pow)
                                                  for phas_pow in self.phase_pulse4_cat_turns_list], dtype=np.int32)
        else:
            time_pulse4_cat_mu_list =   np.array([0], dtype=np.int64)
            phase_pulse4_cat_pow_list = np.array([0], dtype=np.int32)

        if self.target_pulse4_cat_phase == 'RSB':
            self.phases_pulse4_cat_update_dir = np.array([1, 0], dtype=np.int32)
        elif self.target_pulse4_cat_phase == 'BSB':
            self.phases_pulse4_cat_update_dir = np.array([0, 1], dtype=np.int32)
        elif self.target_pulse4_cat_phase == 'RSB-BSB':
            self.phases_pulse4_cat_update_dir = np.array([1, -1], dtype=np.int32)
        elif self.target_pulse4_cat_phase == 'RSB+BSB':
            self.phases_pulse4_cat_update_dir = np.array([1, 1], dtype=np.int32)

        # pulse 5 - readout
        self.ampl_pulse5_readout_asf =  self.qubit.amplitude_to_asf(self.ampl_pulse5_readout_pct / 100.)
        self.att_pulse5_readout_mu =    att_to_mu(self.att_pulse5_readout_db * dB)

        if self.enable_pulse5_readout:
            freq_pulse5_readout_ftw_list =  np.array([self.qubit.frequency_to_ftw(freq_mhz * MHz)
                                                      for freq_mhz in self.freq_pulse5_readout_mhz_list], dtype=np.int32)
            time_pulse5_readout_mu_list =   np.array([self.core.seconds_to_mu(time_us * us)
                                                      for time_us in self.time_pulse5_readout_us_list], dtype=np.int64)
        else:
            freq_pulse5_readout_ftw_list =  np.array([0], dtype=np.int32)
            time_pulse5_readout_mu_list =   np.array([0], dtype=np.int64)

        '''
        CREATE EXPERIMENT CONFIG
        '''
        # create an array of values for the experiment to sweep
        self.config_experiment_list = np.array([
            vals
            for vals in product(freq_cat_center_ftw_list, freq_cat_secular_ftw_list,
                                time_pulse4_cat_mu_list, phase_pulse4_cat_pow_list,
                                freq_pulse5_readout_ftw_list, time_pulse5_readout_mu_list)
        ], dtype=np.int64)
        np.random.shuffle(self.config_experiment_list)

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_experiment_list),
                8)


    '''
    MAIN SEQUENCE
    '''
    @kernel(flags={"fast-math"})
    def initialize_experiment(self) -> TNone:
        self.core.break_realtime()

        # set up beam parameters
        self.qubit.set_att_mu(self.att_doublepass_default_mu)
        self.core.break_realtime()

        # ensure phase_autoclear disabled on all beams to prevent phase accumulator reset
        # enable RAM mode and clear DDS phase accumulator
        self.qubit.set_cfr1()
        self.singlepass0.set_cfr1()
        self.singlepass1.set_cfr1()
        self.qubit.cpld.io_update.pulse_mu(8)
        self.core.break_realtime()

        # set up singlepass AOMs to default values (b/c AOM thermal drift) on ALL profiles
        for i in range(8):
            self.singlepass0.set_mu(self.freq_singlepass_default_ftw_list[0],
                                      asf=self.ampl_singlepass_default_asf_list[0],
                                      profile=i)
            self.singlepass1.set_mu(self.freq_singlepass_default_ftw_list[1],
                                      asf=self.ampl_singlepass_default_asf_list[1],
                                      profile=i)
            self.singlepass0.cpld.io_update.pulse_mu(8)
            delay_mu(8000)
        self.core.break_realtime()

        self.singlepass0.set_att_mu(self.att_singlepass_default_mu_list[0])
        self.singlepass1.set_att_mu(self.att_singlepass_default_mu_list[1])
        self.core.break_realtime()

        self.singlepass0.sw.on()
        self.singlepass1.sw.off()
        self.core.break_realtime()

        # record general subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.sidebandcool_subsequence.record_dma()
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
                freq_cat_center_ftw =       np.int32(config_vals[0])
                freq_cat_secular_ftw =      np.int32(config_vals[1])
                time_pulse4_cat_mu =        config_vals[2]
                phase_pulse4_cat_pow =      np.int32(config_vals[3])
                freq_pulse5_readout_ftw =   np.int32(config_vals[4])
                time_pulse5_readout_mu =    config_vals[5]
                self.core.break_realtime()

                # create phase list
                cat4_phases = [
                    self.phases_pulse4_cat_pow[0] + self.phases_pulse4_cat_update_dir[0] * phase_pulse4_cat_pow,
                    self.phases_pulse4_cat_pow[1] + self.phases_pulse4_cat_update_dir[1] * phase_pulse4_cat_pow,
                ]

                '''INITIALIZE'''
                # initialize ion in S-1/2 state
                self.initialize_subsequence.run_dma()
                # sideband cool
                self.sidebandcool_subsequence.run_dma()

                # set target profile to ensure we run correctly
                self.qubit.set_profile(self.profile_729_target)
                self.qubit.cpld.io_update.pulse_mu(8)

                # synchronize start time to coarse RTIO clock
                time_start_mu = now_mu() & ~0x7

                '''CAT #1'''
                # pulse 0: sigma_x #1
                if self.enable_pulse0_sigmax:
                    self.pulse_sigmax(time_start_mu, 0)

                # pulse 1: cat 1
                if self.enable_pulse1_cat:
                    self.pulse_bichromatic(time_start_mu, self.time_pulse1_cat_mu,
                                           self.phases_pulse1_cat_pow,
                                           freq_cat_center_ftw, freq_cat_secular_ftw)

                # pulse 2a: herald via 397nm
                if self.enable_pulse2_herald:
                    self.readout_subsequence.run_dma()
                    self.pump.off()

                # pulse 2b: repump via 854
                if self.enable_pulse2_quench:
                    self.pump.readout()
                    self.repump_qubit.on()
                    delay_mu(self.initialize_subsequence.time_repump_qubit_mu)
                    self.repump_qubit.off()

                '''CAT #2'''
                # pulse 3: sigma_x #2
                if self.enable_pulse3_sigmax:
                    self.pulse_sigmax(time_start_mu, self.phase_pulse3_sigmax_pow)

                # pulse 4a: cat #2
                if self.enable_pulse4_cat:
                    self.pulse_bichromatic(time_start_mu, time_pulse4_cat_mu,
                                           cat4_phases,
                                           freq_cat_center_ftw, freq_cat_secular_ftw)

                # pulse 4a: herald via 397nm
                if self.enable_pulse4_herald:
                    self.readout_subsequence.run_dma()
                    self.pump.off()

                # pulse 4b: repump via 854
                if self.enable_pulse4_quench:
                    self.pump.readout()
                    self.repump_qubit.on()
                    delay_mu(self.initialize_subsequence.time_repump_qubit_mu)
                    self.repump_qubit.off()

                '''READOUT & STORE RESULTS'''
                # pulse 5: rabiflop/readout
                if self.enable_pulse5_readout:
                    self.pulse_readout(time_pulse5_readout_mu, freq_pulse5_readout_ftw)

                # read out
                self.readout_subsequence.run_dma()

                # get heralded measurement and actual results
                if self.enable_pulse2_herald or self.enable_pulse4_herald:
                    counts_her = self.readout_subsequence.fetch_count()
                else:
                    counts_her = 0
                counts_res = self.readout_subsequence.fetch_count()
                self.update_results(freq_cat_center_ftw,
                                    counts_res,
                                    counts_her,
                                    freq_cat_secular_ftw,
                                    time_pulse4_cat_mu,
                                    phase_pulse4_cat_pow,
                                    freq_pulse5_readout_ftw,
                                    time_pulse5_readout_mu)
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

        # set up singlepass AOMs to default values (b/c AOM thermal drift) on ALL profiles
        for i in range(8):
            self.singlepass0.set_mu(self.freq_singlepass_default_ftw_list[0],
                                    asf=self.ampl_singlepass_default_asf_list[0],
                                    profile=i)
            self.singlepass1.set_mu(self.freq_singlepass_default_ftw_list[1],
                                    asf=self.ampl_singlepass_default_asf_list[1],
                                    profile=i)
            self.singlepass0.cpld.io_update.pulse_mu(8)
            delay_mu(8000)
        self.core.break_realtime()

        self.singlepass0.set_att_mu(self.att_singlepass_default_mu_list[0])
        self.singlepass1.set_att_mu(self.att_singlepass_default_mu_list[1])
        self.core.break_realtime()

        self.singlepass0.sw.on()
        self.singlepass1.sw.off()
        self.core.break_realtime()


    '''
    HELPER FUNCTIONS
    '''
    @kernel(flags={"fast-math"})
    def pulse_sigmax(self, time_start_mu: TInt64, phas_pow: TInt32) -> TNone:
        """
        Run a phase-coherent sigma_x pulse on the qubit.
        Arguments:
            time_start_mu: fiducial timestamp for initial start reference (in machine units).
            phas_pow: relative phase offset for the beam.
        """
        # set up relevant beam waveforms
        self.qubit.set_mu(
            self.freq_sigmax_ftw, asf=self.ampl_sigmax_asf, pow_=phas_pow,
            profile=self.profile_729_target, phase_mode=ad9910.PHASE_MODE_TRACKING, ref_time_mu=time_start_mu
        )
        self.singlepass0.set_mu(
            self.freq_singlepass_default_ftw_list[0], asf=self.ampl_singlepass_default_asf_list[0], pow_=0,
            profile=self.profile_729_target, phase_mode=ad9910.PHASE_MODE_TRACKING, ref_time_mu=time_start_mu
        )
        self.singlepass1.set_mu(
            self.freq_singlepass_default_ftw_list[1], asf=self.ampl_singlepass_default_asf_list[1], pow_=0,
            profile=self.profile_729_target, phase_mode=ad9910.PHASE_MODE_TRACKING, ref_time_mu=time_start_mu
        )
        self.qubit.cpld.io_update.pulse_mu(8)

        # set all attenuators together
        a = self.qubit.cpld.att_reg & ~(
                (0xFF << (0 * 8)) |
                (0xFF << (1 * 8)) |
                (0xFF << (2 * 8))
        )
        a |= (
                (self.att_sigmax_mu << (0 * 8)) |
                (self.att_singlepass_default_mu_list[0] << (1 * 8)) |
                (self.att_singlepass_default_mu_list[1] << (1 * 8))
        )
        self.qubit.cpld.set_all_att_mu(a)

        # run sigmax pulse
        self.singlepass0.sw.on()
        self.singlepass1.sw.off()
        self.qubit.on()
        delay_mu(self.time_sigmax_mu)
        self.qubit.off()

    @kernel(flags={"fast-math"})
    def pulse_bichromatic(self, time_start_mu: TInt64, time_pulse_mu: TInt64, phas_pow_list: TList(TInt32),
                          freq_carrier_ftw: TInt32, freq_secular_ftw: TInt32) -> TNone:
        """
        Run a phase-coherent bichromatic pulse on the qubit.
        Arguments:
            time_start_mu: fiducial timestamp for initial start reference (in machine units).
            time_pulse_mu: length of pulse (in machine units).
            phas_pow_list: relative phase offset for the beams (RSB, BSB) (in pow).
            freq_carrier_ftw: carrier frequency (set by the double pass) in FTW.
            freq_secular_ftw: bichromatic separation frequency (from central frequency) in FTW.
        """
        # set up relevant beam waveforms
        self.qubit.set_mu(
            freq_carrier_ftw, asf=self.ampl_sigmax_asf, pow_=0,
            profile=self.profile_729_target, phase_mode=ad9910.PHASE_MODE_TRACKING, ref_time_mu=time_start_mu
        )
        self.singlepass0.set_mu(
            self.freq_singlepass_default_ftw_list[0] - freq_secular_ftw, asf=self.ampls_cat_asf[0],
            pow_=phas_pow_list[0], profile=self.profile_729_target,
            phase_mode=ad9910.PHASE_MODE_TRACKING, ref_time_mu=time_start_mu
        )
        self.singlepass1.set_mu(
            self.freq_singlepass_default_ftw_list[1] + freq_secular_ftw, asf=self.ampls_cat_asf[1],
            pow_=phas_pow_list[1], profile=self.profile_729_target,
            phase_mode=ad9910.PHASE_MODE_TRACKING, ref_time_mu=time_start_mu
        )

        # set all attenuators together
        a = self.qubit.cpld.att_reg & ~(
                (0xFF << (0 * 8)) |
                (0xFF << (1 * 8)) |
                (0xFF << (2 * 8))
        )
        a |= (
                (self.att_doublepass_default_mu << (0 * 8)) |
                (self.atts_cat_mu[0] << (1 * 8)) |
                (self.atts_cat_mu[1] << (2 * 8))
        )
        self.qubit.cpld.set_all_att_mu(a)

        # run bichromatic pulse
        self.singlepass0.sw.on()
        self.singlepass1.sw.on()
        self.qubit.on()
        delay_mu(time_pulse_mu)
        self.qubit.off()
        self.singlepass1.sw.off()

    @kernel(flags={"fast-math"})
    def pulse_readout(self, time_pulse_mu: TInt64, freq_readout_ftw: TInt32) -> TNone:
        """
        Run a readout pulse.
        Arguments:
            time_pulse_mu: length of pulse (in machine units).
            freq_readout_ftw: readout frequency (set by the double pass) in FTW.
        """
        # set up relevant beam waveforms
        self.qubit.set_mu(freq_readout_ftw, asf=self.ampl_pulse5_readout_asf, pow_=0,
                          profile=self.profile_729_target)
        self.singlepass0.set_mu(self.freq_singlepass_default_ftw_list[0],
                                asf=self.ampl_singlepass_default_asf_list[0],
                                pow_=0, profile=self.profile_729_target)
        self.singlepass1.set_mu(self.freq_singlepass_default_ftw_list[1],
                                asf=self.ampl_singlepass_default_asf_list[1],
                                pow_=0, profile=self.profile_729_target)

        # set all attenuators together
        a = self.qubit.cpld.att_reg & ~(
                (0xFF << (0 * 8)) |
                (0xFF << (1 * 8)) |
                (0xFF << (2 * 8))
        )
        a |= (
                (self.att_pulse5_readout_mu << (0 * 8)) |
                (self.att_singlepass_default_mu_list[0] << (1 * 8)) |
                (self.att_singlepass_default_mu_list[1] << (1 * 8))
        )
        self.qubit.cpld.set_all_att_mu(a)

        # run readout pulse
        self.singlepass0.sw.on()
        self.singlepass1.sw.off()
        self.qubit.on()
        delay_mu(time_pulse_mu)
        self.qubit.off()


    '''
    ANALYSIS
    '''
    def analyze_experiment(self):
        pass

