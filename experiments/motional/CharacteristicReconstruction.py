import numpy as np
from artiq.experiment import *

from LAX_exp.analysis import *
from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import (
    InitializeQubit, SidebandCoolContinuousRAM, Readout, ReadoutAdaptive, RescueIon
)

import math
from itertools import product
from collections.abc import Iterable
from artiq.coredevice import ad9910


class CharacteristicReconstruction(LAXExperiment, Experiment):
    """
    Experiment: Characteristic Reconstruction

    Directly reconstruct the Characteristic function of a given motional state using the Fluhmann technique
    (https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.125.043602).
    """
    name = 'Characteristic Reconstruction'
    kernel_invariants = {
        # subsequences
        'initialize_subsequence', 'sidebandcool_subsequence', 'readout_subsequence', 'readout_adaptive_subsequence',
        'rescue_subsequence',

        # hardware objects
        'singlepass0', 'singlepass1',

        # hardware parameters
        'freq_singlepass_default_ftw_list', 'ampl_singlepass_default_asf_list', 'att_singlepass_default_mu_list',
        'ampl_doublepass_default_asf', 'freq_sigmax_ftw', 'ampl_sigmax_asf', 'time_sigmax_mu',
        'time_force_herald_slack_mu', 'att_reg_bichromatic', 'att_reg_sigmax',

        # cat state parameters
        'ampls_cat_asf', 'time_motion_cat_mu', 'phases_motion_cat_pow', 'phase_char_axis_pow',
        'phases_char_cat_pow', 'phases_char_cat_update_dir',

        # experiment parameters
        'profile_729_SBC', 'profile_729_target', 'config_experiment_list'
    }

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions", NumberValue(default=50, precision=0, step=1, min=1, max=100000))

        # set target profile for stuff
        self.profile_729_SBC =      5
        self.profile_729_target =   6

        # get subsequences
        self.sidebandcool_subsequence = SidebandCoolContinuousRAM(
            self, profile_729=self.profile_729_SBC, profile_854=3,
            ram_addr_start_729=0, ram_addr_start_854=0,
            num_samples=200
        )
        self.initialize_subsequence =   InitializeQubit(self)
        self.readout_subsequence =      Readout(self)
        self.readout_adaptive_subsequence = ReadoutAdaptive(self, time_bin_us=10, error_threshold=1e-2)
        self.rescue_subsequence =       RescueIon(self)


        '''DEFAULT CONFIG ARGUMENTS'''
        # defaults - beam values
        self.max_ampl_singlepass_pct, self.min_att_singlepass_db = (58., 3.)
        self.setattr_argument("freq_singlepass_default_mhz_list",   PYONValue([120.339, 120.339]), group='defaults.beams', tooltip="[rsb_mhz, bsb_mhz]")
        self.setattr_argument("ampl_singlepass_default_pct_list",   PYONValue([50., 0.01]), group='defaults.beams', tooltip="[rsb_pct, bsb_pct]")
        self.setattr_argument("att_singlepass_default_db_list",     PYONValue([7., 7.]), group='defaults.beams', tooltip="[rsb_db, bsb_db]")
        self.setattr_argument("ampl_doublepass_default_pct",    NumberValue(default=50., precision=3, step=5, min=0.01, max=50), group="defaults.beams")
        self.setattr_argument("att_doublepass_default_db",      NumberValue(default=8., precision=1, step=0.5, min=8., max=31.5), group="defaults.beams")

        # defaults - sigma_x
        self.setattr_argument("freq_sigmax_mhz",    NumberValue(default=101.1013, precision=6, step=1, min=50., max=400.), group="defaults.sigmax")
        self.setattr_argument("ampl_sigmax_pct",    NumberValue(default=50., precision=3, step=5, min=0.01, max=50), group="defaults.sigmax")
        self.setattr_argument("att_sigmax_db",      NumberValue(default=8., precision=1, step=0.5, min=8., max=31.5), group="defaults.sigmax")
        self.setattr_argument("time_sigmax_us",     NumberValue(default=1.59, precision=2, step=5, min=0.1, max=10000), group="defaults.sigmax")

        # defaults - bichromatic
        self.setattr_argument("freq_cat_center_mhz_list",   Scannable(
                                                                default=[
                                                                    ExplicitScan([101.1013]),
                                                                    CenterScan(101.1013, 0.01, 0.0001, randomize=True),
                                                                    RangeScan(101.1005, 101.1021, 50, randomize=True),
                                                                ],
                                                                global_min=60., global_max=400, global_step=1,
                                                                unit="MHz", scale=1, precision=6
                                                            ), group='defaults.bichromatic')
        self.setattr_argument("freq_cat_secular_khz_list",  Scannable(
                                                                default=[
                                                                    ExplicitScan([703.101]),
                                                                    CenterScan(703.1, 4, 0.1, randomize=True),
                                                                    RangeScan(701.0, 704.0, 50, randomize=True),
                                                                ],
                                                                global_min=0, global_max=10000, global_step=1,
                                                                unit="kHz", scale=1, precision=3
                                                            ), group='defaults.bichromatic')
        self.setattr_argument("ampls_cat_pct",  PYONValue([50., 50.]), group='defaults.bichromatic', tooltip="[rsb_pct, bsb_pct]")
        self.setattr_argument("atts_cat_db",    PYONValue([13., 13.]), group='defaults.bichromatic', tooltip="[rsb_db, bsb_db]")

        '''ARGUMENTS - MOTIONAL STATE PREPARATION'''
        self.setattr_argument("enable_motion_sigmax",   BooleanValue(default=False), group='motional.config',
                              tooltip='sigma_x selects whether motional state is coherent (True), or cat (False).')
        self.setattr_argument("enable_motion_cat",      BooleanValue(default=False), group='motional.config')
        self.setattr_argument("enable_motion_herald",   BooleanValue(default=True), group='motional.config')
        self.setattr_argument("time_motion_cat_us",     NumberValue(default=100, precision=2, step=5, min=0.1, max=10000), group="motional.config")
        self.setattr_argument("phases_motion_cat_turns",    PYONValue([0., 0.]), group='motional.config', tooltip="[rsb_turns, bsb_turns]")

        '''ARGUMENTS - CHARACTERISTIC FUNCTION MEASUREMENT'''
        # sigma_x: select real/imag part of characteristic function
        self.setattr_argument("characteristic_axis",        EnumerationValue(['Both', 'Real', 'Imaginary'], default='Both'), group='characteristic.axis',
                              tooltip="Selects the real/imag component of the characteristic function by either applying a sigma_x operation (Imag), or not (Real)."
                                      "The 'Both' option enables measurement of both real and imag components within a single experiment.")
        self.setattr_argument("phase_char_axis_turns",  NumberValue(default=0.125, precision=3, step=0.1, min=-1.0, max=1.0), group='characteristic.axis',
                              tooltip="Sets the relative phase of the sigma_x operation used to define the real/imag axis of the characteristic function.")

        # bichromatic: characteristic readout protocol
        self.setattr_argument("phases_char_cat_turns",    PYONValue([0., 0.]), group='characteristic.read', tooltip="[rsb_turns, bsb_turns]")
        self.setattr_argument("target_char_cat_phase",    EnumerationValue(['RSB', 'BSB', 'RSB-BSB', 'RSB+BSB'], default='RSB-BSB'), group='characteristic.read')
        self.setattr_argument("time_char_cat_x_us_list",      Scannable(
                                                                default=[
                                                                    RangeScan(-500, 500, 10, randomize=True),
                                                                    ExplicitScan([100]),
                                                                ],
                                                                global_min=-100000, global_max=100000, global_step=1,
                                                                unit="us", scale=1, precision=5
                                                            ), group='characteristic.read')
        self.setattr_argument("time_char_cat_y_us_list",      Scannable(
                                                                default=[
                                                                    RangeScan(-500, 500, 10, randomize=True),
                                                                    ExplicitScan([100]),
                                                                ],
                                                                global_min=-100000, global_max=100000, global_step=1,
                                                                unit="us", scale=1, precision=5
                                                            ), group='characteristic.read')

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
        self._prepare_argument_checks()

        '''
        CONVERT VALUES TO MACHINE UNITS - DEFAULTS
        '''
        # defaults - singlepass AOMs
        self.singlepass0 = self.get_device("urukul0_ch1")
        self.singlepass1 = self.get_device("urukul0_ch2")
        self.freq_singlepass_default_ftw_list = [self.singlepass0.frequency_to_ftw(freq_mhz * MHz)
                                                 for freq_mhz in self.freq_singlepass_default_mhz_list]
        self.ampl_singlepass_default_asf_list = [self.singlepass0.amplitude_to_asf(ampl_asf / 100.)
                                                 for ampl_asf in self.ampl_singlepass_default_pct_list]
        self.ampl_doublepass_default_asf =     self.qubit.amplitude_to_asf(self.ampl_doublepass_default_pct / 100.)

        # defaults - sigma_x waveform
        self.freq_sigmax_ftw =  self.qubit.frequency_to_ftw(self.freq_sigmax_mhz * MHz)
        self.ampl_sigmax_asf =  self.qubit.amplitude_to_asf(self.ampl_sigmax_pct / 100.)
        self.time_sigmax_mu =   self.core.seconds_to_mu(self.time_sigmax_us * us)

        # defaults - cat
        freq_cat_center_ftw_list =  np.array([self.qubit.frequency_to_ftw(freq_mhz * MHz)
                                              for freq_mhz in self.freq_cat_center_mhz_list])
        freq_cat_secular_ftw_list = np.array([self.singlepass0.frequency_to_ftw(freq_khz * kHz)
                                              for freq_khz in self.freq_cat_secular_khz_list])
        self.ampls_cat_asf =    np.array([self.singlepass0.amplitude_to_asf(ampl_pct / 100.)
                                          for ampl_pct in self.ampls_cat_pct], dtype=np.int32)

        '''
        CONVERT VALUES TO MACHINE UNITS
        '''
        # motional state - cat/bichromatic
        self.time_motion_cat_mu =       self.core.seconds_to_mu(self.time_motion_cat_us * us)
        self.phases_motion_cat_pow =    [self.singlepass0.turns_to_pow(phas_pow)
                                         for phas_pow in self.phases_motion_cat_turns]

        # define characteristic axis (via sigma_x)
        self.phase_char_axis_pow =  self.qubit.turns_to_pow(self.phase_char_axis_turns)
        # configure whether real/imag/both axes of the characteristic function are to be measured
        if self.characteristic_axis == "Real":          characteristic_axis_list = [False]
        elif self.characteristic_axis == "Imaginary":   characteristic_axis_list = [True]
        elif self.characteristic_axis == "Both":        characteristic_axis_list = [True, False]

        # characteristic readout - bichromatic
        self.phases_char_cat_pow =    np.array([self.singlepass0.turns_to_pow(phas_pow)
                                                  for phas_pow in self.phases_char_cat_turns], dtype=np.int32)
        if self.target_char_cat_phase == 'RSB':
            self.phases_char_cat_update_dir = np.array([1, 0], dtype=np.int32)
        elif self.target_char_cat_phase == 'BSB':
            self.phases_char_cat_update_dir = np.array([0, 1], dtype=np.int32)
        elif self.target_char_cat_phase == 'RSB-BSB':
            self.phases_char_cat_update_dir = np.array([1, -1], dtype=np.int32)
        elif self.target_char_cat_phase == 'RSB+BSB':
            self.phases_char_cat_update_dir = np.array([1, 1], dtype=np.int32)

        '''
        CREATE ATTENUATION REGISTERS
        '''
        atts_cat_mu = np.array([att_to_mu(att_db * dB) for att_db in self.atts_cat_db], dtype=np.int32)
        self.att_singlepass_default_mu_list =   [att_to_mu(att_db * dB)
                                                 for att_db in self.att_singlepass_default_db_list]
        self.att_reg_sigmax = 0x00000000 | (
                (att_to_mu(self.att_sigmax_db * dB) << ((self.qubit.beam.chip_select - 4) * 8)) |
                (self.att_singlepass_default_mu_list[0] << ((self.singlepass0.chip_select - 4) * 8)) |
                (self.att_singlepass_default_mu_list[1] << ((self.singlepass1.chip_select - 4) * 8))
        )
        self.att_reg_bichromatic = 0x00000000 | (
                (att_to_mu(self.att_doublepass_default_db * dB) << ((self.qubit.beam.chip_select - 4) * 8)) |
                (atts_cat_mu[0] << ((self.singlepass0.chip_select - 4) * 8)) |
                (atts_cat_mu[1] << ((self.singlepass1.chip_select - 4) * 8))
        )

        '''CREATE EXPERIMENT CONFIG'''
        # create sampling grid in radial coordinates
        vals_char_mu_pow_list = np.array([
            [
                self.core.seconds_to_mu(math.sqrt(x_us ** 2. + y_us ** 2.) * us),
                self.singlepass0.turns_to_pow(np.arctan2(y_us, x_us) / (2. * np.pi))
            ]
            for x_us in self.time_char_cat_x_us_list
            for y_us in self.time_char_cat_y_us_list
        ], dtype=np.int64)

        # heralding values
        self.time_force_herald_slack_mu = self.core.seconds_to_mu(150 * us)

        '''
        CREATE EXPERIMENT CONFIG
        '''
        # use generator to flatten list with some tuples
        def flatten(xs):
            for x in xs:
                if isinstance(x, Iterable) and not isinstance(x, (str, bytes)):
                    yield from flatten(x)
                else:
                    yield x

        # create an array of values for the experiment to sweep
        self.config_experiment_list = np.array([
            list(flatten(vals))
            for vals in product(
                freq_cat_center_ftw_list, freq_cat_secular_ftw_list,
                vals_char_mu_pow_list,
                characteristic_axis_list
            )
        ], dtype=np.int64)
        np.random.shuffle(self.config_experiment_list)

        # # tmp remove
        # self.freq_bsb_ftw = self.qubit.frequency_to_ftw(101.5885 * MHz)
        # # tmp remove

    def _prepare_argument_checks(self) -> TNone:
        """
        Check experiment arguments for validity.
        """
        # ensure singlepass values are safe and valid
        if any((ampl_pct > self.max_ampl_singlepass_pct or ampl_pct < 0.
                for ampl_pct in self.ampl_singlepass_default_pct_list)):
            raise ValueError("Singlepass amplitude outside valid range - [0., {:f}].".format(self.max_ampl_singlepass_pct))
        if any((att_db > 31.5 or att_db < self.min_att_singlepass_db
                for att_db in self.att_singlepass_default_db_list)):
            raise ValueError("Singlepass attenuation outside valid range - [{:.1f}, 31.5].".format(self.min_att_singlepass_db))

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_experiment_list),
                6)


    '''
    MAIN SEQUENCE
    '''
    @kernel(flags={"fast-math"})
    def initialize_experiment(self) -> TNone:
        # ensure phase_autoclear disabled on all beams to prevent phase accumulator reset
        # enable RAM mode and clear DDS phase accumulator
        self.qubit.set_cfr1()
        self.singlepass0.set_cfr1()
        self.singlepass1.set_cfr1()
        self.qubit.cpld.io_update.pulse_mu(8)
        delay_mu(25000)

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

        self.singlepass0.set_att_mu(self.att_singlepass_default_mu_list[0])
        self.singlepass1.set_att_mu(self.att_singlepass_default_mu_list[1])
        self.singlepass0.sw.on()
        self.singlepass1.sw.off()
        delay_mu(25000)

        # record general subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.sidebandcool_subsequence.record_dma()
        self.readout_subsequence.record_dma()
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        # instantiate relevant variables
        time_start_mu = now_mu() & ~0x7
        ion_state = (-1, 0, np.int64(0))

        for trial_num in range(self.repetitions):
            self.core.break_realtime()

            # sweep exp config
            for config_vals in self.config_experiment_list:

                '''PREPARE & CONFIGURE'''
                # extract values from config list
                freq_cat_center_ftw =       np.int32(config_vals[0])
                freq_cat_secular_ftw =      np.int32(config_vals[1])
                time_char_cat_mu =          config_vals[2]
                phase_char_cat_pow =        np.int32(config_vals[3])
                characteristic_axis_bool =  bool(config_vals[4])

                # prepare variables for execution
                char_read_phases = [
                    self.phases_char_cat_pow[0] + self.phases_char_cat_update_dir[0] * phase_char_cat_pow,
                    self.phases_char_cat_pow[1] + self.phases_char_cat_update_dir[1] * phase_char_cat_pow,
                ]
                self.core.break_realtime()

                '''PREPARE MOTIONAL STATE'''
                while True:
                    self.core.break_realtime()

                    '''INITIALIZE'''
                    # initialize ion in S-1/2 state & SBC to ground state
                    self.initialize_subsequence.run_dma()
                    self.sidebandcool_subsequence.run_dma()

                    # set target profile to ensure we run correctly
                    self.qubit.set_profile(self.profile_729_target)
                    self.qubit.cpld.io_update.pulse_mu(8)

                    # synchronize start time to coarse RTIO clock
                    time_start_mu = now_mu() & ~0x7

                    '''APPLY MOTIONAL INTERACTION'''
                    # sigma_x: select cat vs coherent state
                    if self.enable_motion_sigmax:
                        self.pulse_sigmax(time_start_mu, 0, True)

                    # bichromatic interaction to generate motional state
                    if self.enable_motion_cat:
                        self.pulse_bichromatic(time_start_mu, self.time_motion_cat_mu,
                                               self.phases_motion_cat_pow,
                                               freq_cat_center_ftw, freq_cat_secular_ftw)

                    # herald ion via state-dependent fluorescence (to projectively disentangle spin/motion)
                    if self.enable_motion_herald:
                        ion_state = self.readout_adaptive_subsequence.run()
                        delay_mu(20000)
                        self.pump.off()

                        # ensure dark state (flag is 0)
                        if ion_state[0] != 0: continue
                        # otherwise, add minor slack and proceed
                        at_mu(self.core.get_rtio_counter_mu() + self.time_force_herald_slack_mu)

                    # force break loop by default
                    break

                '''DIRECT CHARACTERISTIC MEASUREMENT'''
                # prepare spin state for characteristic readout
                # note: need to set correct profile for normal quenching
                # otherwise might be stuck in SBC quench params)
                self.pump.readout()
                self.repump_qubit.on()
                delay_mu(self.initialize_subsequence.time_repump_qubit_mu)
                self.repump_qubit.off()

                # sigma_x to select axis (does dummy if characteristic_axis_bool is False)
                self.pulse_sigmax(time_start_mu, self.phase_char_axis_pow, characteristic_axis_bool)

                # char read: bichromatic
                self.pulse_bichromatic(time_start_mu, time_char_cat_mu,
                                       char_read_phases,
                                       freq_cat_center_ftw, freq_cat_secular_ftw)

                '''READOUT'''
                # read out fluorescence & save results
                self.readout_subsequence.run_dma()
                counts_res = self.readout_subsequence.fetch_count()
                self.update_results(freq_cat_center_ftw,
                                    counts_res,
                                    freq_cat_secular_ftw,
                                    time_char_cat_mu,
                                    phase_char_cat_pow,
                                    characteristic_axis_bool)
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

        self.singlepass0.set_att_mu(self.att_singlepass_default_mu_list[0])
        self.singlepass1.set_att_mu(self.att_singlepass_default_mu_list[1])
        self.singlepass0.sw.on()
        self.singlepass1.sw.off()
        delay_mu(25000)


    '''
    HELPER FUNCTIONS
    '''
    @kernel(flags={"fast-math"})
    def pulse_sigmax(self, time_start_mu: TInt64 = -1, phas_pow: TInt32 = 0x0, is_real: TBool = False) -> TNone:
        """
        Run a phase-coherent sigma_x pulse on the qubit.
        Arguments:
            time_start_mu: fiducial timestamp for initial start reference (in machine units).
            phas_pow: relative phase offset for the beam.
            is_real: whether to actually run the pulse (True) or a dummy pulse (False).
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
        self.qubit.cpld.set_all_att_mu(self.att_reg_sigmax)

        # run pulse
        self.singlepass0.sw.on()
        self.singlepass1.sw.on()
        if is_real:
            self.qubit.on()
        else:
            self.qubit.off()
        delay_mu(self.time_sigmax_mu)
        self.qubit.off()
        self.singlepass1.sw.off()

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
            freq_carrier_ftw, asf=self.ampl_doublepass_default_asf, pow_=0,
            profile=self.profile_729_target, phase_mode=ad9910.PHASE_MODE_TRACKING, ref_time_mu=time_start_mu
        )
        self.singlepass0.set_mu(
            self.freq_singlepass_default_ftw_list[0]-freq_secular_ftw, asf=self.ampls_cat_asf[0],
            pow_=phas_pow_list[0], profile=self.profile_729_target,
            phase_mode=ad9910.PHASE_MODE_TRACKING, ref_time_mu=time_start_mu
        )
        self.singlepass1.set_mu(
            self.freq_singlepass_default_ftw_list[1]+freq_secular_ftw, asf=self.ampls_cat_asf[1],
            pow_=phas_pow_list[1], profile=self.profile_729_target,
            phase_mode=ad9910.PHASE_MODE_TRACKING, ref_time_mu=time_start_mu
        )
        self.qubit.cpld.set_all_att_mu(self.att_reg_bichromatic)

        # run bichromatic pulse
        self.singlepass0.sw.on()
        self.singlepass1.sw.on()
        self.qubit.on()
        delay_mu(time_pulse_mu)
        self.qubit.off()
        self.singlepass1.sw.off()

