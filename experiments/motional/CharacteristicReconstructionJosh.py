import numpy as np
from artiq.experiment import *

from LAX_exp.analysis import *
from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import (
    InitializeQubit, SidebandCoolContinuousRAM, Readout, RescueIon
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
        # hardware objects
        'initialize_subsequence', 'sidebandcool_subsequence', 'readout_subsequence', 'rescue_subsequence',
        'singlepass0', 'singlepass1',

        # hardware parameters
        'freq_singlepass_default_ftw_list', 'ampl_singlepass_default_asf_list', 'att_singlepass_default_mu_list',
        'ampl_doublepass_default_asf', 'att_doublepass_default_mu',
        'freq_sigmax_ftw', 'ampl_sigmax_asf', 'att_sigmax_mu', 'time_sigmax_mu',
        'time_force_herald_slack_mu',

        # cat-state parameters
        'ampls_cat_asf', 'atts_cat_mu', 'time_pulse1_cat_mu', 'phases_pulse1_cat_pow', 'phase_characteristic_axis_pow',
        'phases_pulse5_cat_pow', 'phases_pulse5_cat_update_dir',

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
        self.setattr_argument("freq_sigmax_mhz",    NumberValue(default=101.3341, precision=6, step=1, min=50., max=400.), group="defaults.sigmax")
        self.setattr_argument("ampl_sigmax_pct",    NumberValue(default=50., precision=3, step=5, min=0.01, max=50), group="defaults.sigmax")
        self.setattr_argument("att_sigmax_db",      NumberValue(default=8., precision=1, step=0.5, min=8., max=31.5), group="defaults.sigmax")
        self.setattr_argument("time_sigmax_us",     NumberValue(default=3.05, precision=2, step=5, min=0.1, max=10000), group="defaults.sigmax")

        # defaults - cat
        self.setattr_argument("freq_cat_center_mhz_list",   Scannable(
                                                                default=[
                                                                    ExplicitScan([101.3341]),
                                                                    CenterScan(101.3341, 0.01, 0.0001, randomize=True),
                                                                    RangeScan(101.3200, 101.3500, 50, randomize=True),
                                                                ],
                                                                global_min=60., global_max=400, global_step=1,
                                                                unit="MHz", scale=1, precision=6
                                                            ), group='defaults.cat')
        self.setattr_argument("freq_cat_secular_khz_list",  Scannable(
                                                                default=[
                                                                    ExplicitScan([701.8]),
                                                                    CenterScan(701.8, 4, 0.1, randomize=True),
                                                                    RangeScan(699.0, 704.2, 50, randomize=True),
                                                                ],
                                                                global_min=0, global_max=10000, global_step=1,
                                                                unit="kHz", scale=1, precision=3
                                                            ), group='defaults.cat')
        self.setattr_argument("ampls_cat_pct",  PYONValue([58., 58.]), group='defaults.cat', tooltip="[rsb_pct, bsb_pct]")
        self.setattr_argument("atts_cat_db",    PYONValue([6., 6.]), group='defaults.cat', tooltip="[rsb_db, bsb_db]")

        '''PULSE ARGUMENTS - MOTIONAL STATE PREPARATION'''
        # pulse 0 - sigma_x0: select motional state type (cat state vs coherent state)
        self.setattr_argument("enable_pulse0_sigmax", BooleanValue(default=False), group='pulse0.sigmax',
                              tooltip='sigma_x (pulse #0) selects whether motional state is coherent (True), or cat (False)')

        # pulse 1 - cat 0
        self.setattr_argument("enable_pulse1_cat",          BooleanValue(default=False), group='pulse1.cat')
        self.setattr_argument("time_pulse1_cat_us",         NumberValue(default=100, precision=2, step=5, min=0.1, max=10000), group="pulse1.cat")
        self.setattr_argument("phases_pulse1_cat_turns",    PYONValue([0., 0.]), group='pulse1.cat', tooltip="[rsb_turns, bsb_turns]")

        #pulse 2 sigmax
        self.setattr_argument("enable_pulse2_sigmax", BooleanValue(default=False), group='pulse2.sigmax',
                              tooltip='sigma_x (pulse #2)')

        # pulse 3 - herald
        self.setattr_argument("enable_pulse3_herald",       BooleanValue(default=True), group='pulse3.herald')
        self.setattr_argument("enable_force_herald",        BooleanValue(default=True), group='pulse3.herald')
        self.setattr_argument("force_herald_threshold",     NumberValue(default=46, precision=0, step=10, min=0, max=10000), group='pulse3.herald')

        '''PULSE ARGUMENTS - CHARACTERISTIC FUNCTION MEASUREMENT'''
        # pulse 4 - sigma_x1: select real/imag part of characteristic function
        self.setattr_argument("characteristic_axis", EnumerationValue(['Both', 'Real', 'Imaginary'], default='Both'), group='characteristic_reconstruction',
                              tooltip="Selects the real/imag component of the characteristic function by either applying a sigma_x operation (Imag), or not (Real)."
                                      "The 'Both' option enables measurement of both real and imag components within a single experiment.")
        self.setattr_argument("phase_characteristic_axis_turns", NumberValue(default=0., precision=3, step=0.1, min=-1.0, max=1.0), group='characteristic_reconstruction',
                              tooltip="Sets the relative phase of the sigma_x operation used to define the real/imag axis of the characteristic function.")

        # pulse 5 - cat 1: characteristic readout protocol
        self.setattr_argument("enable_pulse5_cat",          BooleanValue(default=True), group='characteristic_reconstruction')
        self.setattr_argument("phases_pulse5_cat_turns",    PYONValue([0., 0.]), group='characteristic_reconstruction', tooltip="[rsb_turns, bsb_turns]")
        self.setattr_argument("target_pulse5_cat_phase",    EnumerationValue(['RSB', 'BSB', 'RSB-BSB', 'RSB+BSB'], default='RSB-BSB'), group="characteristic_reconstruction")
        self.setattr_argument("time_pulse5_cat_x_us_list",      Scannable(
                                                                default=[
                                                                    RangeScan(-500, 500, 10, randomize=True),
                                                                    ExplicitScan([100]),
                                                                ],
                                                                global_min=-100000, global_max=100000, global_step=1,
                                                                unit="us", scale=1, precision=5
                                                            ), group="characteristic_reconstruction")
        self.setattr_argument("time_pulse5_cat_y_us_list",      Scannable(
                                                                default=[
                                                                    RangeScan(-500, 500, 10, randomize=True),
                                                                    ExplicitScan([100]),
                                                                ],
                                                                global_min=-100000, global_max=100000, global_step=1,
                                                                unit="us", scale=1, precision=5
                                                            ), group="characteristic_reconstruction")

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

        # pulse 4 - define characteristic axis (via sigma_x)
        self.phase_characteristic_axis_pow =  self.qubit.turns_to_pow(self.phase_characteristic_axis_turns)
        # configure whether real/imag/both axes of the characteristic function are to be measured
        if self.characteristic_axis == "Real":          characteristic_axis_list = [False]
        elif self.characteristic_axis == "Imaginary":   characteristic_axis_list = [True]
        elif self.characteristic_axis == "Both":        characteristic_axis_list = [True, False]

        # pulse 5 - cat 1
        self.phases_pulse5_cat_pow =    np.array([self.singlepass0.turns_to_pow(phas_pow)
                                                  for phas_pow in self.phases_pulse5_cat_turns], dtype=np.int32)

        if self.target_pulse5_cat_phase == 'RSB':
            self.phases_pulse5_cat_update_dir = np.array([1, 0], dtype=np.int32)
        elif self.target_pulse5_cat_phase == 'BSB':
            self.phases_pulse5_cat_update_dir = np.array([0, 1], dtype=np.int32)
        elif self.target_pulse5_cat_phase == 'RSB-BSB':
            self.phases_pulse5_cat_update_dir = np.array([1, -1], dtype=np.int32)
        elif self.target_pulse5_cat_phase == 'RSB+BSB':
            self.phases_pulse5_cat_update_dir = np.array([1, 1], dtype=np.int32)

        # create sampling grid in radial coordinates
        if self.enable_pulse5_cat:
            vals_pulse5_mu_pow_list = np.array([
                [
                    self.core.seconds_to_mu(math.sqrt(x_us ** 2. + y_us ** 2.) * us),
                    self.singlepass0.turns_to_pow(np.arctan2(y_us, x_us) / (2. * np.pi))
                ]
                for x_us in self.time_pulse5_cat_x_us_list
                for y_us in self.time_pulse5_cat_y_us_list
            ], dtype=np.int64)
        else:
            vals_pulse5_mu_pow_list = np.array([[0, 0]], dtype=np.int64)

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
                vals_pulse5_mu_pow_list,
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
        # ensure we only herald once
        if self.enable_force_herald and not self.enable_pulse3_herald:
            raise ValueError("Cannot force_herald if enable_pulse3_herald is disabled. Check input arguments.")

        # ensure single pass values are safe and valid
        if any((ampl_pct > self.max_ampl_singlepass_pct or ampl_pct < 0.
                for ampl_pct in self.ampl_singlepass_default_pct_list)):
            raise ValueError("Singlepass amplitude outside valid range - [0., {:f}].".format(self.max_ampl_singlepass_pct))
        if any((att_db > 31.5 or att_db < self.min_att_singlepass_db
                for att_db in self.att_singlepass_default_db_list)):
            raise ValueError("Singlepass attenuation outside valid range - [{:.1f}, 31.5].".format(self.min_att_singlepass_db))

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_experiment_list),
                7)


    '''
    MAIN SEQUENCE
    '''
    @kernel(flags={"fast-math"})
    def initialize_experiment(self) -> TNone:
        self.core.break_realtime()

        # set up qubit beam for DMA sequences
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
                time_pulse5_cat_mu =        config_vals[2]
                phase_pulse5_cat_pow =      np.int32(config_vals[3])
                characteristic_axis_bool =  bool(config_vals[4])
                self.core.break_realtime()

                # prepare variables for execution
                counts_her = -1
                time_start_mu = now_mu() & ~0x7
                cat5_phases = [
                    self.phases_pulse5_cat_pow[0] + self.phases_pulse5_cat_update_dir[0] * phase_pulse5_cat_pow,
                    self.phases_pulse5_cat_pow[1] + self.phases_pulse5_cat_update_dir[1] * phase_pulse5_cat_pow,
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
                    # pulse 0 (sigma_x #1): select cat vs coherent state
                    if self.enable_pulse0_sigmax:
                        self.pulse_sigmax(time_start_mu, 0)

                    # pulse 1 (cat 1)
                    if self.enable_pulse1_cat:
                        self.pulse_bichromatic(time_start_mu, self.time_pulse1_cat_mu,
                                               self.phases_pulse1_cat_pow,
                                               freq_cat_center_ftw, freq_cat_secular_ftw)

                    # tmp remove
                    if self.enable_pulse2_sigmax:
                        self.pulse_sigmax2(time_start_mu, 0)
                    # tmp remove

                    # pulse 3: herald via 397nm fluorescence
                    if self.enable_pulse3_herald:
                        self.readout_subsequence.run_dma()
                        self.pump.off()

                        # pulse 3: force heralding
                        if self.enable_force_herald:
                            counts_her = self.readout_subsequence.fetch_count()
                            if counts_her > self.force_herald_threshold:
                                continue

                            # add minor slack if we proceed
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


                # pulse 4: sigma_x #2
                if characteristic_axis_bool:
                    self.pulse_sigmax(time_start_mu, self.phase_characteristic_axis_pow)

                # pulse 5: cat #2
                if self.enable_pulse5_cat:
                    self.pulse_bichromatic(time_start_mu, time_pulse5_cat_mu,
                                           cat5_phases,
                                           freq_cat_center_ftw, freq_cat_secular_ftw)

                '''READOUT'''
                # read out fluorescence
                self.readout_subsequence.run_dma()

                # retrieve heralded measurement
                if self.enable_pulse3_herald and not self.enable_force_herald:
                    counts_her = self.readout_subsequence.fetch_count()

                # retrieve results & update dataset
                counts_res = self.readout_subsequence.fetch_count()
                self.update_results(freq_cat_center_ftw,
                                    counts_res,
                                    counts_her,
                                    freq_cat_secular_ftw,
                                    time_pulse5_cat_mu,
                                    phase_pulse5_cat_pow,
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
                (self.att_singlepass_default_mu_list[1] << (2 * 8))
        )
        self.qubit.cpld.set_all_att_mu(a)

        # run sigmax pulse
        self.singlepass0.sw.on()
        self.singlepass1.sw.off()
        self.qubit.on()
        delay_mu(self.time_sigmax_mu)
        self.qubit.off()

    @kernel(flags={"fast-math"})
    def pulse_sigmax2(self, time_start_mu: TInt64, phas_pow: TInt32) -> TNone:
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
                (self.att_singlepass_default_mu_list[1] << (2 * 8))
        )
        self.qubit.cpld.set_all_att_mu(a)

        # run sigmax pulse
        self.singlepass0.sw.on()
        self.singlepass1.sw.off()
        self.qubit.on()
        delay_mu(2800)
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
            self.freq_singlepass_default_ftw_list[0]-freq_secular_ftw, asf=self.ampls_cat_asf[0],
            pow_=phas_pow_list[0], profile=self.profile_729_target,
            phase_mode=ad9910.PHASE_MODE_TRACKING, ref_time_mu=time_start_mu
        )
        self.singlepass1.set_mu(
            self.freq_singlepass_default_ftw_list[1]+freq_secular_ftw, asf=self.ampls_cat_asf[1],
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

    # @kernel(flags={"fast-math"})
    # def pulse_yzde(self) -> TNone:
    #     # set up relevant beam waveforms
    #     self.qubit.set_mu(
    #         self.freq_bsb_ftw, asf=self.ampl_sigmax_asf, pow_=0,
    #         profile=self.profile_729_target, phase_mode=ad9910.PHASE_MODE_ABSOLUTE
    #     )
    #     self.singlepass0.set_mu(
    #         self.freq_singlepass_default_ftw_list[0], asf=self.ampl_singlepass_default_asf_list[0], pow_=0,
    #         profile=self.profile_729_target, phase_mode=ad9910.PHASE_MODE_ABSOLUTE
    #     )
    #     self.singlepass1.set_mu(
    #         self.freq_singlepass_default_ftw_list[1], asf=self.ampl_singlepass_default_asf_list[1], pow_=0,
    #         profile=self.profile_729_target, phase_mode=ad9910.PHASE_MODE_ABSOLUTE
    #     )
    #     self.qubit.cpld.io_update.pulse_mu(8)
    #
    #     # set all attenuators together
    #     a = self.qubit.cpld.att_reg & ~(
    #             (0xFF << (0 * 8)) |
    #             (0xFF << (1 * 8)) |
    #             (0xFF << (2 * 8))
    #     )
    #     a |= (
    #             (self.att_sigmax_mu << (0 * 8)) |
    #             (self.att_singlepass_default_mu_list[0] << (1 * 8)) |
    #             (self.att_singlepass_default_mu_list[1] << (2 * 8))
    #     )
    #     self.qubit.cpld.set_all_att_mu(a)
    #
    #     # run sigmax pulse
    #     self.singlepass0.sw.on()
    #     self.singlepass1.sw.off()
    #     self.qubit.on()
    #     delay_mu(66980)
    #     self.qubit.off()


    '''
    ANALYSIS
    '''
    def analyze_experiment(self):
        pass

