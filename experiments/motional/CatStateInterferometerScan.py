import numpy as np
from LAX_exp.analysis.artiq_conversions import ftw_to_frequency_khz
from LAX_exp.analysis.processing import findThresholdScikit
from artiq.experiment import *
from artiq.coredevice import ad9910

from numpy import array, int32, int64, zeros
from numpy import shape as array_shape
from scipy.optimize import curve_fit

from LAX_exp.language import *
from LAX_exp.system.subsequences import (
    InitializeQubit, SidebandCoolContinuousRAM, Readout, ReadoutAdaptive,
    RescueIon,
)

from LAX_exp.system.objects.PulseShaper import available_pulse_shapes
from LAX_exp.system.objects.dds_pulse_shaper import DDSPulseShaper
from LAX_exp.system.objects.dds_ramper import DDSRamper



class CatStateInterferometerScan(LAXExperiment, Experiment):
    """
    Experiment: Cat State Interferometer Scam

    Create and characterize cat states with projective state preparation.
    Uses adaptive readout to reduce timing overheads and extend available coherence times.
    """
    name = 'Cat State Inteferometer Scan'
    kernel_invariants = {
        # subsequences & objects
        'initialize_subsequence', 'sidebandcool_subsequence', 'readout_subsequence',
        'rescue_subsequence',
    
        # hardware values - ms - bichromatic
        'enable_ms_gate', 'ampls_ms_asf',

        # hardware values - cat1 - bichromatic
        'time_cat_bichromatic_mu', 'phases_pulse1_cat_pow',

        # hardware values - tickle
        'time_tickle_mu',

        #   hardware values - dynamical decoupling
        'enable_dynamical_decoupling', 'ampl_dynamical_decoupling_asf',

        # hardware values - intensity servo
        'enable_servo_relock', 'time_servo_relock_mu',

        # configs
        'profile_729_SBC',
        'profile_729_bichromatic',
        'profile_tickle_RAM',
        'att_reg_ms_gate',
        'config_experiment_list',

        # extras
        'urukul_setup_time_mu', 'indices', 'urukul_dd_reset_time',
    }

    def build_experiment(self):

        # core arguments
        self.setattr_argument("repetitions", NumberValue(default=50, precision=0, step=1, min=1, max=100000))
        self.setattr_argument("enable_linetrigger", BooleanValue(default=False),
                              tooltip="Trigger the beginning of each shot from the AC line.")

        # allocate relevant beam profiles
        self.profile_729_RAP = 0
        self.profile_729_SBC = 1
        self.profile_729_bichromatic = 2

        self.index_729_cat1 = 0
        self.index_729_cat2 = 1
        self.index_729_ms = 2
        self.index_729_parity = 3

        # allocate profiles for dds tickle
        self.profile_tickle_RAM = 0

        # get subsequences
        self.sidebandcool_subsequence = SidebandCoolContinuousRAM(
            self, profile_729=self.profile_729_SBC, profile_854=3,
            ram_addr_start_729=0, ram_addr_start_854=0, num_samples=200
        )
        self.initialize_subsequence = InitializeQubit(self)
        self.readout_subsequence = Readout(self)
        self.rescue_subsequence = RescueIon(self)

        # relevant devices
        self.setattr_device('qubit')
        self.setattr_device('pump')
        self.setattr_device('repump_qubit')
        self.setattr_device('dds_dipole')
        self.setattr_device('trigger_line')

        # set build arguments
        self._build_arguments_ion_parameters()
        self._build_arguments_dynamical_decoupling()
        self._build_arguments_ms_gate()
        self._build_arguments_cat_default()
        self._build_arguments_cat_modes()
        self._build_arguments_tickle_default()
        self._build_arguments_tickle_modes()
        self._build_arguments_intensity_servo()

    def _build_arguments_ion_parameters(self):
        """
        Build arguments for ion frequencies
        """
        _argstr = "ion_parameters"

        self.setattr_argument("freq_carrier_mhz", NumberValue(
            default=100.4,
            min=60., max=400, step=1,
            unit="MHz", scale=1, precision=6
        ), group=_argstr,tooltip="Carrier frequency of the ion.\n"
                                      "Note: this is applied via the main doublepass DDS.\n")

        self.setattr_argument("freq_secular_khz_list", PYONValue([700, 1700]),
                             group=_argstr,
                             tooltip="Secular frequencies of tickle pulse (in kHz) applied via the urukul dds.")

    def _build_arguments_ms_gate(self):
        """
        Build arguments for ms gate beam parameters.
        """
        _argstr = "ms_gate"
        self.setattr_argument("enable_ms_gate", BooleanValue(default=True),
                              group=_argstr,
                              tooltip="Enables application of the ms gate before catting.")

        # ms gate - pulse parameters
        self.setattr_argument("ampls_ms_pct", PYONValue([50., 50.]), group=_argstr,
                              tooltip="DDS amplitudes for the singlepass DDSs during the ms gate.\n"
                                      "Should be a list of [rsb_ampl_pct, bsb_ampl_pct], which are applied to [singlepass1, singlepass2].")
        self.setattr_argument("atts_ms_db", PYONValue([14., 14.]), group=_argstr,
                              tooltip="DDS attenuations for the singlepass DDSs during the ms gate.\n"
                                      "Should be a list of [rsb_att_db, bsb_att_db], which are applied to [singlepass1, singlepass2].")

        # scanning options
        self.setattr_argument("time_ms_gate_us", NumberValue(
            default=100.,
            min=1, max=10000, step=1,
            unit="us", scale=1, precision=5
        ),
                              group=_argstr,
                              tooltip="Pulse time for the the ms gate")

        self.setattr_argument("freq_ms_mode_khz", NumberValue(
            default=700,
            min=500, max=3000, step=0.001,
            unit="kHz", scale=1, precision=6
        ), group=_argstr,
                              tooltip='frequency of mode used for MS gate')

        self.setattr_argument("freq_ms_secular_detuning_khz", NumberValue(
            default=10.,
            min=-100, max=10000, step=1,
            unit="kHz", scale=1, precision=3
        ), group=_argstr,
                             tooltip="Single-sided detuning from the secular frequency for the ms gate, applied via singlepass DDSs.\n"
                                     "The singlepass1 DDS is treated as the RSB, and will thus have its frequency from the secular DECREASED by this amount.\n"
                                     "Similarly, the singlepass2 DDS is treated as the BSB, and will thus have its frequency from the secular INCREASED by this amount.\n"
                                     "i.e. frequencies for [singlepass1, singlepass2] is set as [beams.freq_mhz.freq_singlepass1_mhz - freq_secular_khz - freq_ms_secular_detuning, "
                                     "beams.freq_mhz.freq_singlepass2_mhz + freq_secular_khz + freq_ms_secular_detuning].")

        self.setattr_argument("phase_ms_turns", PYONValue([0., 0.]), group=_argstr,
                              tooltip="Phase sweep values applied to the singlepass DDSs during the ms gate.\n")

        self.setattr_argument('phase_ms_dynamical_decoupling_turns',
                              NumberValue(default=0., unit='turns',
                                  min=0., max=2., step=0.1, precision=3, scale=1.0),
                              tooltip='Phase of the third rf tone applied to the 729 single pass AOM during the MS gate\n'
                                      'For dynamical decoupling to work the phase must be configured correctly as we \n'
                                      'need [H_{bi}, H_{dd}]=0 and this only occurs if these Hamiltonians contain \n'
                                      'the same linear combination of c1*sigma_{x}+c2*sigma_{y} which is determined '
                                      'by \n the laser phase.',
                              group=_argstr)

    def _build_arguments_cat_default(self):
        _argstr = 'default.cat'
        self.setattr_argument("time_cat_bichromatic_us",
                              NumberValue(default=50, precision=2, step=5, min=0.1, max=10000000, scale=1., unit="us"),
                              group=_argstr,
                              tooltip="Pulse time for the bichromatic pulses.")

        self.setattr_argument('phase_cat_dynamical_decoupling_turns',
                              NumberValue(default=0.0, unit='turns',
                                  min=0., max=2., step=0.1, precision=3, scale=1.0),
                              tooltip='Phase of the third rf tone applied to the 729 single pass AOM during the cats\n'
                                      'For dynamical decoupling to work the phase must be configured correctly as we \n'
                                      'need [H_{bi}, H_{dd}]=0 and this only occurs if these Hamiltonians contain \n'
                                      'the same linear combination of c1*sigma_{x}+c2*sigma_{y} which is determined '
                                      'by \n the laser phase.',
                              group=_argstr)

        self.setattr_argument("phases_pulse1_cat_turns", PYONValue([0., 0.]), group=_argstr,
                              tooltip="Relative phases for the singlepass DDSs during the 1st bichromatic pulse.\n"
                                      "Should be a list of [rsb_phase_turns, bsb_phase_turns].\n"
                                      "Note: these phases are applied to the singlepass DDSs, so do not need to be halved, "
                                      "unlike phases applied to the main doublepass.")

    def _build_arguments_tickle_default(self):
        _argstr = 'default.tickle'

        self.setattr_argument("time_tickle_us",
                              NumberValue(default=50, precision=2, step=500, min=0.04, max=10000000, unit="us",
                                          scale=1.),
                              group=_argstr,
                              tooltip="Time for the total pulse (including pulse shape).")

        # waveform - pulse shaping
        self.setattr_argument("enable_pulse_shaping", BooleanValue(default=False),
                              group=_argstr,
                              tooltip="Applies pulse shaping to the edges of the tickle pulse.")
        self.setattr_argument("type_pulse_shape",
                              EnumerationValue(list(available_pulse_shapes.keys()), default='sine_squared'),
                              group=_argstr,
                              tooltip="Pulse shape type to be used.")


    def _build_arguments_cat_modes(self):
        """
        Build arguments for default cat beam parameters.
        """
        _argstr = 'cat.mode_parameters'
        self.setattr_argument("ampls_cat_mode_pct", PYONValue([[50., 50.], [50., 50.]]), group=_argstr,
                              tooltip="DDS amplitudes for the singlepass DDSs during the bichromatic pulses.\n"
                                      "Should be a list of [rsb_ampl_pct, bsb_ampl_pct], which are applied to [singlepass1, singlepass2].")
        self.setattr_argument("atts_cat_mode_db", PYONValue([[11., 11.], [11., 11.]]), group=_argstr,
                              tooltip="DDS attenuations for the singlepass DDSs during the bichromatic pulses.\n"
                                      "Should be a list of [rsb_att_db, bsb_att_db], which are applied to [singlepass1, singlepass2].")


    def _build_arguments_tickle_modes(self):
        """
        Build core sweep arguments for the tickle pulse.
        """
        _argstr = "tickle.mode_parameters"  # string to use for arguments

        # waveform - parameter sweeps
        self.setattr_argument("att_tickle_modes_db",
                              PYONValue([31.5, 31.5]),
                              group=_argstr,
                              tooltip="Attenuation to be used for the urukul channel used for generating the tickle.")
        self.setattr_argument("ampl_tickle_modes_pct",
                              PYONValue([1., 1.]),
                              group=_argstr,
                              tooltip='Amplitude of tickle pulse.')

        self.setattr_argument("phase_tickle_turns", PYONValue([0., 0.]),
                              group=_argstr,
                              tooltip="Phase of tickle pulse (in turns) applied via the urukul dds.")

        self.setattr_argument("freq_tickle_detunings_mode_khz_list", Scannable(
                              default=[
                                  ExplicitScan([0]),
                                  CenterScan(0, 4, 0.1, randomize=True),
                                  RangeScan(-10, 10, 50, randomize=True),
                              ],
                              global_max = 2000, global_min = -2000,
                                global_step=0.01, precision=3, scale=1.0),
                              group=_argstr,
                              tooltip="Detuning from secular frequency of tickle pulse (in kHz) applied via the urukul dds.")

    def _build_arguments_dynamical_decoupling(self):
        """
        Get arguments for dynamical decoupling
        """
        _argstr = 'dynamical_decoupling'
        self.setattr_argument('enable_dynamical_decoupling', BooleanValue(True),
                              tooltip='Indicate whether to apply a third rf tone on the 729 single pass AOM. This '
                                      'tone \n'
                                      'will produce a linear combination of sigma_x and sigma_y \n'
                                      '(dependent on laser phase) which will rotate away sigma_z errors \n',
                              group=_argstr)

        self.setattr_argument('ampl_dynamical_decoupling_pct', NumberValue(default=50., max=50., min=0.01,
                                                                           step=5, precision=2, scale=1., unit='%'),
                              tooltip='Amplitude of third rf tone applied to the 729 single pass AOM. \n'
                                      'Stronger tones will further weaken any sigma_z errors BUT \n'
                                      'NOTE: Changing the strength will change AC stark shifts to the sidebands.',
                              group=_argstr)

        self.setattr_argument('att_dynamical_decoupling_dB', NumberValue(default=28., max=31.5, min=5.,
                                                                         step=0.5, precision=1, scale=1., unit='dB'),
                              tooltip='Attenuation of the third rf tone applied to the 729 single pass AOM.',
                              group=_argstr)

    def _build_arguments_intensity_servo(self):
        """
        Build arguements for intensity servo hold between shots
        """
        _argstr = 'intensity_servo_relock'
        self.setattr_argument('enable_servo_relock', BooleanValue(default=False), group=_argstr,
                              tooltip='Enables the servo to relock the intensity of the 729 beam after every shot')
        self.setattr_argument('time_servo_relock_us', NumberValue(default=2000, precision=3, step=1, min=1,
                                                                  max=10000, scale=1., unit='us'),
                              group=_argstr, tooltip='Length of time to let the servo relock before each shot')

    def prepare_experiment(self):
        """
        Prepare & precompute experimental values.
        """
        self._prepare_argument_checks()
        self.time_urukul_reset_mu = int64(1330) # half the time for urukul to implement set_mu for DD phase shifting

        ### Build Pulse Shaper
        if self.enable_pulse_shaping:
            self.pulse_shape = self.type_pulse_shape
        else:
            self.pulse_shape = 'square'

        ### MAGIC NUMBERS ###
        self.urukul_setup_time_mu = int64(8)  # extra delay between calls to ad9910.set_mu
        self.urukul_dd_reset_time = int64(650) # half the time for urukul to implement set_mu for DD phase shifting

        # run component preparation
        self._prepare_experiment_dynamical_decoupling()
        self._prepare_experiment_cat_default()
        self._prepare_experiment_cat_modes()
        self._prepare_experiment_tickle_default()
        freq_detuning_tickle_ftw_list = self._prepare_experiment_tickle_modes()
        self._prepare_experiment_ion_parameters()
        self._prepare_experiment_ms_gate()

        self.phase_dd_phase_shift_pow = self.qubit.turns_to_pow(0.5)

        num_secular_freqs = len(self.freq_secular_khz_list)

        secular_freq_idx_list = range(num_secular_freqs)

        # create experiment config
        self.config_experiment_list = create_experiment_config(
            # tickle sweeps
            secular_freq_idx_list,
            freq_detuning_tickle_ftw_list,
            config_type=float, shuffle_config=True
        )

        # create index list for later use
        self.indices = [self.index_729_cat1, self.index_729_cat2, self.index_729_ms]

        # placeholder arrays for urukul values (first index is  profile number, second index is channel on urukul 0)
        self.phase_beams_pow_list = zeros((8, 4), dtype=int32)
        self.freq_beams_ftw_list = zeros((8, 4), dtype=int32)
        self.ampl_beams_asf_list = zeros((8, 4), dtype=int32)

        # set timing for intensity servo between shots
        self.time_servo_relock_mu = self.core.seconds_to_mu(self.time_servo_relock_us * us)

        self.phase_cat2_offset = self.qubit.turns_to_pow(0.5)
        self.set_default_configuration()


    def _prepare_experiment_ion_parameters(self):
        """
        Prepare general ion parameters
        :return: list of carrier frequencies for cat and MS
        """
        self.freq_carrier_ftw = self.qubit.frequency_to_ftw(self.freq_carrier_mhz * MHz)

        self.freq_secular_ftw_list = array([self.dds_pulse_shaper_tickle_list[idx].dds_target.frequency_to_ftw(self.freq_secular_khz_list[idx]*kHz)
         for  idx in range(len(self.freq_secular_khz_list))])


    def _prepare_experiment_ms_gate(self):
        """
        Prepare experimental values for ms gate

        :return: tuple of ms gate time, detunings from secular frequency, and ms beam phases
        """

        """
        Build Pulse Shaper
        """
        self.dds_ramper_ms = DDSRamper(self, self.qubit.singlepass1,
                                       num_samples=200, # number of points in the ramp
                                       ramp_dest=2, # amplitude ramp
                                       data_high=self.ampls_ms_pct[0] / 100., # max value of the ramp
                                       data_low=0) # min value of the ramp

        self.dds_ramper_ms.add_dds_target(self.qubit.singlepass2,
                                          ramp_dest=2, # ampltiude ramp
                                          data_high=self.ampls_ms_pct[1] / 100., # max value of the ramp
                                          data_low=0) # min value of the ramp

        self.time_ms_gate_mu = self.core.seconds_to_mu(self.time_ms_gate_us/2*us)
        self.freq_ms_mode_ftw = self.qubit.frequency_to_ftw(self.freq_ms_mode_khz * kHz)
        self.freq_ms_gate_secular_detuning_ftw = self.qubit.frequency_to_ftw(self.freq_ms_secular_detuning_khz * kHz)
        self.phase_ms_pow = array(
            [self.qubit.turns_to_pow(phase_ms_turns) for phase_ms_turns in self.phase_ms_turns])
        self.phase_ms_dynamical_decoupling_pow_list = self.qubit.singlepass0.turns_to_pow(
            self.phase_ms_dynamical_decoupling_turns)


        self.phase_ms_dynamical_decoupling_pow = self.qubit.singlepass0.turns_to_pow(self.phase_ms_dynamical_decoupling_turns)

        self.ampls_ms_asf = array([self.qubit.amplitude_to_asf(ampl_ms_pct/100.) for ampl_ms_pct in self.ampls_ms_pct])


        self.phase_ms_update_dir = array([1, 1], dtype=int32)

        self.att_reg_ms_gate = 0x00000000 | (
                (self.qubit.att_qubit_mu << ((self.qubit.beam.chip_select - 4) * 8)) |
                (self.att_dynamical_decoupling_mu << ((self.qubit.singlepass0.chip_select - 4) * 8)) |
                (att_to_mu(self.atts_ms_db[0] * dB) << ((self.qubit.singlepass1.chip_select - 4) * 8)) |
                (att_to_mu(self.atts_ms_db[1] * dB) << ((self.qubit.singlepass2.chip_select - 4) * 8))
        )

    def _prepare_experiment_dynamical_decoupling(self):
        """
        Prepare experiment for dynamical decoupling.
        """
        self.ampl_dynamical_decoupling_asf = self.qubit.singlepass0.amplitude_to_asf(
            self.ampl_dynamical_decoupling_pct / 100.)
        self.att_dynamical_decoupling_mu = self.qubit.singlepass0.cpld.att_to_mu(
            self.att_dynamical_decoupling_dB * dB)

    def _prepare_experiment_cat_default(self):
        """
        Prepare experiment values for cat/bichromatic.
        :return:
        """
        '''CONVERT VALUES TO MACHINE UNITS - BICHROMATIC/CAT DEFAULTS '''
        # defaults - main doublepass (near chamber)
        self.time_cat_bichromatic_mu = self.core.seconds_to_mu(self.time_cat_bichromatic_us/2 * us)
        self.phases_pulse1_cat_pow = [self.qubit.singlepass0.turns_to_pow(phas_pow)
                                      for phas_pow in self.phases_pulse1_cat_turns]

        '''CONVERT VALUES TO MACHINE UNITS - DD PHASE'''
        self.phase_cat_dynamical_decoupling_pow = self.qubit.singlepass0.turns_to_pow(self.phase_cat_dynamical_decoupling_turns)

        # specify phase update array based on user arguments
        self.phase_cat_update_dir = array([1, -1], dtype=int32)


    def _prepare_experiment_cat_modes(self):
        """
        Convert cat parameters for each mode and store in a list
        """
        self.ampls_cat_asf_list = [0]*len(self.ampls_cat_mode_pct)
        for idx in range(len(self.ampls_cat_asf_list)):
            self.ampls_cat_asf_list[idx] = array([self.qubit.singlepass0.amplitude_to_asf(ampl_pct / 100.)
                                    for ampl_pct in self.ampls_cat_mode_pct[idx]])

        self.atts_cat_reg_mu_list = [0]*len(self.atts_cat_mode_db)
        for idx in range(len(self.atts_cat_mode_db)):
            att_reg_cat_interferometer = 0x00000000 | (
                    (self.qubit.att_qubit_mu << ((self.qubit.beam.chip_select - 4) * 8)) |
                    (self.att_dynamical_decoupling_mu << ((self.qubit.singlepass0.chip_select - 4) * 8)) |
                    (att_to_mu(self.atts_cat_mode_db[idx][0] * dB) << ((self.qubit.singlepass1.chip_select - 4) * 8)) |
                    (att_to_mu(self.atts_cat_mode_db[idx][1] * dB) << ((self.qubit.singlepass2.chip_select - 4) * 8))
            )
            self.atts_cat_reg_mu_list[idx] = att_reg_cat_interferometer
        '''
        CREATE ATTENUATION REGISTERS
        '''

    def _prepare_experiment_tickle_default(self):
        """
        Prepare general experiment values for the tickle pulse.
        """
        # convert values to convenience units
        self.time_tickle_mu = self.core.seconds_to_mu(self.time_tickle_us * us)
        self.time_tickle_dd_phase_flip_mu = self.core.seconds_to_mu(self.time_tickle_us * us/2)

    def _prepare_experiment_tickle_modes(self):
        """
        Prepare general experiment values for the tickle pulse for each individual mode.
        """
        self.dds_pulse_shaper_tickle_list = [0]*len(self.ampl_tickle_modes_pct)
        for idx in range(len(self.ampl_tickle_modes_pct)):
            dds_pulse_shaper_tickle = DDSPulseShaper(self, dds_target=self.dds_dipole.dds,
                           ram_profile=self.profile_tickle_RAM,
                           ram_addr_start=202, num_samples=250,
                           ampl_max_pct=self.ampl_tickle_modes_pct[idx],
                           pulse_shape=self.pulse_shape,
                           phase_autoclear=1)
            self.dds_pulse_shaper_tickle_list[idx] = dds_pulse_shaper_tickle

        self.att_tickle_mu_list = [att_to_mu(self.att_tickle_modes_db[idx] * dB)
                                   for idx in range(len(self.att_tickle_modes_db))]
        self.phase_tickle_pow_list = array([self.dds_pulse_shaper_tickle_list[idx].dds_target.turns_to_pow(self.phase_tickle_turns[idx])
        for idx in range(len(self.phase_tickle_turns))])


        freq_tickle_detuning_ftw_list = [
        self.dds_pulse_shaper_tickle_list[0].dds_target.frequency_to_ftw(freq_tickle_detuning_khz*kHz)
        for freq_tickle_detuning_khz in self.freq_tickle_detunings_mode_khz_list]

        return freq_tickle_detuning_ftw_list


    def _prepare_argument_checks(self) -> TNone:
        """
        Check experiment arguments for validity.
        """


        '''
        BICHROMATIC/CAT CHECKS
        '''

        #ensure same number of parameters are given for att and ampl
        if array_shape(self.ampls_cat_mode_pct) != array_shape(self.atts_cat_mode_db):
            raise ValueError('Must provide the same number of amplitude parameters as attenuation parameters for '
                             'the cat pulses for each mode\n'
                             'Please check the "ampls_cat_mode_pct" and "atts_cat_mode_db" parameters.')

        if (array_shape(self.freq_secular_khz_list) != array_shape(self.att_tickle_modes_db)
            or array_shape(self.att_tickle_modes_db) != array_shape(self.ampl_tickle_modes_pct)
            or array_shape(self.ampl_tickle_modes_pct)!= array_shape(self.phase_tickle_turns)):

            raise ValueError('Must provide the same number of secular frequencies, amplitudes, and attenuations for '
                             'the tickle pulses for each mode\n'
                             'Please check the "freq_secular_khz_list", "att_tickle_modes_db", '
                             '"ampl_tickle_modes_pct", and "phase_tickle_turns"parameters.')

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_experiment_list),
                3)

    '''
    MAIN SEQUENCE
    '''

    @kernel(flags={"fast-math"})
    def initialize_experiment(self) -> TNone:
        # record general subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.sidebandcool_subsequence.record_dma()
        self.readout_subsequence.record_dma()
        self.core.break_realtime()

        # initialize ramper
        self.dds_ramper_ms.sequence_initialize()
        delay_mu(50000)


    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        # # predeclare variables ahead of time
        _loop_iter = 0  # used to check_termination more frequently

        # MAIN LOOP
        for trial_num in range(self.repetitions):
            for config_vals in self.config_experiment_list:
                '''
                PREPARE & CONFIGURE
                '''
                # extract values from config list
                idx_freq_secular = int32(config_vals[0])
                freq_tickle_detuning_ftw = int32(config_vals[1])

                freq_secular_ftw = int32(self.freq_secular_ftw_list[idx_freq_secular])

                self.update_configuration(idx_freq_secular)
                '''
                BEGIN MAIN SEQUENCE
                '''
                self.core.break_realtime()
                dds_pulse_shaper_tickle = self.dds_pulse_shaper_tickle_list[idx_freq_secular]
                dds_pulse_shaper_tickle.sequence_initialize()

                self.core.break_realtime()  # add slack for execution
                delay_mu(125000)  # add even more slack lol

                # calculate for ms gate timing
                if self.enable_dynamical_decoupling:
                    time_ms_gate_dd_mu = self.time_ms_gate_mu - self.time_urukul_reset_mu - self.dds_ramper_ms.ramp_firing_delay
                    if time_ms_gate_dd_mu < 0:
                        time_ms_gate_dd_mu = 8
                else:
                    time_ms_gate_dd_mu = self.time_ms_gate_mu - (self.dds_ramper_ms.ramp_firing_delay >> 1)

                """
                Wait for Line Trigger
                """
                if self.enable_linetrigger:
                    self.trigger_line.trigger(self.trigger_line.time_timeout_mu, self.trigger_line.time_holdoff_mu)

                '''
                Relock Intensity Servo
                '''
                if self.enable_servo_relock:
                    self.qubit.relock_intensity_servo(self.time_servo_relock_mu)

                '''
                INITIALIZE ION STATE
                '''
                # initialize ion in S-1/2 state & SBC to ground state
                self.initialize_subsequence.run_dma()
                self.sidebandcool_subsequence.run_dma()

                # set cfr1 so we clear phases of all urukul0 channels on next io_update
                self.qubit.off()
                self.qubit.set_cfr1(phase_autoclear=1)
                self.qubit.singlepass0.set_cfr1(phase_autoclear=1)
                self.qubit.singlepass1.set_cfr1(phase_autoclear=1)
                self.qubit.singlepass2.set_cfr1(phase_autoclear=1)

                # set tickle frequency/phases

                dds_pulse_shaper_tickle.dds_target.set_ftw(freq_secular_ftw + freq_tickle_detuning_ftw)
                dds_pulse_shaper_tickle.dds_target.set_pow(self.phase_tickle_pow_list[idx_freq_secular])
                dds_pulse_shaper_tickle.dds_target.set_att_mu(self.att_tickle_mu_list[idx_freq_secular])
                dds_pulse_shaper_tickle.dds_target.sw.off()
                delay_mu(8)

                # set up config of shaped pulses to be fired for tickling
                # also sets up phase autoclear
                time_actual_tickle_list_mu = dds_pulse_shaper_tickle.configure_train_all_dds([self.time_tickle_mu])

                # set correct profile
                self.qubit.set_profile(self.profile_729_bichromatic)

                # set cfr1 so we clear phases of all urukul0 channels on next io_update
                self.qubit.set_cfr1(phase_autoclear=1)
                self.qubit.singlepass0.set_cfr1(phase_autoclear=1)
                self.dds_ramper_ms.set_autoclear_phase_accumulator_all_dds()

                ref_time_mu = (now_mu() + 8) & ~7
                at_mu(ref_time_mu)
                dds_pulse_shaper_tickle.dds_targets[0].cpld.io_update.pulse_mu(8)
                at_mu(ref_time_mu)
                self.qubit.io_update()
                # for ururuk channel used for tickling keep RAM enabled but ensure we don't clear phase on io_update
                dds_pulse_shaper_tickle.dds_targets[0].set_cfr1(ram_enable=1, phase_autoclear=0,
                                                                     ram_destination=ad9910.RAM_DEST_ASF)
                dds_pulse_shaper_tickle.dds_targets[0].cpld.io_update.pulse_mu(8)

                # reset cfr1 so we no longer clear phases on io_update
                self.qubit.set_cfr1()
                self.qubit.singlepass0.set_cfr1()
                self.dds_ramper_ms.reset_cfr1_all_dds()
                self.qubit.io_update()

                '''
                MS Gate
                '''
                if self.enable_ms_gate:
                    self.qubit.cpld.set_all_att_mu(self.att_reg_ms_gate)
                    self.pulse_ms(self.index_729_ms,
                                time_ms_gate_dd_mu,
                                phase_track=True,
                                ref_time_mu=ref_time_mu)

                '''
                CAT #1
                '''
                self.qubit.cpld.set_all_att_mu(self.atts_cat_reg_mu_list[idx_freq_secular])
                # cat1 - bichromatic cat pulse
                self.pulse_bichromatic(self.index_729_cat1,
                                       self.time_cat_bichromatic_mu,
                                       self.phase_cat_dynamical_decoupling_pow,
                                       phase_track=True,
                                       ref_time_mu=ref_time_mu
                                       )

                '''
                TICKLE PULSE
                '''
                # ensure we do not have a negative tickle time when DDing
                # with DD minimuim tickle time is ~1.3us, i.e time it takes to switch DD phase
                time_tickle_mu = time_actual_tickle_list_mu[0] >> 1
                if self.enable_dynamical_decoupling:
                    time_tickle_dd_mu = time_tickle_mu - self.time_urukul_reset_mu - self.dds_ramper_ms.ramp_firing_delay
                    if time_tickle_dd_mu < 0:
                        time_tickle_dd_mu = 8
                else:
                    time_tickle_dd_mu = time_tickle_mu - (self.dds_ramper_ms.ramp_firing_delay >> 1)

                if self.enable_dynamical_decoupling:
                    self.set_carrier_phase(self.phase_cat_dynamical_decoupling_pow - self.phase_dd_phase_shift_pow,
                                           phase_track=True, ref_time_mu=ref_time_mu)
                with parallel:
                    dds_pulse_shaper_tickle.run_train_all_dds()
                    if self.enable_dynamical_decoupling:
                        # see time to set profile
                        delay_mu(dds_pulse_shaper_tickle.ram_firing_delay)
                        self.qubit.singlepass0_on()
                        self.qubit.on()
                        delay_mu(time_tickle_dd_mu)
                        self.qubit.singlepass0_off()
                        self.set_carrier_phase(self.phase_cat_dynamical_decoupling_pow,
                                               phase_track=True,
                                               ref_time_mu=ref_time_mu)
                        self.qubit.singlepass0_on()
                self.qubit.off()

                '''
                CAT #2
                '''
                self.pulse_bichromatic(self.index_729_cat2,
                                       self.time_cat_bichromatic_mu,
                                       self.phase_cat_dynamical_decoupling_pow,
                                       phase_track=True,
                                       ref_time_mu=ref_time_mu)


                '''
                READ OUT & STORE RESULTS
                '''
                # read out fluorescence & clean up loop
                self.readout_subsequence.run_dma()
                counts_res = self.readout_subsequence.fetch_count()


                # cleanup dds_pulse_shaper_tickle
                dds_pulse_shaper_tickle.sequence_cleanup()

                # store results
                self.update_results(freq_secular_ftw,
                                    counts_res,
                                    freq_tickle_detuning_ftw)


                # check termination more frequently in case reps are low
                if _loop_iter % 100 == 0:
                    self.check_termination()
                _loop_iter += 1

            # rescue ion as needed & support graceful termination
            self.check_termination()

    @kernel(flags={'fast-math'})
    def pulse_bichromatic(self, index: TInt32,
                          time_pulse_mu: TInt64,
                          phase_dd_pow: TInt32 = 0,
                          phase_track: TBool = False,
                          ref_time_mu: TInt64 = 0) -> TNone:
        """
        Bichromatic interaction to produce a cat state
        :param index: urukul index with the proper settings
        :param time_pulse_mu: how long to apply the bichromatic interaction (which creates the cat state) for
        :param phase_dd_pow: phase (in pow) of dynamical decoupling pulse
        """
        # set everything to correct index
        self.qubit.off()
        self.setup_beam_profile(index,
                                phase_track,
                                ref_time_mu)
        at_mu((now_mu() + 8) & ~7)
        self.qubit.io_update()

        if self.enable_dynamical_decoupling:
            # correct time to account (and ensure it is non-negative) for call to set_mu to change dd phase
            # minimium CAT time with DD is ~1.3us, i.e time for DD phase shift
            if time_pulse_mu - self.time_urukul_reset_mu < 8:
                time_pulse_mu = 8
            else:
                time_pulse_mu = time_pulse_mu - self.time_urukul_reset_mu
            self.qubit.singlepass0_on()
        else:
            self.qubit.singlepass0_off()

        # turn on sidebands
        self.qubit.singlepass1_on()
        self.qubit.singlepass2_on()
        # turn on main doublepass and begin bichromatic
        self.qubit.on()
        delay_mu(time_pulse_mu)
        if self.enable_dynamical_decoupling:
            self.qubit.singlepass0_off()
            self.set_carrier_phase(phase_dd_pow - self.phase_dd_phase_shift_pow,
                                   phase_track=True, ref_time_mu=ref_time_mu)
            self.qubit.singlepass0_on()
        delay_mu(time_pulse_mu)

        # turn off all beams except singlepass 0 to prevent thermal fluctuations on the singlepass
        self.qubit.off()
        self.qubit.singlepass1_off()
        self.qubit.singlepass2_off()
        self.qubit.singlepass0_on()

    @kernel(flags={'fast-math'})
    def pulse_ms(self,
                 index,
                 time_ms_gate_dd_mu: TInt64,
                 ref_time_mu: TInt64,
                 phase_track=True):
        if phase_track == False:
            phase_mode = ad9910.PHASE_MODE_CONTINUOUS
        else:
            phase_mode = ad9910.PHASE_MODE_TRACKING

        # set everything to correct profile
        self.qubit.off()
        self.setup_beam_profile(index,
                                phase_track,
                                ref_time_mu)
        self.dds_ramper_ms.configure_ramp_all_dds(1)

        '''RISING PORITION OF PULSE'''
        self.qubit.on()
        self.dds_ramper_ms.run_ramp_all_dds()
        delay_mu(self.dds_ramper_ms.drg_time_ramp_mu[0])
        self.dds_ramper_ms.run_ramp_all_dds()

        '''FLAT TOP PORITION OF PULSE'''
        if self.enable_dynamical_decoupling:
            delay_mu(self.dds_ramper_ms.ramp_firing_delay)
            self.qubit.singlepass0_on()
        delay_mu(time_ms_gate_dd_mu)
        if self.enable_dynamical_decoupling:
            self.qubit.singlepass0_off()
            self.set_carrier_phase(self.phase_ms_dynamical_decoupling_pow + self.phase_dd_phase_shift_pow,
                                   phase_track=True, ref_time_mu=ref_time_mu)
            self.qubit.singlepass0_on()
        with parallel:
            # this must go after set_mu because set_mu calls io_update
            delay_mu(time_ms_gate_dd_mu)
            self.dds_ramper_ms.configure_falling_ramp_all_dds()

        '''FALLING EDGE OF PULSE'''
        self.qubit.singlepass0_off()
        self.dds_ramper_ms.run_ramp_all_dds()
        delay_mu(self.dds_ramper_ms.drg_time_ramp_mu[0])
        self.dds_ramper_ms.switch_off_all_dds()
        self.qubit.off()

        '''ENSURE ALL CFR REGISTERS ARE RESET TO NORMAL VALUES FOR SINGLE TONE OPERATION'''
        self.dds_ramper_ms.reset_cfrs_all_dds()

    @kernel(flags={"fast-math"})
    def set_carrier_phase(self, phase_dd_pow_: TInt32,
                          phase_track: TBool = False,
                          ref_time_mu: TInt64 = 0) -> TNone:
        """
        Change phase of singlepass0 for spin refocusing during DD
        :param phase_dd_pow_: phase of DD pulse used
        """
        if phase_track == False:
            phase_mode = ad9910.PHASE_MODE_CONTINUOUS
        else:
            phase_mode = ad9910.PHASE_MODE_TRACKING
        self.qubit.singlepass0.set_mu(
            self.qubit.freq_singlepass0_default_ftw,
            asf=self.ampl_dynamical_decoupling_asf,
            pow_=phase_dd_pow_,
            profile=self.profile_729_bichromatic,
            phase_mode=phase_mode,
            ref_time_mu=ref_time_mu
        )

    @kernel(flags={"fast-math"})
    def setup_beam_profile(self,
                           index: TInt32,
                           phase_track: TBool = False,
                           ref_time_mu: TInt64 = 0) -> TNone:
        """
        Configure parameters for relevant indexs on urukul
        """
        # ensure all beams are off
        self.qubit.off()
        self.qubit.singlepass0_off()
        self.qubit.singlepass1_off()
        self.qubit.singlepass2_off()
        # set up relevant beam waveforms
        if phase_track == False:
            phase_mode = ad9910.PHASE_MODE_CONTINUOUS
        else:
            phase_mode = ad9910.PHASE_MODE_TRACKING
        self.qubit.set_mu(
            self.freq_beams_ftw_list[index][0],
            asf=self.ampl_beams_asf_list[index][0],
            pow_=self.phase_beams_pow_list[index][0],
            profile=self.profile_729_bichromatic,
            phase_mode=phase_mode,
            ref_time_mu=ref_time_mu
        )
        self.qubit.singlepass0.set_mu(
            self.freq_beams_ftw_list[index][1],
            asf=self.ampl_beams_asf_list[index][1],
            pow_=self.phase_beams_pow_list[index][1],
            profile=self.profile_729_bichromatic,
            phase_mode=phase_mode,
            ref_time_mu=ref_time_mu
        )
        self.qubit.singlepass1.set_mu(
            self.freq_beams_ftw_list[index][2],
            asf=self.ampl_beams_asf_list[index][2],
            pow_=self.phase_beams_pow_list[index][2],
            profile=self.profile_729_bichromatic,
            phase_mode=phase_mode,
            ref_time_mu=ref_time_mu
        )

        self.qubit.singlepass2.set_mu(
            self.freq_beams_ftw_list[index][3],
            asf=self.ampl_beams_asf_list[index][3],
            pow_=self.phase_beams_pow_list[index][3],
            profile=self.profile_729_bichromatic,
            phase_mode=phase_mode,
            ref_time_mu=ref_time_mu
        )

    @rpc
    def set_default_configuration(self):
        """
        Store the values common across all indexs for all shots
        """
        for index in self.indices:
            self.ampl_beams_asf_list[index][0] = self.qubit.ampl_qubit_asf
            self.phase_beams_pow_list[index][0] = 0
            self.freq_beams_ftw_list[index][1] = self.qubit.freq_singlepass0_default_ftw
            self.freq_beams_ftw_list[index][0] = self.freq_carrier_ftw

        # set up values consistent with bichromatic operations
        for index in [self.index_729_cat1, self.index_729_cat2, self.index_729_ms]:
            self.ampl_beams_asf_list[index][1] = self.ampl_dynamical_decoupling_asf

        self.phase_beams_pow_list[self.index_729_cat1][2] = self.phases_pulse1_cat_pow[0]
        self.phase_beams_pow_list[self.index_729_cat1][3] = self.phases_pulse1_cat_pow[1]

        # prepare phase arrays for bichromatic
        cat4_phases = [
            self.phases_pulse1_cat_pow[0] + self.phase_cat_update_dir[0] * self.phase_cat2_offset,
            self.phases_pulse1_cat_pow[1] + self.phase_cat_update_dir[1] * self.phase_cat2_offset,
        ]

        # set up values for cat
        self.phase_beams_pow_list[self.index_729_cat1][1] = self.phase_cat_dynamical_decoupling_pow

        self.phase_beams_pow_list[self.index_729_cat2][1] = self.phase_cat_dynamical_decoupling_pow
        self.phase_beams_pow_list[self.index_729_cat2][2] = cat4_phases[0]
        self.phase_beams_pow_list[self.index_729_cat2][3] = cat4_phases[1]

        # ms gate parameters
        self.ampl_beams_asf_list[self.index_729_ms][2] = self.ampls_ms_asf[0]
        self.ampl_beams_asf_list[self.index_729_ms][3] = self.ampls_ms_asf[1]

        ms_phases = [
            self.phase_ms_pow[0],
            self.phase_ms_pow[1],
        ]

        # set up values for ms gate
        self.freq_beams_ftw_list[self.index_729_ms][2] = self.qubit.freq_singlepass1_default_ftw - self.freq_ms_mode_ftw - self.freq_ms_gate_secular_detuning_ftw
        self.freq_beams_ftw_list[self.index_729_ms][3] = self.qubit.freq_singlepass2_default_ftw + self.freq_ms_mode_ftw + self.freq_ms_gate_secular_detuning_ftw

        self.phase_beams_pow_list[self.index_729_ms][1] = self.phase_ms_dynamical_decoupling_pow
        self.phase_beams_pow_list[self.index_729_ms][2] = self.phase_ms_pow[0]
        self.phase_beams_pow_list[self.index_729_ms][3] = self.phase_ms_pow[1]

    @kernel(flags={"fast-math"})
    def update_configuration(self, secular_freq_idx):
        # set up values consistent across cat indices
        for index in [self.index_729_cat1, self.index_729_cat2]:
            self.ampl_beams_asf_list[index][2] = self.ampls_cat_asf_list[secular_freq_idx][0]
            self.ampl_beams_asf_list[index][3] = self.ampls_cat_asf_list[secular_freq_idx][1]

            self.freq_beams_ftw_list[index][2] = self.qubit.freq_singlepass1_default_ftw - self.freq_secular_ftw_list[secular_freq_idx]
            self.freq_beams_ftw_list[index][3] = self.qubit.freq_singlepass2_default_ftw + self.freq_secular_ftw_list[secular_freq_idx]


    def analyze_experiment(self):
        results, num_states = self._process_results()

        try:
            fit_funcs = self._get_fit_funcs(num_states)
        except NotImplementedError as e:
            print("Unable to find fit funcs")
            fit_funcs = None

        for secular_idx in results.keys():
            detuning_list = []
            population_list = []
            fit_x_list = []
            fit_y_list = []

            for state in range(num_states):
                secular = results[secular_idx]['secular']
                detunings =  results[secular_idx]['detuning']
                detuning_list.append(detunings)
                populations = results[secular_idx]['populations'][:, state]
                population_list.append(populations)

                fit_x = np.linspace(np.min(detunings), np.max(detunings), 1000)
                fit_x_list.append(fit_x)

                try:
                    popt, pcov = curve_fit(fit_funcs[state], detunings, populations)
                    fit_y = fit_funcs[state](fit_x, *popt)
                    fit_y_list.append(fit_y)
                except Exception as e:
                    fit_y_list.append([None] *len(fit_x))
                    print("Unable to fit functions")

            # format dictionary for applet plotting
            plotting_results = {'x': detuning_list,
                                'y': population_list,
                                'fit_x': fit_x_list,
                                'fit_y': fit_y_list,
                                'subplot_titles': f'Cat Linescan {secular:.2f}',
                                'subplot_x_labels': 'Tickle Detuning (kHz)',
                                'subplot_y_labels': 'State Population',
                                'rid': self.scheduler.rid,
                                'ylims': [[0, 1], [0, 1]],
                                'legend_labels': self._make_legend_labels(num_states)
                                }

            self.create_matplotlib_applet(plotting_results,
                                          name=f'Cat State Interferometer {secular:.2f}',
                                          group=['plotting', 'motional'],
                                          num_subplots=1)

    def _make_legend_labels(self, num_states):
        num_ions = num_states - 1
        return ["d" * (num_ions - i) + "b" * i for i in range(num_states)]


    def _process_results(self):
        # get results
        results_tmp = array(self.results)
        secular_arr = ftw_to_frequency_khz(array(results_tmp[:,0]))
        counts_arr = array(results_tmp[:, 1])
        detuning_arr = ftw_to_frequency_khz(array(results_tmp[:, 2]))

        secular_vals = np.unique(secular_arr)
        detuning_vals = np.unique(detuning_arr)

        # calculate fluorescence detection threshold
        threshold_list = np.sort(findThresholdScikit(counts_arr))
        num_states = len(threshold_list) + 1

        results_storer = {}

        for sec_idx, sec in enumerate(secular_vals):
            mask = secular_arr==sec

            detuning_subset = detuning_arr[mask]
            counts_subset = counts_arr[mask]

            count_states = np.digitize(counts_subset, threshold_list)

            # group all shots by identical detuning
            detunings, detuning_idx = np.unique(detuning_subset, return_inverse=True)
            population_vals = np.zeros((len(np.unique(detunings)), num_states))

            for det_idx in range(len(detunings)):
                det_mask = det_idx == detuning_idx
                for state in range(num_states):
                    population_vals[det_idx, state] = np.mean(count_states[det_mask] == state)

            results_storer[sec_idx] = {
                "secular": sec,
                "detuning": detunings,
                "populations": population_vals,
            }

        return results_storer, num_states

    def _get_fit_funcs(self, num_states):
        if num_states == 2:
            fit_funcs = self._get_single_ion_cat_lineshape()
        elif num_states == 3 and not self.enable_ms_gate:
            fit_funcs = self._get_unentangled_two_ion_cat_lineshape()
        elif num_states == 3 and self.enable_ms_gate:
            fit_funcs = self._get_entangled_two_ion_cat_lineshape()
        else:
            raise NotImplementedError
        return fit_funcs


    def _get_single_ion_cat_lineshape(self):
        return [self._fit_func_single_ion_d, self._fit_func_single_ion_b]

    def _fit_func_single_ion_b(self, delta, d, alpha, t,t0, delta_0, phi):
        delta_shift = 2 * np.pi * (delta - delta_0)
        x = delta_shift * t / 2
        mass_spec_phi = 2 * alpha * np.sin(x) / (x + 1e-12) * np.sin(x + delta_shift * t0 - phi)
        return (1 - 2 * d) * np.cos(mass_spec_phi) ** 2 + d

    def _fit_func_single_ion_d(self, delta, d, alpha, t,t0, delta_0, phi):
        delta_shift = 2 * np.pi * (delta - delta_0)
        x = delta_shift * t / 2
        mass_spec_phi = 2 * alpha * np.sin(x) / (x + 1e-12) * np.sin(x + delta_shift * t0 - phi)
        return (1 - 2 * d) * np.sin(mass_spec_phi) ** 2 + d

    def _get_unentangled_two_ion_cat_lineshape(self):
        return [self._fit_func_unentangled_two_ion_dd, self._fit_func_unentangled_two_ion_bd,
        self._fit_func_unentangled_two_ion_bb]

    def _fit_func_unentangled_two_ion_bb(self, delta, d, alpha, t, t0, delta_0, phi):
        delta_shift = 2 * np.pi * (delta - delta_0)
        x = delta_shift * t / 2
        mass_spec_phi = 2 * alpha * np.sin(x) / (x + 1e-12) * np.sin(x + delta_shift * t0 - phi)
        return (1 - 2 * d) * np.cos(mass_spec_phi) ** 4 + d

    def _fit_func_unentangled_two_ion_dd(self, delta, d, alpha, t, t0, delta_0, phi):
        delta_shift = 2 * np.pi * (delta - delta_0)
        x = delta_shift * t / 2
        mass_spec_phi = 2 * alpha * np.sin(x) / (x + 1e-12) * np.sin(x + delta_shift * t0 - phi)
        return (1 - 2 * d) * np.sin(mass_spec_phi) ** 4 + d

    def _fit_func_unentangled_two_ion_bd(self, delta, d, alpha, t, t0, delta_0, phi):
        delta_shift = 2 * np.pi * (delta - delta_0)
        x = delta_shift * t / 2
        mass_spec_phi = 2 * alpha * np.sin(x) / (x + 1e-12) * np.sin(x + delta_shift * t0 - phi)
        fit_bb = (1 - 2 * d) * np.cos(mass_spec_phi) ** 4 + d
        fit_dd =(1 - 2 * d) * np.sin(mass_spec_phi) ** 4 + d
        return 1 - fit_bb - fit_dd

    def _get_entangled_two_ion_cat_lineshape(self):
        return [self._fit_func_entangled_two_ion_dd, self._fit_func_entangled_two_ion_bd,
                self._fit_func_entangled_two_ion_bb]

    def _fit_func_entangled_two_ion_bb(self, delta, d, alpha, t, t0, delta_0, phi):
        delta_shift = 2 * np.pi * (delta - delta_0)
        x = delta_shift * t / 2
        mass_spec_phi = 4 * alpha * np.sin(x) / (x + 1e-12) * np.sin(x + delta_shift * t0 - phi)
        return (1 - 2 * d) * np.cos(mass_spec_phi) ** 2 + d

    def _fit_func_entangled_two_ion_dd(self, delta, d, alpha, t, t0, delta_0, phi):
        delta_shift = 2 * np.pi * (delta - delta_0)
        x = delta_shift * t / 2
        mass_spec_phi = 4 * alpha * np.sin(x) / (x + 1e-12) * np.sin(x + delta_shift * t0 - phi)
        return (1 - 2 * d) * np.sin(mass_spec_phi) ** 2 + d

    def _fit_func_entangled_two_ion_bd(self, delta, d, alpha, t, t0, delta_0, phi):
        delta_shift = 2 * np.pi * (delta - delta_0)
        x = delta_shift * t / 2
        mass_spec_phi = 4 * alpha * np.sin(x) / (x + 1e-12) * np.sin(x + delta_shift * t0 - phi)
        fit_bb = (1 - 2 * d) * np.cos(mass_spec_phi) ** 2 + d
        fit_dd = (1 - 2 * d) * np.sin(mass_spec_phi) ** 2 + d
        return 1 - fit_bb - fit_dd



