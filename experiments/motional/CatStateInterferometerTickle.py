from artiq.experiment import *
from artiq.coredevice import ad9910

from numpy import copy as np_copy
from numpy import array, int32, int64, arange, zeros, mean

from LAX_exp.language import *
from LAX_exp.system.subsequences import (
    InitializeQubit, SidebandCoolContinuousRAM, Readout, ReadoutAdaptive,
    RescueIon, QubitRAP
)

from LAX_exp.system.objects.PulseShaper import available_pulse_shapes
from LAX_exp.system.objects.dds_pulse_shape import DDSPulseShaper


class CatStateInterferometerTickle(LAXExperiment, Experiment):
    """
    Experiment: Cat State Interferometer Tickle

    Create and characterize cat states with projective state preparation.
    Uses adaptive readout to reduce timing overheads and extend available coherence times.
    """
    name = 'Cat State Inteferometer Tickle'
    kernel_invariants = {
        # subsequences & objects
        'initialize_subsequence', 'sidebandcool_subsequence', 'readout_subsequence',
        'readout_adaptive_subsequence', 'rescue_subsequence', 'rap_subsequence', 'dds_pulse_shaper',

        # hardware values - cat - default
        'ampl_doublepass_default_asf', 'time_herald_slack_mu', 'time_adapt_read_slack_mu', 'max_herald_attempts',

        # hardware values - cat - bichromatic
        'ampls_cat_asf', 'time_cat1_bichromatic_mu', 'phases_pulse1_cat_pow', 'phases_cat2_cat_update_dir',

        # hardware values - tickle
        'att_tickle_mu',

        # hardware values - readout
        'readout_config', 'ampl_729_readout_asf', 'freq_rap_center_ftw', 'freq_rap_dev_ftw', 'time_rap_mu',

        #   hardware values - dynamical decoupling
        'ampl_dynamical_decoupling_asf', 'att_dynamical_decoupling_mu',

        # configs
        'profile_729_SBC', 'profile_729_readout',
        'profile_729_cat1a', 'profile_729_cat1b','profile_729_cat2a', 'profile_729_cat2b',
        'profile_729_pi_pulse', 'profile_tickle_RAM'
        'att_reg_cat_interferometer', 'att_reg_readout_sbr', 'att_reg_readout_rap',
        'config_experiment_list',

        # extras
        'urukul_setup_time_mu', 'ampl_dynamical_decoupling_pi_asf', 'profiles',
        'phase_cat_shift_pow','time_dynamical_decoupling_pi_pulse_mu', 'att_reg_dynamical_decoupling_pi_pulse'
    }

    def build_experiment(self):

        # core arguments
        self.setattr_argument("repetitions", NumberValue(default=50, precision=0, step=1, min=1, max=100000))
        self.setattr_argument("readout_type", EnumerationValue(["None", "SBR", "RAP"], default="RAP"),
                              tooltip="None: NO 729nm pulses are applied before state-selective readout.\n"
                                      "SBR (Sideband Ratio): compares RSB and BSB amplitudes.\n"
                                      "RAP (Rapid Adiabatic Passage): Does RAP to measure fock state overlap.\n"
                                      "Note: readout pulses are NOT phase coherent with any bichromatic/sigma_x pulses.")

        # allocate relevant beam profiles
        self.profile_729_RAP = 0
        self.profile_729_SBC = 1
        self.profile_729_readout = 2
        self.profile_729_cat1a = 3
        self.profile_729_cat1b = 4
        self.profile_729_cat2a = 5
        self.profile_729_cat2b = 6
        self.profile_729_pi_pulse = 7

        # allocate profiles for dds tickle
        self.profile_tickle_RAM = 0

        # get subsequences
        self.sidebandcool_subsequence = SidebandCoolContinuousRAM(
            self, profile_729=self.profile_729_SBC, profile_854=3,
            ram_addr_start_729=0, ram_addr_start_854=0, num_samples=200
        )
        self.initialize_subsequence = InitializeQubit(self)
        self.readout_subsequence = Readout(self)
        self.readout_adaptive_subsequence = ReadoutAdaptive(self, time_bin_us=10, error_threshold=1e-2)
        self.rescue_subsequence = RescueIon(self)

        # relevant devices
        self.setattr_device('qubit')
        self.setattr_device('pump')
        self.setattr_device('repump_qubit')


        #todo: create LAX device for tickle channel
        self.setattr_device('urukul1_ch1')
        self.setattr_device('ttl1')

        # set build arguments
        self._build_arguments_default()
        self._build_arguments_dynamical_decoupling()
        self._build_arguments_cat1()
        self._build_arguments_cat2()
        self._build_arguments_tickle_waveform()
        self._build_arguments_tickle_pulseshape()
        self._build_arguments_readout()

        # # instantiate RAP here since it relies on experiment arguments
        self.rap_subsequence = QubitRAP(
            self, ram_profile=self.profile_729_RAP, ram_addr_start=202, num_samples=250,
            ampl_max_pct=self.ampl_rap_pct, pulse_shape="blackman"
        )

    def _build_arguments_default(self):
        """
        Build arguments for default beam parameters.
        """
        # defaults - beam values - doublepass (main)
        self.setattr_argument("ampl_doublepass_default_pct",
                              NumberValue(default=50., precision=3, step=5, min=0.01, max=50, scale=1., unit="%"),
                              group="default.cat",
                              tooltip="DDS amplitude for the main doublepass during the bichromatic pulse.")
        self.setattr_argument("att_doublepass_default_db",
                              NumberValue(default=8., precision=1, step=0.5, min=8., max=31.5, scale=1., unit="dB"),
                              group="default.cat",
                              tooltip="DDS attenuation for the main doublepass during the bichromatic pulse.")
        self.setattr_argument("ampls_cat_pct", PYONValue([50., 50.]), group='default.cat',
                              tooltip="DDS amplitudes for the singlepass DDSs during the bichromatic pulses.\n"
                                      "Should be a list of [rsb_ampl_pct, bsb_ampl_pct], which are applied to [singlepass1, singlepass2].")
        self.setattr_argument("atts_cat_db", PYONValue([11., 11.]), group='default.cat',
                              tooltip="DDS attenuations for the singlepass DDSs during the bichromatic pulses.\n"
                                      "Should be a list of [rsb_att_db, bsb_att_db], which are applied to [singlepass1, singlepass2].")
        self.setattr_argument("freq_cat_center_mhz_list", Scannable(
            default=[
                ExplicitScan([101.075]),
                CenterScan(101.075, 0.01, 0.0001, randomize=True),
                RangeScan(101.075, 101.1018, 50, randomize=True),
            ],
            global_min=60., global_max=400, global_step=1,
            unit="MHz", scale=1, precision=6
        ), group='default.cat',
                              tooltip="Center frequency for the bichromatic pulses.\n"
                                      "Note: this is applied via the main doublepass DDS.\n"
                                      "The singlepass DDSs center frequencies are set as their default values from the dataset manager "
                                      "(e.g. beams.freq_mhz.freq_singlepass0_mhz).")
        self.setattr_argument("freq_cat_secular_khz_list", Scannable(
            default=[
                ExplicitScan([701.8485]),
                CenterScan(702.01, 4, 0.1, randomize=True),
                RangeScan(699.0, 704.2, 50, randomize=True),
            ],
            global_min=0, global_max=10000, global_step=1,
            unit="kHz", scale=1, precision=3
        ), group='default.cat',
                              tooltip="Single-sided detuning frequency for the bichromatic pulses, applied via singlepass DDSs.\n"
                                      "The singlepass1 DDS is treated as the RSB, and will thus have its frequency DECREASED by this amount.\n"
                                      "Similarly, the singlepass2 DDS is treated as the BSB, and will thus have its frequency INCREASED by this amount.\n"
                                      "i.e. frequencies for [singlepass1, singlepass2] is set as [beams.freq_mhz.freq_singlepass1_mhz - freq_cat_secular_khz, "
                                      "beams.freq_mhz.freq_singlepass2_mhz + freq_cat_secular_khz].")

    def _build_arguments_dynamical_decoupling(self):
        # get arguments for dynamical decoupling
        self.setattr_argument('enable_dynamical_decoupling', BooleanValue(True),
                              tooltip='Indicate whether to apply a third rf tone on the 729 single pass AOM. This '
                                      'tone \n'
                                      'will produce a combination of linear combination of sigma_x and sigma_y \n'
                                      '(dependent on laser phase) which will rotate away any sigma_z errors \n',
                              group='default.dynamical_decoupling')

        self.setattr_argument('ampl_dynamical_decoupling_pct', NumberValue(default=50., max=50., min=0.01,
                                                                           step=5, precision=2, scale=1., unit='%'),
                              tooltip='Amplitude of third rf tone applied to the 729 single pass AOM. \n'
                                      'Stronger tones will further weaken any sigma_z errors BUT \n'
                                      'NOTE: Changing the strength will change AC stark shifts to the sidebands.',
                              group='default.dynamical_decoupling')

        self.setattr_argument('att_dynamical_decoupling_dB', NumberValue(default=28., max=31.5, min=7.,
                                                                         step=0.5, precision=1, scale=1., unit='dB'),
                              tooltip='Attenuation applied to the third rf tone applied to the 729 single pass AOM.',
                              group='default.dynamical_decoupling')

        self.setattr_argument('phase_dynamical_decoupling_cat_turns_list',
                              Scannable(default=[
                                  ExplicitScan([0.]),
                                  RangeScan(0, 1, 11),
                                  CenterScan(0.5, 1, 0.1)], unit='turns',
                                  global_min=0., global_max=2., global_step=0.1, precision=3, scale=1.0),
                              tooltip='Phase of the third rf tone applied to the 729 single pass AOM during the '
                                      'first catting operation. \n'
                                      'For dynamical decoupling to work the phase must be configured correctly as we \n'
                                      'need [H_{bi}, H_{dd}]=0 and this only occurs if these Hamiltonians contain \n'
                                      'the same linear combination of c1*sigma_{x}+c2*sigma_{y} which is determined '
                                      'by \n the laser phase.',
                              group='default.dynamical_decoupling')

        self.setattr_argument('att_dynamical_decoupling_pi_pulse_dB', NumberValue(default=5., max=31.5, min=5.,
                                                                                  step=0.5, precision=1, scale=1.,
                                                                                  unit='dB'),
                              tooltip='Attenuation applied to the carrier tone when \n'
                                      'performing dynamical decoupling pi pulse',
                              group='default.dynamical_decoupling')

        self.setattr_argument("time_dynamical_decoupling_pi_pulse_us",
                              NumberValue(default=2.2, precision=2, step=5, min=0.1, max=1000000, scale=1., unit="us"),
                              group='default.dynamical_decoupling',
                              tooltip="Pulse time for dynamical decoupling pi pulse.")

        self.setattr_argument("phase_dynamical_decoupling_pi_pulse_turns_list",
                             PYONValue([0.25, -0.25, 0.25, 0.25, -0.25]),
                              group='default.dynamical_decoupling',
                              tooltip='phase of dynamical decoupling pi pulse offset from phase of \n'
                                      'continuous dynamical decoupling. \n'
                                      'With respect to the dynamical decoupling axis \n'
                                      'a phase of 0 (0 turns) rotates about the x axis of the Bloch sphere \n'
                                      'and a phase of pi/2 (0.25 turns) rotates about the x axis of the Bloch sphere. \n'
                                      'A phase of pi/2 can help correct secular frequency errors but then the second '
                                      'part of the cat pulse must have rsb-bsb phase changed by pi.')


        self.setattr_argument('ampl_dynamical_decoupling_pi_pct', NumberValue(default=50., max=50., min=0.01,
                                                                              step=5, precision=2, scale=1., unit='%'),
                              tooltip='Amplitude of rf tone applied to the 729 single pass AOM for the \n'
                                      'dynamical decoupling pi pulse',
                              group='default.dynamical_decoupling')

        self.setattr_argument("phase_cat_shift_turns",
                              NumberValue(default=0.5, precision=2, step=0.05, min=-1., max=1., scale=1,
                                          unit='turns'),
                              group='default.dynamical_decoupling',
                              tooltip='Phase shift between first and second cat pulse within a block. \n'
                                      'A phase of keeps the first and second cat pulse with a block the same rsb '
                                      'and bsb phases. \n'
                                      'Shifting the phase is important as "phase_dynamical_decoupling_pi_pulse_turns" if'
                                      'set to a value other than 0. Then we have to shift the phase of rsb-bsb')

    def _build_arguments_cat1(self):
        """
        Build arguments for bichromatic/cat pulse #1.
        """
        # cat #1 config
        self.setattr_argument("enable_cat1_bichromatic", BooleanValue(default=True), group='cat1.config',
                              tooltip="Enables application of the 1st bichromatic pulse.")
        self.setattr_argument("enable_cat1_herald", BooleanValue(default=False), group='cat1.config',
                              tooltip="Enables spin-state heralding via state-selective fluorescence. "
                                      "Heralding only progresses if the state is dark, since otherwise, the motional state is destroyed.\n"
                                      "Pulses are applied as [bichromatic, herald, quench].\n"
                                      "Note: uses adaptive readout - ensure adaptive readout arguments are correctly set in the dataset manager.")
        self.setattr_argument("enable_cat1_quench", BooleanValue(default=False), group='cat1.config',
                              tooltip="Enables quenching via 854nm to return the spin-state to the S-1/2 state.\n"
                                      "Note: if quench is applied to a superposition state, then the result is a mixed state, not a pure state.\n"
                                      "Pulses are applied as [bichromatic, herald, quench].")
        self.setattr_argument("time_cat1_bichromatic_us",
                              NumberValue(default=50, precision=2, step=5, min=0.1, max=10000000, scale=1., unit="us"),
                              group="cat1.config",
                              tooltip="Pulse time for the 1st bichromatic pulse.")
        self.setattr_argument("phases_pulse1_cat_turns", PYONValue([0., 0.]), group='cat1.config',
                              tooltip="Relative phases for the singlepass DDSs during the 1st bichromatic pulse.\n"
                                      "Should be a list of [rsb_phase_turns, bsb_phase_turns].\n"
                                      "Note: these phases are applied to the singlepass DDSs, so do not need to be halved, "
                                      "unlike phases applied to the main doublepass.")

        # ramsey delay between cat1 & cat2
        self.setattr_argument("enable_ramsey_delay", BooleanValue(default=False), group='cat1.ramsey',
                              tooltip="Enables a Ramsey delay between the 1st and 2nd bichromatic pulses. "
                                      "Useful for doing motional coherence tests.")
        self.setattr_argument("time_ramsey_delay_us_list", Scannable(
            default=[
                ExplicitScan([100]),
                RangeScan(0, 500, 50, randomize=True),
            ],
            global_min=1, global_max=100000, global_step=1,
            unit="us", scale=1, precision=5
        ),
                              group="cat1.ramsey",
                              tooltip="Ramsey delay time between 1st and 2nd bichromatic pulses.")

    def _build_arguments_cat2(self):
        """
        Build arguments for bichromatic/cat pulse #2.
        """
        # cat2 - config
        self.setattr_argument("enable_cat2_bichromatic", BooleanValue(default=True),
                              group='cat2.bichromatic',
                              tooltip="Enables application of the 2nd bichromatic pulse.")
        self.setattr_argument("enable_cat2_herald", BooleanValue(default=False), group='cat2.bichromatic',
                              tooltip="Enables spin-state heralding via state-selective fluorescence. "
                                      "Heralding only progresses if the state is dark, since otherwise, the motional state is destroyed.\n"
                                      "Pulses are applied as [bichromatic, herald, quench].\n"
                                      "Note: uses adaptive readout - ensure adaptive readout arguments are correctly set in the dataset manager.")
        self.setattr_argument("enable_cat2_quench", BooleanValue(default=False), group='cat2.bichromatic',
                              tooltip="Enables quenching via 854nm to return the spin-state to the S-1/2 state.\n"
                                      "Note: if quench is applied to a superposition state, then the result is a mixed state, not a pure state.\n"
                                      "Pulses are applied as [bichromatic, herald, quench].")

        # cat2 - pulse parameters
        self.setattr_argument("time_cat2_cat_us_list", Scannable(
            default=[
                ExplicitScan([50.]),
                RangeScan(0, 500, 50, randomize=True),
            ],
            global_min=1, global_max=10000, global_step=1,
            unit="us", scale=1, precision=5
        ),
                              group="cat2.bichromatic",
                              tooltip="Pulse time for the 2nd bichromatic pulse.")

        self.setattr_argument("target_cat2_cat_phase",
                              EnumerationValue(['RSB', 'BSB', 'RSB-BSB', 'RSB+BSB'], default='RSB-BSB'),
                              group="cat2.bichromatic",
                              tooltip="Phase update array for the singlepass DDSs during the 2nd bichromatic pulse.\n"
                                      "This configures how phase_cat2_turns_list are to be applied to the DDSs.")
        self.setattr_argument("phase_cat2_turns_list", Scannable(
            default=[
                ExplicitScan([0.]),
                RangeScan(0, 1.0, 11, randomize=True),
            ],
            global_min=-1.0, global_max=1.0, global_step=0.1,
            unit="turns", scale=1, precision=3
        ), group="cat2.bichromatic",
                              tooltip="Phase sweep values applied to the singlepass DDSs during the 2nd bichromatic pulse.\n"
                                      "These values are multiplied/scaled by the array specified by target_cat2_cat_phase.")

    def _build_arguments_readout(self):
        """
        Build arguments for readout pulse.
        """
        # sideband-type readout (SBR)
        self.setattr_argument("ampl_729_readout_pct",
                              NumberValue(default=50., precision=3, step=5, min=0.01, max=50, unit="%", scale=1.),
                              group="read.SBR",
                              tooltip="729nm DDS amplitude (in percent of full scale) to use for readout.\n"
                                      "This is applied via the main doublepass.")
        self.setattr_argument("att_729_readout_db",
                              NumberValue(default=8., precision=1, step=0.5, min=8., max=31.5, unit="dB", scale=1.),
                              group="read.SBR",
                              tooltip="729nm DDS attenuation (in dB) to use for readout. "
                                      "This is applied via the main doublepass.")
        self.setattr_argument("freq_729_readout_mhz_list", Scannable(
            default=[
                ExplicitScan([101.4308]),
                CenterScan(101.4308, 0.01, 0.0002, randomize=True),
                RangeScan(101.4288, 101.4328, 50, randomize=True),
            ],
            global_min=60., global_max=400, global_step=1,
            unit="MHz", scale=1, precision=6
        ),
                              group="read.SBR",
                              tooltip="729nm DDS frequencies to use for readout. "
                                      "This is applied via the main doublepass.")
        self.setattr_argument("time_729_readout_us_list", Scannable(
            default=[
                ExplicitScan([27.35]),
                RangeScan(0.01, 800, 200, randomize=True),
            ],
            global_min=0.01, global_max=100000, global_step=1,
            unit="us", scale=1, precision=5
        ),
                              group="read.SBR",
                              tooltip="729nm readout pulse times.")

        # RAP-based readout
        self.setattr_argument("att_rap_db",
                              NumberValue(default=8, precision=1, step=0.5, min=8, max=31.5, unit="dB", scale=1.),
                              group="read.RAP")
        self.setattr_argument("ampl_rap_pct",
                              NumberValue(default=50., precision=3, step=5, min=1, max=50, unit="%", scale=1.),
                              group="read.RAP")
        self.setattr_argument("freq_rap_center_mhz",
                              NumberValue(default=100.755, precision=6, step=1e-2, min=60, max=200, unit="MHz",
                                          scale=1.),
                              group='read.RAP')
        self.setattr_argument("freq_rap_dev_khz",
                              NumberValue(default=72., precision=2, step=0.01, min=1, max=1e4, unit="kHz", scale=1.),
                              group='read.RAP')
        self.setattr_argument("time_rap_us",
                              NumberValue(default=400., precision=3, min=1, max=1e7, step=1, unit="us", scale=1.),
                              group="read.RAP")

    def _build_arguments_tickle_waveform(self):
        """
        Build core sweep arguments for the tickle pulse.
        """
        _argstr = "tickle"  # string to use for arguments

        # waveform - parameter sweeps
        self.setattr_argument("enable_tickle_pulse", BooleanValue(default=True),
                              group='{}.waveform'.format(_argstr),
                              tooltip="Enables the tickle pulse.")
        self.setattr_argument("freq_tickle_mhz_list", Scannable(
            default=[
                ExplicitScan([1.2]),
                CenterScan(1.2, 0.1, 0.001, randomize=True),
                RangeScan(1.2, 1.3, 26, randomize=True),
            ],
            global_min=1, global_max=400, global_step=0.001,
            unit="MHz", scale=1, precision=6),
                              group="{}.waveform".format(_argstr),
                              tooltip="Frequency of tickle pulse (in MHz) applied via the urukul dds.")

        self.setattr_argument("phase_tickle_turns_list", Scannable(
            default=[
                ExplicitScan([0.]),
                RangeScan(0, 1.0, 26, randomize=True),
            ],
            global_min=0.0, global_max=1.0, global_step=1,
            unit="turns", scale=1, precision=5),
                              group="{}.waveform".format(_argstr),
                              tooltip="Phase of tickle pulse (in turns) applied via the urukul dds.")

        self.setattr_argument("ampl_tickle_pct",
                              NumberValue(default=50., precision=2, min=0., max=50., unit="%", scale=1.),
                              group="{}.waveform".format(_argstr),
                              tooltip='Amplitude of tickle pulse.')

        self.setattr_argument("time_heating_us",
                              NumberValue(default=50, precision=2, step=500, min=0.04, max=10000000, unit="us",
                                          scale=1.),
                              group="{}.waveform".format(_argstr),
                              tooltip="Time for the total pulse (including pulse shape).")
        self.setattr_argument("att_tickle_db",
                              NumberValue(default=25., precision=1, step=0.5, min=8., max=31.5, unit="dB", scale=1.),
                              group="{}.waveform".format(_argstr),
                              tooltip="Attenuation to be used for the urukul channel used for generating the tickle.")

    def _build_arguments_tickle_pulseshape(self):
        """
        Build core modulation arguments for the tickle pulse.
        """
        _argstr = "tickle"  # string to use for arguments

        # waveform - pulse shaping
        self.setattr_argument("enable_pulse_shaping", BooleanValue(default=False),
                              group='{}.shape'.format(_argstr),
                              tooltip="Applies pulse shaping to the edges of the tickle pulse.")
        self.setattr_argument("type_pulse_shape",
                              EnumerationValue(list(available_pulse_shapes.keys()), default='sine_squared'),
                              group='{}.shape'.format(_argstr),
                              tooltip="Pulse shape type to be used.")


    def prepare_experiment(self):
        """
        Prepare & precompute experimental values.
        """
        self._prepare_argument_checks()

        ### Build Pulse Shaper
        self.dds_pulse_shaper = DDSPulseShaper(dds_target= self.urukul1_ch1,
                                              ram_profile=self.profile_tickle_RAM,
                                              ram_addr_start=202, num_samples=250,
                                              ampl_max_pct=self.ampl_tickle_pct)

        ### MAGIC NUMBERS ###
        self.time_adapt_read_slack_mu = self.core.seconds_to_mu(20 * us)  # post-heralding slack to fix RTIOUnderflows
        self.time_herald_slack_mu = self.core.seconds_to_mu(150 * us)  # add slack only if herald success
        self.max_herald_attempts = 200  # max herald attempts before config skipped
        self.urukul_setup_time_mu = int64(8)  # extra delay between calls to ad9910.set_mu

        # run component preparation
        freq_729_readout_ftw_list, time_729_readout_mu_list = self._prepare_experiment_readout()
        phase_dynamical_decoupling_cat_pow_list = self._prepare_experiment_dynamical_decoupling()
        (freq_cat_center_ftw_list, freq_cat_secular_ftw_list, time_cat2_cat_mu_list,
         phase_cat2_cat_pow_list, time_ramsey_delay_mu_list) = self._prepare_experiment_cat_general()
        freq_tickle_ftw_list, phase_tickle_pow_list = self._prepare_experiment_tickle()


        # create experiment config
        self.config_experiment_list = create_experiment_config(
            # cat sweeps
            freq_cat_center_ftw_list, freq_cat_secular_ftw_list,
            time_cat2_cat_mu_list, phase_cat2_cat_pow_list,
            freq_729_readout_ftw_list, time_729_readout_mu_list,
            time_ramsey_delay_mu_list,

            # tickle sweeps
            freq_tickle_ftw_list, phase_tickle_pow_list,

            # dynamical decoupling sweeps
            phase_dynamical_decoupling_cat_pow_list,
            config_type=float, shuffle_config=True
        )

        self.phase_dynamical_decoupling_pi_pulse_pow_list = [self.qubit.singlepass0.turns_to_pow(
            phase_dynamical_decoupling_pi_pulse_turns) for phase_dynamical_decoupling_pi_pulse_turns in
                                                            self.phase_dynamical_decoupling_pi_pulse_turns_list]
        self.ampl_dynamical_decoupling_pi_asf = self.qubit.singlepass0.amplitude_to_asf(
            self.ampl_dynamical_decoupling_pi_pct / 100.)

        # create profile list for later use
        self.profiles = [self.profile_729_cat1a, self.profile_729_cat1b, self.profile_729_cat2a, self.profile_729_cat2b,
                         self.profile_729_pi_pulse, self.profile_729_readout]

        # create placeholder arrays for later uses
        self.phase_beams_pow_list = zeros((8, 4), dtype=int32)
        self.freq_beams_ftw_list = zeros((8, 4), dtype=int32)
        self.ampl_beams_asf_list = zeros((8, 4), dtype=int32)
        # schedule when to turn on tickle during experiment

        self.set_default_profile_configuration()

    def _prepare_experiment_dynamical_decoupling(self):
        if self.enable_dynamical_decoupling:
            self.ampl_dynamical_decoupling_asf = self.qubit.singlepass0.amplitude_to_asf(
                self.ampl_dynamical_decoupling_pct / 100.)
            self.att_dynamical_decoupling_mu = self.qubit.singlepass0.cpld.att_to_mu(
                self.att_dynamical_decoupling_dB * dB)
            phase_dynamical_decoupling_cat_pow_list = array(
                [self.qubit.singlepass0.turns_to_pow(phase_dynamical_decoupling_turns) for
                 phase_dynamical_decoupling_turns in self.phase_dynamical_decoupling_cat_turns_list])
            self.phase_cat_shift_pow = self.qubit.singlepass0.turns_to_pow(
                self.phase_cat_shift_turns)
            self.num_dynamical_decoupling_pi_pulses = len(list(self.phase_dynamical_decoupling_pi_pulse_turns_list)) - 2
            print(self.num_dynamical_decoupling_pi_pulses)
        else:
            self.phase_cat_shift_pow = 0
            self.att_dynamical_decoupling_mu = self.qubit.att_singlepass0_default_mu
            phase_dynamical_decoupling_cat_pow_list = array([0])
            self.ampl_dynamical_decoupling_asf = self.qubit.ampl_singlepass0_default_asf
            self.num_dynamical_decoupling_pi_pulses = 0

        self.time_dynamical_decoupling_pi_pulse_mu = (
            self.core.seconds_to_mu(self.time_dynamical_decoupling_pi_pulse_us * us))

        return phase_dynamical_decoupling_cat_pow_list

    def _prepare_experiment_readout(self):
        """
        Prepare experiment values for state readout.
        :return: tuple of (freq_729_readout_ftw_list, time_729_readout_mu_list)
        """
        # configure readout method
        if 'SBR' in self.readout_type:  # 1 is SBR
            self.readout_config = 1
            freq_729_readout_ftw_list = array([self.qubit.frequency_to_ftw(freq_mhz * MHz)
                                               for freq_mhz in self.freq_729_readout_mhz_list], dtype=int32)
            time_729_readout_mu_list = array([self.core.seconds_to_mu(time_us * us)
                                              for time_us in self.time_729_readout_us_list], dtype=int64)

        elif 'RAP' in self.readout_type:  # 2 is RAP
            self.readout_config = 2
            # return dummy readout sweep arrays if using RAP
            freq_729_readout_ftw_list = array([1], dtype=int32)
            time_729_readout_mu_list = array([256], dtype=int64)

        elif 'None' in self.readout_type:  # -1 is None
            self.readout_config = -1
            # return dummy readout sweep arrays if no 729nm-based readout
            freq_729_readout_ftw_list = array([1], dtype=int32)
            time_729_readout_mu_list = array([256], dtype=int64)

        # do a final check of the readout_type argument
        if not any(kw in self.readout_type for kw in ('None', 'RAP', 'SBR')):
            raise ValueError("Invalid readout type. Must be one of (None, SBR, RAP).")

        # prepare RAP arguments
        self.freq_rap_center_ftw = self.qubit.frequency_to_ftw(self.freq_rap_center_mhz * MHz)
        self.freq_rap_dev_ftw = self.qubit.frequency_to_ftw(self.freq_rap_dev_khz * kHz)
        self.time_rap_mu = self.core.seconds_to_mu(self.time_rap_us * us)

        # readout ampl for SBR
        self.ampl_729_readout_asf = self.qubit.amplitude_to_asf(self.ampl_729_readout_pct / 100.)

        # return sweep arrays for building the expconfig
        return freq_729_readout_ftw_list, time_729_readout_mu_list

    def _prepare_experiment_cat_general(self):
        """
        Prepare experiment values for cat/bichromatic.
        :return: tuple of (freq_cat_center_ftw_list, freq_cat_secular_ftw_list, time_cat2_cat_mu_list,
            phase_cat2_cat_pow_list, time_ramsey_delay_mu_list)
        """
        '''
        CONVERT VALUES TO MACHINE UNITS - BICHROMATIC/CAT DEFAULTS
        '''
        # defaults - main doublepass (near chamber)
        self.ampl_doublepass_default_asf = self.qubit.amplitude_to_asf(self.ampl_doublepass_default_pct / 100.)

        # defaults - cat
        freq_cat_center_ftw_list = array([self.qubit.frequency_to_ftw(freq_mhz * MHz)
                                          for freq_mhz in self.freq_cat_center_mhz_list])
        freq_cat_secular_ftw_list = array([self.qubit.singlepass0.frequency_to_ftw(freq_khz * kHz)
                                           for freq_khz in self.freq_cat_secular_khz_list])
        self.ampls_cat_asf = array([self.qubit.singlepass0.amplitude_to_asf(ampl_pct / 100.)
                                    for ampl_pct in self.ampls_cat_pct])

        '''
        CONVERT VALUES TO MACHINE UNITS - BICHROMATIC/CAT PULSES
        '''
        # cat1 values
        self.time_cat1_bichromatic_mu = self.core.seconds_to_mu(self.time_cat1_bichromatic_us / 2 * us)
        self.phases_pulse1_cat_pow = [self.qubit.singlepass0.turns_to_pow(phas_pow)
                                      for phas_pow in self.phases_pulse1_cat_turns]

        # inter-cat ramsey delay
        if self.enable_ramsey_delay:
            time_ramsey_delay_mu_list = [self.core.seconds_to_mu(time_delay_us / (self.num_dynamical_decoupling_pi_pulses+1) * us)
                                         for time_delay_us in self.time_ramsey_delay_us_list]
        else:
            time_ramsey_delay_mu_list = [0]

        # cat2 values
        if self.enable_cat2_bichromatic:
            time_cat2_cat_mu_list = [self.core.seconds_to_mu(time_us / 2 * us)
                                     for time_us in self.time_cat2_cat_us_list]
            phase_cat2_cat_pow_list = [self.qubit.singlepass0.turns_to_pow(phas_pow)
                                       for phas_pow in self.phase_cat2_turns_list]
        else:
            time_cat2_cat_mu_list = [0]
            phase_cat2_cat_pow_list = [0]

        # specify phase update array based on user arguments
        if self.target_cat2_cat_phase == 'RSB':
            self.phases_cat2_cat_update_dir = array([1, 0], dtype=int32)
        elif self.target_cat2_cat_phase == 'BSB':
            self.phases_cat2_cat_update_dir = array([0, 1], dtype=int32)
        elif self.target_cat2_cat_phase == 'RSB-BSB':
            self.phases_cat2_cat_update_dir = array([1, -1], dtype=int32)
        elif self.target_cat2_cat_phase == 'RSB+BSB':
            self.phases_cat2_cat_update_dir = array([1, 1], dtype=int32)

        '''
        CREATE ATTENUATION REGISTERS
        '''
        # attenuation register - bichromatic: main doublepass set to specified experiment argument value
        self.att_reg_cat_interferometer = 0x00000000 | (
                (att_to_mu(self.att_doublepass_default_db * dB) << ((self.qubit.beam.chip_select - 4) * 8)) |
                (self.att_dynamical_decoupling_mu << ((self.qubit.singlepass0.chip_select - 4) * 8)) |
                (att_to_mu(self.atts_cat_db[0] * dB) << ((self.qubit.singlepass1.chip_select - 4) * 8)) |
                (att_to_mu(self.atts_cat_db[1] * dB) << ((self.qubit.singlepass2.chip_select - 4) * 8))
        )

        self.att_reg_dynamical_decoupling_pi_pulse = 0x00000000 | (
                (att_to_mu(8 * dB) << ((self.qubit.beam.chip_select - 4) * 8)) |
                (att_to_mu(self.att_dynamical_decoupling_pi_pulse_dB * dB) << ((self.qubit.singlepass0.chip_select - 4) * 8)) |
                (att_to_mu(31.5 * dB) << ((self.qubit.singlepass1.chip_select - 4) * 8)) |
                (att_to_mu(31.5 * dB) << ((self.qubit.singlepass2.chip_select - 4) * 8))
        )

        # attenuation register - readout (SBR): singlepasses set to default
        self.att_reg_readout_sbr = 0x00000000 | (
                (att_to_mu(self.att_729_readout_db * dB) << ((self.qubit.beam.chip_select - 4) * 8)) |
                (self.qubit.att_singlepass0_default_mu << ((self.qubit.singlepass0.chip_select - 4) * 8)) |
                (self.qubit.att_singlepass1_default_mu << ((self.qubit.singlepass1.chip_select - 4) * 8)) |
                (self.qubit.att_singlepass2_default_mu << ((self.qubit.singlepass2.chip_select - 4) * 8))
        )

        # attenuation register - readout (RAP): singlepasses set to default
        self.att_reg_readout_rap = 0x00000000 | (
                (att_to_mu(self.att_rap_db * dB) << ((self.qubit.beam.chip_select - 4) * 8)) |
                (self.qubit.att_singlepass0_default_mu << ((self.qubit.singlepass0.chip_select - 4) * 8)) |
                (self.qubit.att_singlepass1_default_mu << ((self.qubit.singlepass1.chip_select - 4) * 8)) |
                (self.qubit.att_singlepass2_default_mu << ((self.qubit.singlepass2.chip_select - 4) * 8))
        )

        # return sweep lists required for
        return (freq_cat_center_ftw_list, freq_cat_secular_ftw_list, time_cat2_cat_mu_list,
                phase_cat2_cat_pow_list, time_ramsey_delay_mu_list)

    def _prepare_experiment_tickle(self):
        """
        Prepare general experiment values for the tickle pulse.
        :return: tuple of (freq_tickle_hz_list, phase_tickle_list)
        """
        # convert values to convenience units
        self.att_tickle_mu = att_to_mu(self.att_tickle_db * dB)
        freq_tickle_ftw_list = [mhz_to_ftw(freq_tickle_mhz) for freq_tickle_mhz in self.freq_tickle_mhz_list]
        phase_tickle_pow_list = [self.urukul1_ch1.turns_to_pow(phase_tickle_turns) for phase_tickle_turns in self.phase_tickle_turns_list]
        self.time_tickle_mu = us_to_mu(self.time_heating_us / (self.num_dynamical_decoupling_pi_pulses+1))


        # don't apply sweep if tickle is disabled
        if self.enable_tickle_pulse:
            return freq_tickle_ftw_list, phase_tickle_pow_list
        else:
            return [-1], [0]

    def _prepare_argument_checks(self) -> TNone:
        """
        Check experiment arguments for validity.
        """
        # todo
        '''
        BICHROMATIC/CAT CHECKS
        '''
        pass

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_experiment_list),
                10)

    '''
    MAIN SEQUENCE
    '''

    @kernel(flags={"fast-math"})
    def initialize_experiment(self) -> TNone:
        # note: no need to remove phase_autoclear from CFR1 b/c PHASE_MODE_TRACKING does it for us
        # record general subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.sidebandcool_subsequence.record_dma()
        self.readout_subsequence.record_dma()
        self.core.break_realtime()

        # configure RAP pulse
        if self.readout_config == 2:
            self.rap_subsequence.configure(self.time_rap_mu, self.freq_rap_center_ftw, self.freq_rap_dev_ftw)
            delay_mu(50000)

        # configure tickle
        self.dds_pulse_shaper.sequence_initialize()
        self.urukul1_ch1.set_att_mu(self.att_tickle_mu)

    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        # predeclare variables ahead of time
        ion_state = (-1, 0, int64(0))  # store ion state for adaptive readout
        herald_counter = 0  # store herald attempts
        _loop_iter = 0  # used to check_termination more frequently

        # MAIN LOOP
        for trial_num in range(self.repetitions):
            for config_vals in self.config_experiment_list:
                '''
                PREPARE & CONFIGURE
                '''
                # extract values from config list
                freq_cat_center_ftw = int32(config_vals[0])
                freq_cat_secular_ftw = int32(config_vals[1])
                time_cat2_cat_mu = int64(config_vals[2])
                phase_cat2_cat_pow = int32(config_vals[3])
                freq_729_readout_ftw = int32(config_vals[4])
                time_729_readout_mu = int64(config_vals[5])
                time_ramsey_delay_mu = int64(config_vals[6])
                freq_tickle_ftw = int32(config_vals[7])
                phase_tickle_pow = int32(config_vals[8])
                phase_dynamical_decoupling_cat_pow = int32(config_vals[9])

                herald_counter = 0  # clear herald counter

                self.update_profile_configuration(freq_cat_center_ftw, freq_cat_secular_ftw,
                                     phase_cat2_cat_pow,
                                     phase_dynamical_decoupling_cat_pow)

                # configure tickle
                if self.enable_tickle_pulse:
                    # set urukul frequency/phases
                    self.core.break_realtime()
                    self.urukul1_ch2.set_mu(freq_tickle_ftw,
                                            phase_tickle_pow,
                                            profile=self.profile_tickle_RAM)
                    self.dds_pulse_shaper.configure_train(self.time_tickle_mu)

                '''
                BEGIN MAIN SEQUENCE
                '''
                while True:
                    # check heralding OK (otherwise execution is blocked)
                    if herald_counter >= self.max_herald_attempts:
                        print("\t\tWarning: too many heralds. Moving onto next configuration.")
                        self.check_termination()  # check termination b/c we haven't in a while
                        self.core.break_realtime()  # add slack (not really necessary but whatevs)
                        break

                    self.core.break_realtime()  # add slack for execution
                    delay_mu(125000)  # add even more slack lol

                    '''
                    INITIALIZE ION STATE
                    '''
                    # initialize ion in S-1/2 state & SBC to ground state
                    self.initialize_subsequence.run_dma()
                    self.sidebandcool_subsequence.run_dma()

                    # set cfr1 so we clear phases of all urukul0 channels on next io_update
                    self.qubit.off()
                    self.setup_beam_profiles()
                    self.qubit.set_cfr1(phase_autoclear=1)
                    self.qubit.singlepass0.set_cfr1(phase_autoclear=1)
                    self.qubit.singlepass1.set_cfr1(phase_autoclear=1)
                    self.qubit.singlepass2.set_cfr1(phase_autoclear=1)

                    # synchronize start time to urukul


                    # io_update clears phase of all channels to be cleared
                    self.qubit.io_update()

                    # reset cfr1 so we no longer clear phases on any io_update
                    self.qubit.set_cfr1()
                    self.qubit.singlepass0.set_cfr1()
                    self.qubit.singlepass1.set_cfr1()
                    self.qubit.singlepass2.set_cfr1()
                    self.qubit.io_update()

                    '''
                    CAT #1
                    '''
                    # cat1 - bichromatic cat pulse
                    # note: enable_dynamical_decoupling logic can be moved into pulse_cat,
                    #   since it already does a enable_dynamical_decoupling check
                    with parallel:
                        with sequential:
                            if self.enable_cat1_bichromatic:
                                self.pulse_cat(self.profile_729_cat1a, self.time_cat1_bichromatic_mu)
                                if self.enable_dynamical_decoupling:
                                    self.perform_dynamical_decoupling_pi_pulse()
                                self.pulse_cat(self.profile_729_cat1b, self.time_cat1_bichromatic_mu)

                    # cat1 - force herald (to projectively disentangle spin/motion)
                    if self.enable_cat1_herald:
                        ion_state = self.readout_adaptive_subsequence.run()
                        delay_mu(self.time_adapt_read_slack_mu)  # add slack following completion
                        self.pump.off()

                        # ensure dark state (flag is 0)
                        if ion_state[0] != 0:
                            herald_counter += 1  # increment herald counter to check for errors
                            continue
                        # otherwise, add minor slack TO RTIO COUNTER (not timeline) and proceed
                        # todo: should we be doing now_mu instead? get_rtio_counter isn't very deterministic ...
                        at_mu(self.core.get_rtio_counter_mu() + self.time_herald_slack_mu)

                    # cat1 - quench spin-up to spin-down; can be used to create mixed state
                    if self.enable_cat1_quench:
                        self.pump.readout()
                        self.repump_qubit.on()
                        delay_mu(self.initialize_subsequence.time_repump_qubit_mu)
                        self.repump_qubit.off()

                    '''
                    RAMSEY DELAY
                    '''
                    if self.enable_ramsey_delay:
                        for idx in range(self.num_dynamical_decoupling_pi_pulses):
                            with parallel:
                                delay_mu(time_ramsey_delay_mu)
                                if self.enable_dynamical_decoupling:
                                    self.qubit.on()
                                    self.qubit.singlepass0_on()
                                    self.write_pi_pulse_phase(self.phase_dynamical_decoupling_pi_pulse_pow_list[idx + 1]
                                                              + phase_dynamical_decoupling_cat_pow)
                            self.qubit.off()
                            if self.enable_dynamical_decoupling:
                                self.perform_dynamical_decoupling_pi_pulse()
                                # reset attenuators for continuous DD
                                self.qubit.cpld.set_all_att_mu(self.att_reg_cat_interferometer)
                                # reset profile for continuous DD
                                self.qubit.set_profile(self.profile_729_cat1b)
                                # ensure it aligns to a future clock cycle
                                at_mu((now_mu() + 8) & ~7)
                                self.qubit.io_update()

                        with parallel:
                            delay_mu(time_ramsey_delay_mu)
                            if self.enable_dynamical_decoupling:
                                self.qubit.on()
                                self.qubit.singlepass0_on()
                                self.write_pi_pulse_phase(self.phase_dynamical_decoupling_pi_pulse_pow_list[-1] +
                                                          phase_dynamical_decoupling_cat_pow)
                        self.qubit.off()

                    '''
                    TICKLE PULSE
                    '''
                    if self.enable_tickle_pulse:
                        for idx in range(self.num_dynamical_decoupling_pi_pulses):
                            if self.enable_dynamical_decoupling:
                                self.qubit.on()
                            self.dds_pulse_shaper.run_train_single()
                            with parallel:
                                self.qubit.off()
                                self.write_pi_pulse_phase(self.phase_dynamical_decoupling_pi_pulse_pow_list[idx+1] +
                                                          phase_dynamical_decoupling_cat_pow)
                            if self.enable_dynamical_decoupling:
                                self.perform_dynamical_decoupling_pi_pulse()
                                # reset attenuators for continuous DD
                                self.qubit.cpld.set_all_att_mu(self.att_reg_cat_interferometer)
                                # reset profile for continuous DD
                                self.qubit.set_profile(self.profile_729_cat1b)
                                # align to a future clock cycle
                                at_mu((now_mu() + 8) & ~7)
                                self.qubit.io_update()

                        if self.enable_dynamical_decoupling:
                            self.qubit.on()
                            self.qubit.singlepass0_on()
                            self.write_pi_pulse_phase(self.phase_dynamical_decoupling_pi_pulse_pow_list[-1] +
                                                      phase_dynamical_decoupling_cat_pow)
                        self.dds_pulse_shaper.run_train_single()
                        self.qubit.off()

                    '''
                    CAT #2
                    '''
                    with parallel:
                        # turn off urukul
                        if self.enable_cat2_bichromatic:
                            with sequential:
                                self.pulse_cat(self.profile_729_cat2a, time_cat2_cat_mu)
                                if self.enable_dynamical_decoupling:
                                    self.perform_dynamical_decoupling_pi_pulse()
                                self.pulse_cat(self.profile_729_cat2b, time_cat2_cat_mu)

                    # cat2 - force herald (to projectively disentangle spin/motion)
                    # todo: should we be doing now_mu instead? get_rtio_counter isn't very deterministic ...
                    if self.enable_cat2_herald:
                        ion_state = self.readout_adaptive_subsequence.run()
                        delay_mu(self.time_adapt_read_slack_mu)  # add slack following completion
                        self.pump.off()

                        # ensure dark state (flag is 0)
                        if ion_state[0] != 0:
                            herald_counter += 1  # increment herald counter to check for errors
                            continue
                        # otherwise, add minor slack TO RTIO COUNTER (not timeline) and proceed
                        at_mu(self.core.get_rtio_counter_mu() + self.time_herald_slack_mu)

                    # cat2 - quench spin-up to spin-down; can be used to create mixed state
                    if self.enable_cat2_quench:
                        self.pump.readout()
                        self.repump_qubit.on()
                        delay_mu(self.initialize_subsequence.time_repump_qubit_mu)
                        self.repump_qubit.off()

                    # force break loop by default
                    break

                '''
                READ OUT & STORE RESULTS
                '''
                # only read out if no boooboo
                if herald_counter < self.max_herald_attempts:
                    # 729nm based readout (sideband ratio, rabi flopping)
                    if self.readout_config == 1:
                        self.pulse_readout_sbr(time_729_readout_mu, freq_729_readout_ftw)
                    # RAP based readout
                    elif self.readout_config == 2:
                        self.pulse_readout_rap()

                    # read out fluorescence & clean up loop
                    self.readout_subsequence.run_dma()
                    self.rescue_subsequence.resuscitate()
                    counts_res = self.readout_subsequence.fetch_count()
                else:
                    # return -1 so user knows booboo happened
                    counts_res = -1

                # cleanup dds_pulse_shaper
                # todo: do we need this???
                self.dds_pulse_shaper.sequence_cleanup()

                # store results
                self.update_results(freq_cat_center_ftw,
                                    counts_res,
                                    freq_cat_secular_ftw,
                                    time_cat2_cat_mu << 1,
                                    phase_cat2_cat_pow,
                                    freq_729_readout_ftw,
                                    time_729_readout_mu,
                                    time_ramsey_delay_mu * (self.num_dynamical_decoupling_pi_pulses+1),
                                    freq_sweep_hz,
                                    phase_dynamical_decoupling_cat_pow)


                # check termination more frequently in case reps are low
                if _loop_iter % 100 == 0:
                    self.check_termination()
                _loop_iter += 1

            # rescue ion as needed & support graceful termination
            self.check_termination()

    '''
    HELPER FUNCTIONS
    '''
    @kernel(flags={'fast-math'})
    def perform_dynamical_decoupling_pi_pulse(self) -> TNone:
        """
        Perform pi pulse so refocusing spins after continuous dynamical decoupling
        :param profile: urukul profile with the correct ampltiude, phase, and frequency settings
        :param time_pi_pulse_mu: length of the pi pulse
        """
        # set to pi pulse parameters
        self.qubit.off()
        self.qubit.singlepass0_on()
        self.qubit.singlepass1_off()
        self.qubit.singlepass2_off()
        self.qubit.cpld.set_all_att_mu(self.att_reg_dynamical_decoupling_pi_pulse)
        self.qubit.set_profile(self.profile_729_pi_pulse)
        at_mu((now_mu() + 8) & ~7)
        self.qubit.io_update()
        self.qubit.on()
        delay_mu(self.time_dynamical_decoupling_pi_pulse_mu)
        self.qubit.off()

    @kernel(flags={'fast-math'})
    def pulse_cat(self, profile: TInt32, time_pulse_mu: TInt64) -> TNone:
        """
        Bichromatic interaction to produce a cat state
        :param profile: urukul profile with the proper settings
        :param time_pulse_mu: how long to apply the bichromatic interaction (which creates the cat state) for
        """
        # set everything back
        self.qubit.off()
        self.qubit.cpld.set_all_att_mu(self.att_reg_cat_interferometer)
        self.qubit.set_profile(profile)
        at_mu((now_mu() + 8) & ~7)
        self.qubit.io_update()
        # turn on all beams
        self.qubit.singlepass1_on()
        self.qubit.singlepass2_on()
        if self.enable_dynamical_decoupling:
            self.qubit.singlepass0_on()
        else:
            self.qubit.singlepass0_off()
        self.qubit.on()
        delay_mu(time_pulse_mu)
        # turn off all beams except singlepass 0 to prevent thermal fluctuations on the singlepass
        self.qubit.off()
        self.qubit.singlepass1_off()
        self.qubit.singlepass2_off()
        self.qubit.singlepass0_on()

    @kernel(flags={'fast_math'})
    def write_pi_pulse_phase(self, phase_pow: TInt32):

        self.qubit.singlepass0.set_mu(
            self.freq_beams_ftw_list[self.profile_729_pi_pulse][1],
            asf=self.ampl_beams_asf_list[self.profile_729_pi_pulse][1],
            pow_=phase_pow,
            profile= self.profile_729_pi_pulse,
            phase_mode=ad9910.PHASE_MODE_CONTINUOUS
        )

    @kernel(flags={"fast-math"})
    def setup_beam_profiles(self) -> TNone:
        """
        Configure parameters for relevant profiles on urukul
        :param time_start_mu: fiducial timestamp for initial start reference (in machine units).
        """
        # ensure all beams are off
        self.qubit.off()
        self.qubit.singlepass0_off()
        self.qubit.singlepass1_off()
        self.qubit.singlepass2_off()
        # set up relevant beam waveforms
        for profile in self.profiles:
            self.qubit.set_mu(
                self.freq_beams_ftw_list[profile][0], asf=self.ampl_beams_asf_list[profile][0],
                pow_=self.phase_beams_pow_list[profile][0], profile=profile,
                phase_mode=ad9910.PHASE_MODE_CONTINUOUS
            )
            self.qubit.singlepass0.set_mu(
                self.freq_beams_ftw_list[profile][1],
                asf=self.ampl_beams_asf_list[profile][1],
                pow_=self.phase_beams_pow_list[profile][1],
                profile=profile,
                phase_mode=ad9910.PHASE_MODE_CONTINUOUS
            )
            self.qubit.singlepass1.set_mu(
                self.freq_beams_ftw_list[profile][2],
                asf=self.ampl_beams_asf_list[profile][2],
                pow_=self.phase_beams_pow_list[profile][2],
                profile=profile,
                phase_mode=ad9910.PHASE_MODE_CONTINUOUS,
            )
            self.qubit.singlepass2.set_mu(
                self.freq_beams_ftw_list[profile][3],
                asf=self.ampl_beams_asf_list[profile][3],
                pow_=self.phase_beams_pow_list[profile][3],
                profile=profile,
                phase_mode=ad9910.PHASE_MODE_CONTINUOUS
            )

    @rpc
    def set_default_profile_configuration(self):

        # set up values consistent across profiles
        for profile in self.profiles:
            self.ampl_beams_asf_list[profile][0] = self.ampl_doublepass_default_asf
            self.phase_beams_pow_list[profile][0] = 0

            self.freq_beams_ftw_list[profile][1] = self.qubit.freq_singlepass0_default_ftw
            self.ampl_beams_asf_list[profile][2] = self.ampls_cat_asf[0]
            self.ampl_beams_asf_list[profile][3] = self.ampls_cat_asf[1]

        # set up values consistent across cat profiles
        for profile in [self.profile_729_cat1a, self.profile_729_cat1b, self.profile_729_cat2a, self.profile_729_cat2b]:
            self.ampl_beams_asf_list[profile][1] = self.ampl_dynamical_decoupling_asf

        # set up values consistent across pi pulse profiles

        self.ampl_beams_asf_list[self.profile_729_pi_pulse][1] = self.ampl_dynamical_decoupling_pi_asf
        self.phase_beams_pow_list[self.profile_729_pi_pulse][2] = 0
        self.phase_beams_pow_list[self.profile_729_pi_pulse][3] = 0

        self.phase_beams_pow_list[self.profile_729_cat1a][2] = self.phases_pulse1_cat_pow[0]
        self.phase_beams_pow_list[self.profile_729_cat1a][3] = self.phases_pulse1_cat_pow[1]

        self.phase_beams_pow_list[self.profile_729_cat1b][2] = (self.phases_pulse1_cat_pow[0] + self.phase_cat_shift_pow)
        self.phase_beams_pow_list[self.profile_729_cat1b][3] = (self.phases_pulse1_cat_pow[1] - self.phase_cat_shift_pow)


    @kernel(flags={'fast-math'})
    def update_profile_configuration(self, freq_cat_center_ftw: TInt32, freq_cat_secular_ftw: TInt32,
                                     phase_cat2_cat_pow: TInt32,
                                     phase_dynamical_decoupling_cat_pow: TInt32) -> TNone:

        """
        Fill out list for easy access of values when setting up the profiles on the urukul
        :param freq_cat_center_ftw: center frequency of the red and blue sidebands
        :param freq_cat_secular_ftw: single-sided detuning of the red and blue sidebands from the center cat frequency
        :param phase_cat2_cat_pow: phase of the second cat pulse
        :param phase_dynamical_decoupling_cat_pow: phase of the dynamical decoupling tone
        """
        self.phase_beams_pow_list[self.profile_729_pi_pulse][1] = (
                self.phase_dynamical_decoupling_pi_pulse_pow_list[0] + phase_dynamical_decoupling_cat_pow)

        # prepare phase arrays for bichromatic
        cat4_phases = [
            self.phases_cat2_cat_update_dir[0] * phase_cat2_cat_pow,
            self.phases_cat2_cat_update_dir[1] * phase_cat2_cat_pow,
        ]

        # set up values consistent across profiles
        for profile in self.profiles:
            self.freq_beams_ftw_list[profile][0] = freq_cat_center_ftw

            self.freq_beams_ftw_list[profile][2] = self.qubit.freq_singlepass1_default_ftw - freq_cat_secular_ftw
            self.freq_beams_ftw_list[profile][3] = self.qubit.freq_singlepass2_default_ftw + freq_cat_secular_ftw

        self.phase_beams_pow_list[self.profile_729_cat1a][1] = phase_dynamical_decoupling_cat_pow

        self.phase_beams_pow_list[self.profile_729_cat1b][1] = phase_dynamical_decoupling_cat_pow

        self.phase_beams_pow_list[self.profile_729_cat2a][1] = phase_dynamical_decoupling_cat_pow
        self.phase_beams_pow_list[self.profile_729_cat2a][2] = cat4_phases[0]
        self.phase_beams_pow_list[self.profile_729_cat2a][3] = cat4_phases[1]
        self.phase_beams_pow_list[self.profile_729_cat2b][1] = phase_dynamical_decoupling_cat_pow
        self.phase_beams_pow_list[self.profile_729_cat2b][2] = cat4_phases[0] + self.phase_cat_shift_pow
        self.phase_beams_pow_list[self.profile_729_cat2b][3] = cat4_phases[1] - self.phase_cat_shift_po


    @kernel(flags={"fast-math"})
    def pulse_readout_sbr(self, time_pulse_mu: TInt64, freq_readout_ftw: TInt32) -> TNone:
        """
        Run a sideband readout pulse.
        :param time_pulse_mu: length of pulse (in machine units).
        :param freq_readout_ftw: readout frequency (set by the double pass) in FTW.
        """
        # set up relevant beam waveforms
        self.qubit.off()
        self.qubit.singlepass1_off()
        self.qubit.singlepass2_off()
        self.qubit.set_profile(self.profile_729_readout)
        self.qubit.io_update()
        delay_mu(8)
        self.qubit.set_mu(freq_readout_ftw, asf=self.ampl_729_readout_asf,
                          pow_=0, profile=self.profile_729_readout,
                          phase_mode=ad9910.PHASE_MODE_CONTINUOUS)
        delay_mu(self.urukul_setup_time_mu)
        self.qubit.singlepass0.set_mu(self.qubit.freq_singlepass0_default_ftw,
                                      asf=self.qubit.ampl_singlepass0_default_asf, pow_=0,
                                      profile=self.profile_729_readout,
                                      phase_mode=ad9910.PHASE_MODE_CONTINUOUS)
        delay_mu(self.urukul_setup_time_mu)
        self.qubit.cpld.set_all_att_mu(self.att_reg_readout_sbr)
        # run readout pulse
        self.qubit.singlepass0_on()
        self.qubit.singlepass1_off()
        self.qubit.singlepass2_off()
        self.qubit.on()
        delay_mu(time_pulse_mu)
        self.qubit.off()

    @kernel(flags={"fast-math"})
    def pulse_readout_rap(self) -> TNone:
        """
        Run a RAP readout pulse.
        """
        # set up relevant beam waveforms
        self.qubit.off()
        self.qubit.singlepass0_on()
        self.qubit.singlepass1_off()
        self.qubit.singlepass2_off()
        self.qubit.singlepass0.set_mu(self.qubit.freq_singlepass0_default_ftw,
                                      asf=self.qubit.ampl_singlepass0_default_asf, pow_=0,
                                      profile=self.profile_729_RAP,
                                      phase_mode=ad9910.PHASE_MODE_CONTINUOUS)
        delay_mu(self.urukul_setup_time_mu)
        self.qubit.cpld.set_all_att_mu(self.att_reg_readout_rap)
        # run RAP readout pulse
        # run rap turns on qubit
        self.rap_subsequence.run_rap(self.time_rap_mu)
