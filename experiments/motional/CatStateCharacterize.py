from artiq.experiment import *
from artiq.coredevice import ad9910

from numpy import copy as np_copy
from numpy import array, int32, int64, arange, zeros, mean

from LAX_exp.language import *
from LAX_exp.system.subsequences import (
    InitializeQubit, SidebandCoolContinuousRAM, Readout, ReadoutAdaptive,
    RescueIon, QubitRAP
)

from LAX_exp.system.objects.SpinEchoWizardRDX import SpinEchoWizardRDX
from LAX_exp.system.objects.PulseShaper import available_pulse_shapes
from LAX_exp.system.objects.PhaserPulseShaper import PhaserPulseShaper, PULSESHAPER_MAX_WAVEFORMS


# todo: ensure phaser pulses are phase-tracking somehow


class CatStateCharacterize(LAXExperiment, Experiment):
    """
    Experiment: Cat State Characterize

    Create and characterize cat states with projective state preparation.
    Uses adaptive readout to reduce timing overheads and extend available coherence times.
    """
    name = 'Cat State Characterize'
    kernel_invariants = {
        # subsequences & objects
        'initialize_subsequence', 'sidebandcool_subsequence', 'readout_subsequence',
        'readout_adaptive_subsequence', 'rescue_subsequence', 'rap_subsequence',
        'spinecho_wizard', 'pulse_shaper',

        # hardware values - cat - default
        'ampl_doublepass_default_asf', 'freq_sigmax_ftw', 'ampl_sigmax_asf', 'time_sigmax_mu',
        'time_herald_slack_mu', 'time_adapt_read_slack_mu', 'max_herald_attempts',

        # hardware values - cat - bichromatic
        'ampls_cat_asf', 'time_cat1_bichromatic_mu', 'phases_pulse1_cat_pow', 'phase_cat1_antisigmax_pow',
        'phase_cat2_sigmax_pow', 'phase_cat2_antisigmax_pow', 'phases_cat2_cat_pow', 'phases_cat2_cat_update_dir',

        # hardware values - QVSA
        'att_phaser_mu', 'freq_phaser_carrier_hz', 'freq_osc_base_hz_list',
        'waveform_index_to_compiled_wav', '_waveform_param_list', '_num_phaser_oscs',

        # hardware values - readout
        'readout_config', 'ampl_729_readout_asf', 'freq_rap_center_ftw', 'freq_rap_dev_ftw', 'time_rap_mu',

        #   hardware values - dynamical decoupling
        'ampl_dynamical_decoupling_asf',
        'att_dynamical_decoupling_mu',

        # configs
        'profile_729_SBC', 'profile_729_RAP',
        'profile_729_target',
        'att_reg_sigmax', 'att_reg_bichromatic', 'att_reg_readout_sbr', 'att_reg_readout_rap',
        'config_experiment_list',
    }

    def build_experiment(self):

        self._num_phaser_oscs = 5  # number of phaser oscillators in use

        # core arguments
        self.setattr_argument("repetitions", NumberValue(default=2, precision=0, step=1, min=1, max=100000))
        self.setattr_argument("readout_type", EnumerationValue(["None", "SBR", "RAP"], default="RAP"),
                              tooltip="None: NO 729nm pulses are applied before state-selective readout.\n"
                                      "SBR (Sideband Ratio): compares RSB and BSB amplitudes.\n"
                                      "RAP (Rapid Adiabatic Passage): Does RAP to measure fock state overlap.\n"
                                      "Note: readout pulses are NOT phase coherent with any bichromatic/sigma_x pulses.")

        # allocate relevant beam profiles
        self.profile_729_SBC = 1
        self.profile_729_RAP = 2
        self.profile_729_target = 6

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
        self.setattr_device('phaser_eggs')

        # instantiate helper objects
        self.spinecho_wizard = SpinEchoWizardRDX(self)

        # set build arguments
        self._build_arguments_default()
        self._build_arguments_dynamical_decoupling()
        self._build_arguments_cat1()
        self._build_arguments_cat2()
        self._build_arguments_qvsa_waveform()
        self._build_arguments_qvsa_sweep()
        self._build_arguments_qvsa_modulation()
        self._build_arguments_readout()

        # instantiate RAP here since it relies on experiment arguments
        self.rap_subsequence = QubitRAP(
            self, ram_profile=self.profile_729_RAP, ram_addr_start=202, num_samples=250,
            ampl_max_pct=self.ampl_rap_pct, pulse_shape="blackman"
        )

    def _build_arguments_default(self):
        """
        Build arguments for default beam parameters.
        """
        # defaults - sigma_x
        self.setattr_argument("freq_sigmax_mhz",
                              NumberValue(default=101.08388, precision=6, step=1, min=50., max=400., scale=1.,
                                          unit="MHz"),
                              group="default.sigmax",
                              tooltip="Frequency for both the sigma_x and anti-sigma_x pulses. "
                                      "Applied to main/chamber doublepass.")
        self.setattr_argument("ampl_sigmax_pct",
                              NumberValue(default=50., precision=3, step=5, min=0.01, max=50, scale=1., unit="%"),
                              group="default.sigmax",
                              tooltip="DDS amplitude for both the sigma_x and anti-sigma_x pulses. "
                                      "Applied to main/chamber doublepass.")
        self.setattr_argument("att_sigmax_db",
                              NumberValue(default=8., precision=1, step=0.5, min=8., max=31.5, scale=1., unit="dB"),
                              group="default.sigmax",
                              tooltip="DDS attenuation for both the sigma_x and anti-sigma_x pulses. "
                                      "Applied to main/chamber doublepass.")
        self.setattr_argument("time_sigmax_us",
                              NumberValue(default=2.38, precision=3, step=0.1, min=0.01, max=10000, scale=1., unit="us"),
                              group="default.sigmax",
                              tooltip="Pulse time for both the sigma_x and anti-sigma_x pulses. "
                                      "Applied to main/chamber doublepass.")

        # defaults - beam values - doublepass (main)
        self.setattr_argument("ampl_doublepass_default_pct",
                              NumberValue(default=50., precision=3, step=5, min=0.01, max=50, scale=1., unit="%"),
                              group="default.cat",
                              tooltip="DDS amplitude for the main doublepass during the bichromatic pulse.")
        self.setattr_argument("att_doublepass_default_db",
                              NumberValue(default=16., precision=1, step=0.5, min=8., max=31.5, scale=1., unit="dB"),
                              group="default.cat",
                              tooltip="DDS attenuation for the main doublepass during the bichromatic pulse.")
        self.setattr_argument("ampls_cat_pct", PYONValue([50., 50.]), group='default.cat',
                              tooltip="DDS amplitudes for the singlepass DDSs during the bichromatic pulses.\n"
                                      "Should be a list of [rsb_ampl_pct, bsb_ampl_pct], which are applied to [singlepass0, singlepass1].")
        self.setattr_argument("atts_cat_db", PYONValue([11., 11.]), group='default.cat',
                              tooltip="DDS attenuations for the singlepass DDSs during the bichromatic pulses.\n"
                                      "Should be a list of [rsb_att_db, bsb_att_db], which are applied to [singlepass0, singlepass1].")
        self.setattr_argument("freq_cat_center_mhz_list", Scannable(
            default=[
                ExplicitScan([101.08388]),
                CenterScan(101.0918, 0.01, 0.0001, randomize=True),
                RangeScan(101.0818, 101.1018, 50, randomize=True),
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
                                      "The singlepass0 DDS is treated as the RSB, and will thus have its frequency DECREASED by this amount.\n"
                                      "Similarly, the singlepass1 DDS is treated as the BSB, and will thus have its frequency INCREASED by this amount.\n"
                                      "i.e. frequencies for [singlepass0, singlepass1] is set as [beams.freq_mhz.freq_singlepass0_mhz - freq_cat_secular_khz, "
                                      "beams.freq_mhz.freq_singlepass1_mhz + freq_cat_secular_khz].")

    def _build_arguments_dynamical_decoupling(self):
        # get arguments for dynamical decoupling
        self.setattr_argument('enable_dynamical_decoupling', BooleanValue(True),
                              tooltip='Indicate whether to apply a third rf tone on the 729 single pass AOM. This '
                                      'tone \n'
                                      'will produce a combination of linear combination of sigma_x and sigma_y \n'
                                      '(dependent on laser phase) which will rotate away any sigma_z errors \n',
                              group='default.dynamical_decoupling')

        self.setattr_argument('freq_dynamical_decoupling_mhz_list',
                              Scannable(default=[ExplicitScan([120.339]), RangeScan(120.1, 120.5, 10),
                                                 CenterScan(120.339, 0.01, 0.001)], unit='MHz',
                                        global_min=80., global_max=120., global_step=0.001, precision=4, scale=1.0),
                              group='default.dynamical_decoupling',
                              tooltip='Frequency of third rf tone to apply to the 729 single pass AOM. This frequency '
                                      '\n'
                                      'should be configured to drive the carrier transition between S1/2 and D5/2 \n'
                                      'and produce a c1*sigma_x + c2*sigma_y operator. \n'
                                      'NOTE WELL: this tone is going to the singlepass so generating a carrier tone \n'
                                      'means outputting a tone at the singlepass default freq listed in the parameter \n'
                                      'settings (\pm any shifts or detunings)')

        self.setattr_argument('ampl_dynamical_decoupling_pct', NumberValue(default=50., max=50., min=0.01,
                                                                           step=5, precision=0, scale=1., unit='%'),
                              tooltip='Amplitude of third rf tone applied to the 729 single pass AOM. \n'
                                      'Stronger tones will further weaken any sigma_z errors BUT \n'
                                      'NOTE: Changing the strength will change AC stark shifts to the sidebands.',
                              group='default.dynamical_decoupling')

        self.setattr_argument('att_dynamical_decoupling_dB', NumberValue(default=28., max=31.5, min=7.,
                                                                         step=0.5, precision=1., scale=1., unit='dB'),
                              tooltip='Attenuation applied to the third rf tone applied to the 729 single pass AOM.',
                              group='default.dynamical_decoupling')

        self.setattr_argument('phase_dynamical_decoupling_cat1_turns_list',
                              Scannable(default=[
                                                # ExplicitScan([0.]),
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

        self.setattr_argument('phase_dynamical_decoupling_cat2_turns_list',
                              Scannable(default=[ExplicitScan([0.]), RangeScan(0, 1, 10),
                                                 CenterScan(0.5, 1, 0.1)], unit='turns',
                                        global_min=0., global_max=2., global_step=0.1, precision=3, scale=1.0),
                              tooltip='Phase of the third rf tone applied to the 729 single pass AOM during the '
                                      'second catting operation. \n'
                                      'For dynamical decoupling to work the phase must be configured correctly as we \n'
                                      'need [H_{bi}, H_{dd}]=0 and this only occurs if these Hamiltonians contain \n'
                                      'the same linear combination of c1*sigma_{x}+c2*sigma_{y} which is determined '
                                      'by \n'
                                      'the laser phase.',
                              group='default.dynamical_decoupling')

    def _build_arguments_cat1(self):
        """
        Build arguments for bichromatic/cat pulse #1.
        """
        # cat #1 config (sigma_x only)
        self.setattr_argument("enable_cat1_sigmax", BooleanValue(default=False), group='cat1.config',
                              tooltip="Applies a sigma_x pulse BEFORE the 1st bichromatic pulse.\n"
                                      "If sigma_x is applied (i.e. True), the bichromatic pulse creates a pure eigenstate (e.g. coherent state).\n"
                                      "If sigma_x is disabled (i.e. False), the bichromatic pulse creates a superposition state (e.g. cat state).")
        self.setattr_argument("enable_cat1_antisigmax", BooleanValue(default=False), group='cat1.config',
                              tooltip="Applies a sigma_x pulse AFTER the 1st bichromatic pulse.\n"
                                      "If the sigma_x pulse is APPLIED, then this pulse disentangles spin from motion.\n"
                                      "If the sigma_x pulse is DISABLED, then this pulse simply selects whether an "
                                      "odd or even superposition is associated with the dark state.")
        self.setattr_argument("phase_cat1_antisigmax_turns",
                              NumberValue(default=0.0, precision=3, step=0.1, min=-1.0, max=1.0, scale=1.,
                                          unit="turns"),
                              group='cat1.config',
                              tooltip="Relative phase applied for the anti-sigma_x pulse.\n"
                                      "Note: this phase is applied via the main doublepass DDS, so values should be halved.")

        # cat #1 config
        self.setattr_argument("enable_cat1_bichromatic", BooleanValue(default=True), group='cat1.config',
                              tooltip="Enables application of the 1st bichromatic pulse.\n"
                                      "Pulses are applied as [sigma_x, bichromatic, antisigma_x, herald, quench].")
        self.setattr_argument("enable_cat1_herald", BooleanValue(default=False), group='cat1.config',
                              tooltip="Enables spin-state heralding via state-selective fluorescence. "
                                      "Heralding only progresses if the state is dark, since otherwise, the motional state is destroyed.\n"
                                      "Pulses are applied as [sigma_x, bichromatic, antisigma_x, herald, quench].\n"
                                      "Note: uses adaptive readout - ensure adaptive readout arguments are correctly set in the dataset manager.")
        self.setattr_argument("enable_cat1_quench", BooleanValue(default=False), group='cat1.config',
                              tooltip="Enables quenching via 854nm to return the spin-state to the S-1/2 state.\n"
                                      "Note: if quench is applied to a superposition state, then the result is a mixed state, not a pure state.\n"
                                      "Pulses are applied as [sigma_x, bichromatic, antisigma_x, herald, quench].")
        self.setattr_argument("time_cat1_bichromatic_us",
                              NumberValue(default=5e6, precision=2, step=5, min=0.1, max=10000, scale=1., unit="us"),
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
        # cat #2 config (sigma_x only)
        self.setattr_argument("enable_cat2_sigmax", BooleanValue(default=False),
                              group='cat2.config',
                              tooltip="Applies a sigma_x pulse BEFORE the 2nd bichromatic pulse.\n"
                                      "If sigma_x is applied (i.e. True), the bichromatic pulse creates a pure eigenstate (e.g. coherent state).\n"
                                      "If sigma_x is disabled (i.e. False), the bichromatic pulse creates a superposition state (e.g. cat state).")
        self.setattr_argument("phase_cat2_sigmax_turns",
                              NumberValue(default=0., precision=3, step=0.1, min=-1.0, max=1.0, scale=1., unit="turns"),
                              group='cat2.config',
                              tooltip="Relative phase applied for the sigma_x pulse.\n"
                                      "Note: this phase is applied via the main doublepass DDS, so values should be halved.")
        self.setattr_argument("enable_cat2_antisigmax", BooleanValue(default=False), group='cat2.config',
                              tooltip="Applies a sigma_x pulse AFTER the 2nd bichromatic pulse.\n"
                                      "If the sigma_x pulse is APPLIED, then this pulse disentangles spin from motion.\n"
                                      "If the sigma_x pulse is DISABLED, then this pulse simply selects whether an "
                                      "odd or even superposition is associated with the dark state.")
        self.setattr_argument("phase_cat2_antisigmax_turns",
                              NumberValue(default=0., precision=3, step=0.1, min=-1.0, max=1.0, scale=1., unit="turns"),
                              group='cat2.config',
                              tooltip="Relative phase applied for the anti-sigma_x pulse.\n"
                                      "Note: this phase is applied via the main doublepass DDS, so values should be halved.")

        # cat2 - config
        self.setattr_argument("enable_cat2_bichromatic", BooleanValue(default=False),
                              group='cat2.config',
                              tooltip="Enables application of the 2nd bichromatic pulse.\n"
                                      "Pulses are applied as [sigma_x, bichromatic, antisigma_x, herald, quench].")
        self.setattr_argument("enable_cat2_herald", BooleanValue(default=False), group='cat2.config',
                              tooltip="Enables spin-state heralding via state-selective fluorescence. "
                                      "Heralding only progresses if the state is dark, since otherwise, the motional state is destroyed.\n"
                                      "Pulses are applied as [sigma_x, bichromatic, antisigma_x, herald, quench].\n"
                                      "Note: uses adaptive readout - ensure adaptive readout arguments are correctly set in the dataset manager.")
        self.setattr_argument("enable_cat2_quench", BooleanValue(default=True), group='cat2.config',
                              tooltip="Enables quenching via 854nm to return the spin-state to the S-1/2 state.\n"
                                      "Note: if quench is applied to a superposition state, then the result is a mixed state, not a pure state.\n"
                                      "Pulses are applied as [sigma_x, bichromatic, antisigma_x, herald, quench].\n")

        # cat2 - pulse parameters
        self.setattr_argument("time_cat2_cat_us_list", Scannable(
            default=[
                ExplicitScan([100]),
                RangeScan(0, 500, 50, randomize=True),
            ],
            global_min=1, global_max=100000, global_step=1,
            unit="us", scale=1, precision=5
        ),
                              group="cat2.bichromatic",
                              tooltip="Pulse time for the 2nd bichromatic pulse.")
        # todo: make this more obviously separate from phase_cat2_cat_turns_list - e.g. phases_cat2_cat_offset_turns
        self.setattr_argument("phases_cat2_cat_turns", PYONValue([0., 0.]), group='cat2.bichromatic',
                              tooltip="Relative phase OFFSET for the singlepass DDSs during the 2nd bichromatic pulse. "
                                      "These phases are applied IN ADDITION to the values from phase_cat2_cat_turns_list.\n"
                                      "Should be a list of [rsb_phase_turns, bsb_phase_turns].\n"
                                      "Note: these phases are applied to the singlepass DDSs, so do not need to be halved, "
                                      "unlike phases applied to the main doublepass.")
        self.setattr_argument("target_cat2_cat_phase",
                              EnumerationValue(['RSB', 'BSB', 'RSB-BSB', 'RSB+BSB'], default='RSB-BSB'),
                              group="cat2.bichromatic",
                              tooltip="Phase update array for the singlepass DDSs during the 2nd bichromatic pulse.\n"
                                      "This configures how phase_cat2_cat_turns_list are to be applied to the DDSs.")
        self.setattr_argument("phase_cat2_cat_turns_list", Scannable(
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
                              NumberValue(default=400., precision=3, min=1, max=1e5, step=1, unit="us", scale=1.),
                              group="read.RAP")

    def _build_arguments_qvsa_sweep(self):
        """
        Build core sweep arguments for the QVSA pulse.
        """
        _argstr = "QVSA"  # string to use for arguments

        # waveform - parameter sweeps
        self.setattr_argument("freq_osc_sweep_arr", PYONValue([-1., 1., 0., 0., 0.]),
                              group="{}.sweep".format(_argstr),
                              tooltip="Defines how oscillator freqs should be scaled for values in freq_osc_sweep_khz_list.\n"
                                      "Indices of freq_osc_sweep_arr correspond to the oscillator number. "
                                      "e.g. [1, -1, 0, 0, 0] will adjust osc_0 by +1x the freq value, and osc_1 by -1x the freq value, with the rest untouched.\n"
                                      "Must be a list of length {:d}.".format(self._num_phaser_oscs))
        self.setattr_argument("freq_osc_sweep_khz_list", Scannable(
            default=[
                ExplicitScan([0.]),
                CenterScan(0., 4, 0.1, randomize=True),
                RangeScan(0., 100.0, 26, randomize=True),
            ],
            global_min=-10000, global_max=10000, global_step=10,
            unit="kHz", scale=1, precision=6
        ),
                              group="{}.sweep".format(_argstr),
                              tooltip="Frequency sweep applied via the phaser oscillators.\n"
                                      "Values for each oscillator are adjusted by the array in freq_osc_sweep_arr.")

        # phaser - phase configuration
        self.setattr_argument("phase_osc_sweep_arr", PYONValue([1., 0., 0., 0., 0.]),
                              group="{}.sweep".format(_argstr),
                              tooltip="Defines how oscillator phases should be adjusted for each value in phase_osc_sweep_turns_list. "
                                      "e.g. [1, -1, 0, 0, 0] will adjust osc_0 by +1x the phase value, and osc_1 by -1x the phase value, with the rest untouched. "
                                      "Must be a list of length {:d}.".format(self._num_phaser_oscs))
        self.setattr_argument("phase_osc_sweep_turns_list", Scannable(
            default=[
                ExplicitScan([0.]),
                RangeScan(0, 1.0, 26, randomize=True),
            ],
            global_min=0.0, global_max=1.0, global_step=1,
            unit="turns", scale=1, precision=5
        ),
                              group="{}.sweep".format(_argstr),
                              tooltip="Phase sweep applied via the phaser oscillators.\n"
                                      "Values for each oscillator are adjusted by the array in phase_osc_sweep_arr.")

    def _build_arguments_qvsa_waveform(self):
        """
        Build core waveform arguments for the QVSA pulse.
        """
        _argstr = "QVSA"  # string to use for arguments

        # waveform - global config
        self.setattr_argument("enable_QVSA_pulse", BooleanValue(default=False),
                              group='{}.global'.format(_argstr),
                              tooltip="Enables the QVSA pulse.")
        self.setattr_argument("freq_phaser_carrier_mhz",
                              NumberValue(default=86., precision=7, step=1, min=0.001, max=4800, unit="MHz", scale=1.),
                              group="{}.global".format(_argstr),
                              tooltip="Phaser output center frequency.\n"
                                      "Note: actual center frequency depends on the devices.phaser.freq_center_mhz dataset argument, "
                                      "which should be manually entered into the dataset manager by the user after "
                                      "configuring the TRF and NCO via e.g. the phaser_configure tool.\n"
                                      "Ensure all values are set correctly.")
        self.setattr_argument("freq_global_offset_mhz",
                              NumberValue(default=0., precision=6, step=1., min=-10., max=10., unit="MHz", scale=1.),
                              group="{}.global".format(_argstr),
                              tooltip="Apply a frequency offset via the phaser oscillators to avoid any DUC/NCO/TRF output spurs.\n"
                                      "Range is limited by the phaser oscillator freq range, i.e. [-10, 10] MHz (includes the frequencies in freq_osc_khz_list).")
        self.setattr_argument("phase_global_ch1_turns",
                              NumberValue(default=0., precision=5, step=0.05, min=-1.1, max=1.1, unit="turns",
                                          scale=1.),
                              group="{}.global".format(_argstr),
                              tooltip="Sets a global CH1 phase via the DUC.\n"
                                      "Note: the eggs.phas_ch1_inherent_turns dataset argument is overridden "
                                      "in this experiment.")

        # waveform - custom specification
        self.setattr_argument("time_heating_us",
                              NumberValue(default=100, precision=2, step=500, min=0.04, max=100000000, unit="us",
                                          scale=1.),
                              group="{}.waveform".format(_argstr),
                              tooltip="Time for the total pulse (excluding pulse shaping)."
                                      "e.g. a time_heating_us of 1ms with 2 segments => a segment time of 500us.")
        self.setattr_argument("freq_osc_khz_list", PYONValue([-702.687, 702.687, 0.000, 0., 0.]),
                              group="{}.waveform".format(_argstr),
                              tooltip="Phaser oscillator frequencies.")
        self.setattr_argument("phase_osc_turns_list", PYONValue([0., 0., 0., 0., 0.]),
                              group="{}.waveform".format(_argstr),
                              tooltip="Relative phases between each phaser oscillator. Applied on both CH0 and CH1.")
        self.setattr_argument("att_phaser_db",
                              NumberValue(default=5., precision=1, step=0.5, min=0, max=31.5, unit="dB", scale=1.),
                              group="{}.waveform".format(_argstr),
                              tooltip="Phaser attenuation to be used for both CH0 and CH1.")
        self.setattr_argument("ampl_osc_frac_list", PYONValue([40., 40., 15., 0., 0.]),
                              group="{}.waveform".format(_argstr),
                              tooltip="Phaser oscillator amplitudes. Applied to both CH0 and CH1.\n"
                                      "Note: CH1 amplitudes will be scaled by the amplitude scaling factors in devices.phaser.ch1.ampl_ch1_osc_scale_arr.")
        self.setattr_argument("phase_osc_ch1_offset_turns", PYONValue([0., 0., 0.5, 0.5, 0.5]),
                              group="{}.waveform".format(_argstr),
                              tooltip="Sets the relative CH1 phase via the phaser oscillators.")

    def _build_arguments_qvsa_modulation(self):
        """
        Build core modulation arguments for the QVSA pulse.
        """
        _argstr = "QVSA"  # string to use for arguments

        # waveform - pulse shaping
        self.setattr_argument("enable_pulse_shaping", BooleanValue(default=True),
                              group='{}.shape'.format(_argstr),
                              tooltip="Applies pulse shaping to the edges of the phaser pulse.\n"
                                      "Note: pulse shaping is applied to each constituent PSK block.")
        self.setattr_argument("type_pulse_shape",
                              EnumerationValue(list(available_pulse_shapes.keys()), default='sine_squared'),
                              group='{}.shape'.format(_argstr),
                              tooltip="Pulse shape type to be used.")
        self.setattr_argument("time_pulse_shape_rolloff_us",
                              NumberValue(default=100, precision=1, step=100, min=0.2, max=100000, unit="us",
                                          scale=1.),
                              group='{}.shape'.format(_argstr),
                              tooltip="Time constant of the pulse shape. This is used for both the pulse rollon AND rolloff.\n"
                                      "e.g. a 1ms main pulse time with 100us time_pulse_shape_rolloff_us will result in a 1ms + 2*100us = 1.2ms total pulse time.\n"
                                      "All constituent PSK blocks will have this pulse time applied.\n"
                                      "Note: DMA issues limit the total number of samples (i.e. time_pulse_shape_rolloff_us * freq_pulse_shape_sample_khz).")
        self.setattr_argument("freq_pulse_shape_sample_khz",
                              NumberValue(default=500, precision=0, step=100, min=1, max=5000, unit="kHz",
                                          scale=1.),
                              group='{}.shape'.format(_argstr),
                              tooltip="Sample rate used for pulse shaping.\n"
                                      "This value is inexact and is fixed at multiples of the phaser oscillator update "
                                      "rate (i.e. 40ns) times the number of oscillators in use.")

        # waveform - PSK (Phase-shift Keying) - general
        self.setattr_argument("enable_phase_shift_keying", BooleanValue(default=False),
                              group="{}.psk".format(_argstr),
                              tooltip="Enable PSK-ing: break the main pulse into individual blocks with different phases.\n"
                                      "Number of PSKs is determined by number of phases in phase_osc<x>_psk_turns. "
                                      "All oscillator PSK schedules must have same length.")
        self.setattr_argument("enable_psk_delay", BooleanValue(default=False),
                              group="{}.psk".format(_argstr),
                              tooltip="Add a delay between PSK pulses where oscillator amplitudes are set to 0. "
                                      "Can be used to create e.g. a Ramsey or DD-type pulse sequence.\n"
                                      "Requires enable_phase_shift_keying to be enabled; otherwise, does nothing.\n"
                                      "Note: prepare/cleanup methods (e.g. set phaser atts, set ext switch) are not called for the delay.")
        self.setattr_argument("time_psk_delay_us",
                              NumberValue(default=200., precision=3, min=1, max=1e5, step=1, unit="us", scale=1.),
                              group="{}.psk".format(_argstr),
                              tooltip="Delay time (in us) delay between PSK pulses. Used for ramsey-ing.\n"
                                      "Note: enable_phase_shift_keying AND enable_psk_delay must be True.")

        # waveform - PSK (Phase-shift Keying) - schedule
        self.setattr_argument("phase_osc0_psk_turns", PYONValue([0., 0.5]), group="{}.psk".format(_argstr),
                              tooltip="PSK phase schedule for osc0.")
        self.setattr_argument("phase_osc1_psk_turns", PYONValue([0., 0.5]), group="{}.psk".format(_argstr),
                              tooltip="PSK phase schedule for osc1.")
        self.setattr_argument("phase_osc2_psk_turns", PYONValue([0., 0.]), group="{}.psk".format(_argstr),
                              tooltip="PSK phase schedule for osc2.")
        self.setattr_argument("phase_osc3_psk_turns", PYONValue([0., 0.]), group="{}.psk".format(_argstr),
                              tooltip="PSK phase schedule for osc3.")
        self.setattr_argument("phase_osc4_psk_turns", PYONValue([0., 0.]), group="{}.psk".format(_argstr),
                              tooltip="PSK phase schedule for osc4.")

    def prepare_experiment(self):
        """
        Prepare & precompute experimental values.
        """
        self._prepare_argument_checks()

        ### MAGIC NUMBERS ###
        self.time_adapt_read_slack_mu = self.core.seconds_to_mu(
            20 * us)  # post-adaptive readout (for heralding) slack to prevent RTIOUnderflow errors
        self.time_herald_slack_mu = self.core.seconds_to_mu(
            150 * us)  # add slack to RTIOCounter only if heralding succeeds
        self.max_herald_attempts = 200  # max number of herald attempts before config is skipped

        # run component preparation
        freq_729_readout_ftw_list, time_729_readout_mu_list = self._prepare_experiment_readout()
        (freq_dynamical_decoupling_ftw_list, phase_dynamical_decoupling_cat1_pow_list,
         phase_dynamical_decoupling_cat2_pow_list) = self._prepare_experiment_dynamical_decoupling()
        (freq_cat_center_ftw_list, freq_cat_secular_ftw_list, time_cat2_cat_mu_list,
         phase_cat2_cat_pow_list, time_ramsey_delay_mu_list) = self._prepare_experiment_cat_general()
        freq_osc_sweep_hz_list, waveform_num_list = self._prepare_experiment_qvsa_general()
        self._prepare_experiment_qvsa_waveform()

        # create experiment config
        self.config_experiment_list = create_experiment_config(
            # cat sweeps
            freq_cat_center_ftw_list, freq_cat_secular_ftw_list,
            time_cat2_cat_mu_list, phase_cat2_cat_pow_list,
            freq_729_readout_ftw_list, time_729_readout_mu_list,
            time_ramsey_delay_mu_list,

            # QVSA sweeps
            freq_osc_sweep_hz_list, waveform_num_list,

            # dynamical decoupling sweeps
            freq_dynamical_decoupling_ftw_list, phase_dynamical_decoupling_cat1_pow_list,
            phase_dynamical_decoupling_cat2_pow_list,
            config_type=float, shuffle_config=True
        )

    def _prepare_experiment_dynamical_decoupling(self):
        freq_dynamical_decoupling_ftw_list = array(
            [self.qubit.singlepass2.frequency_to_ftw(freq_dynamical_decoupling_mhz * MHz) for
             freq_dynamical_decoupling_mhz in self.freq_dynamical_decoupling_mhz_list])

        self.ampl_dynamical_decoupling_asf = self.qubit.singlepass2.amplitude_to_asf(
            self.ampl_dynamical_decoupling_pct / 100.)

        if self.enable_dynamical_decoupling:
            self.att_dynamical_decoupling_mu = self.qubit.singlepass2.cpld.att_to_mu(
                self.att_dynamical_decoupling_dB)
        else:
            self.att_dynamical_decoupling_mu = self.qubit.att_singlepass2_default_mu

        phase_dynamical_decoupling_cat1_pow_list = array(
            [self.qubit.singlepass2.turns_to_pow(phase_dynamical_decoupling_turns) for
             phase_dynamical_decoupling_turns in self.phase_dynamical_decoupling_cat1_turns_list])

        phase_dynamical_decoupling_cat2_pow_list = array(
            [self.qubit.singlepass2.turns_to_pow(phase_dynamical_decoupling_turns) for
             phase_dynamical_decoupling_turns in self.phase_dynamical_decoupling_cat2_turns_list])

        return (freq_dynamical_decoupling_ftw_list, phase_dynamical_decoupling_cat1_pow_list,
                phase_dynamical_decoupling_cat2_pow_list,)

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

        # defaults - sigma_x pulses
        self.freq_sigmax_ftw = self.qubit.frequency_to_ftw(self.freq_sigmax_mhz * MHz)
        self.ampl_sigmax_asf = self.qubit.amplitude_to_asf(self.ampl_sigmax_pct / 100.)
        self.time_sigmax_mu = self.core.seconds_to_mu(self.time_sigmax_us * us)

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
        self.phase_cat1_antisigmax_pow = self.qubit.turns_to_pow(self.phase_cat1_antisigmax_turns)
        self.time_cat1_bichromatic_mu = self.core.seconds_to_mu(self.time_cat1_bichromatic_us * us)
        self.phases_pulse1_cat_pow = [self.qubit.singlepass0.turns_to_pow(phas_pow)
                                      for phas_pow in self.phases_pulse1_cat_turns]

        # inter-cat ramsey delay
        if self.enable_ramsey_delay:
            time_ramsey_delay_mu_list = [self.core.seconds_to_mu(time_delay_us * us)
                                         for time_delay_us in self.time_ramsey_delay_us_list]
        else:
            time_ramsey_delay_mu_list = array([0], dtype=int64)

        # cat2 values
        self.phase_cat2_sigmax_pow = self.qubit.turns_to_pow(self.phase_cat2_sigmax_turns)
        self.phase_cat2_antisigmax_pow = self.qubit.turns_to_pow(self.phase_cat2_antisigmax_turns)
        self.phases_cat2_cat_pow = [self.qubit.singlepass0.turns_to_pow(phas_pow)
                                    for phas_pow in self.phases_cat2_cat_turns]

        if self.enable_cat2_bichromatic:
            time_cat2_cat_mu_list = array([self.core.seconds_to_mu(time_us * us)
                                           for time_us in self.time_cat2_cat_us_list], dtype=int64)
            phase_cat2_cat_pow_list = array([self.qubit.singlepass0.turns_to_pow(phas_pow)
                                             for phas_pow in self.phase_cat2_cat_turns_list], dtype=int32)
        else:
            time_cat2_cat_mu_list = array([0], dtype=int64)
            phase_cat2_cat_pow_list = array([0], dtype=int32)

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
        # attenuation register - sigma_x: singlepasses set to default
        self.att_reg_sigmax = 0x00000000 | (
                (att_to_mu(self.att_sigmax_db * dB) << ((self.qubit.beam.chip_select - 4) * 8)) |
                (self.qubit.att_singlepass0_default_mu << ((self.qubit.singlepass0.chip_select - 4) * 8)) |
                (self.qubit.att_singlepass1_default_mu << ((self.qubit.singlepass1.chip_select - 4) * 8)) |
                (self.qubit.att_singlepass2_default_mu << ((self.qubit.singlepass2.chip_select - 4) * 8))
        )

        # attenuation register - bichromatic: main doublepass set to specified experiment argument value
        self.att_reg_bichromatic = 0x00000000 | (
                (att_to_mu(self.att_doublepass_default_db * dB) << ((self.qubit.beam.chip_select - 4) * 8)) |
                (att_to_mu(self.atts_cat_db[0] * dB) << ((self.qubit.singlepass0.chip_select - 4) * 8)) |
                (att_to_mu(self.atts_cat_db[1] * dB) << ((self.qubit.singlepass1.chip_select - 4) * 8)) |
                (self.att_dynamical_decoupling_mu << ((self.qubit.singlepass2.chip_select - 4) * 8))
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

    def _prepare_experiment_qvsa_general(self):
        """
        Prepare general experiment values for the QVSA pulse.
        :return: tuple of (freq_osc_sweep_hz_list, waveform_num_list)
        """
        # set correct phase delays for field geometries (0.5 for osc_2 for dipole)
        # note: sequence blocks stored as [block_num, osc_num] and store [ampl_pct, phase_turns]
        #   e.g. self.sequence_blocks[2, 5, 0] gives ampl_pct of 5th osc in 2nd block
        # note: create object here instead of build since phase_osc_ch1_offset_turns isn't well-defined until prepare
        self.pulse_shaper = PhaserPulseShaper(self, array(self.phase_osc_ch1_offset_turns))

        # convert values to convenience units
        self.att_phaser_mu = att_to_mu(self.att_phaser_db * dB)
        self.freq_phaser_carrier_hz = (self.freq_phaser_carrier_mhz - self.freq_global_offset_mhz) * MHz

        # format build arguments as numpy arrays with appropriate units
        freq_osc_sweep_hz_list = array(list(self.freq_osc_sweep_khz_list)) * kHz
        self.freq_osc_base_hz_list = array(self.freq_osc_khz_list) * kHz + self.freq_global_offset_mhz * MHz
        self.freq_osc_sweep_arr = array(self.freq_osc_sweep_arr, dtype=float)
        self.phase_osc_sweep_turns_list = list(self.phase_osc_sweep_turns_list)
        self.phase_osc_sweep_arr = array(self.phase_osc_sweep_arr, dtype=float)

        # create waveform parameter sweep config (for use by _prepare_experiment_qvsa_waveform)
        self._waveform_param_list = create_experiment_config(
            self.phase_osc_sweep_turns_list,
            shuffle_config=False, config_type=float
        )
        waveform_num_list = arange(len(self._waveform_param_list))

        # don't apply sweep if QVSA is disabled
        if self.enable_QVSA_pulse:
            return freq_osc_sweep_hz_list, waveform_num_list
        else:
            return [-1], [0]

    def _prepare_experiment_qvsa_waveform(self) -> TNone:
        """
        Calculate waveforms and timings for the QVSA pulse.
        Uses SpinEchoWizard and PhaserPulseShaper objects to simplify waveform compilation.
        """
        '''
        PREPARE WAVEFORM COMPILATION
        '''
        self.waveform_index_to_compiled_wav = list()  # store compiled waveform values
        # note: waveform_index_to_pulseshaper_id is NOT kernel_invariant b/c gets updated in phaser_record
        self.waveform_index_to_pulseshaper_id = zeros(len(self._waveform_param_list),
                                                      dtype=int32)  # store waveform ID linked to DMA sequence

        # calculate block timings and scale amplitudes for ramsey-ing
        num_psk_blocks = len(self.phase_osc0_psk_turns)
        if self.enable_phase_shift_keying:

            # case: PSK w/ramsey delay
            if self.enable_psk_delay:
                num_blocks = 2 * num_psk_blocks - 1
                # note: don't create block_time_list_us here b/c we need to adjust for ramsey delay time sweeps
                block_ampl_scale_list = riffle([1] * num_psk_blocks, [0] * (num_psk_blocks - 1))

            # case: PSK only, no ramsey delay
            else:
                num_blocks = num_psk_blocks
                block_time_list_us = [self.time_heating_us / num_psk_blocks] * num_psk_blocks
                block_ampl_scale_list = [1] * num_psk_blocks

        # case: rabi-style (no PSK, no ramsey delay)
        else:
            num_blocks = 1
            block_time_list_us = [self.time_heating_us]
            block_ampl_scale_list = [1]

        '''
        DESIGN WAVEFORM SEQUENCE
        '''
        # create bare waveform block sequence & set amplitudes
        _osc_vals_blocks = zeros((num_blocks, self._num_phaser_oscs, 2), dtype=float)
        _osc_vals_blocks[:, :, 0] = array(self.ampl_osc_frac_list)
        _osc_vals_blocks[:, :, 0] *= array([block_ampl_scale_list]).transpose()

        # set bsb phase and account for oscillator delay time
        # note: use mean of osc freqs since I don't want to record a waveform for each osc freq
        freq_osc_sweep_avg_hz = mean(list(self.freq_osc_sweep_khz_list)) * kHz
        t_update_delay_s_list = array([0, 40e-9, 80e-9, 80e-9, 120e-9])[:self._num_phaser_oscs]
        phase_osc_update_delay_turns_list = (
                (self.freq_osc_base_hz_list +
                 self.freq_osc_sweep_arr * freq_osc_sweep_avg_hz) *
                t_update_delay_s_list
        _osc_vals_blocks[:, :, 1] += array(self.phase_osc_turns_list) + phase_osc_update_delay_turns_list

        # set PSK phase update schedule
        if self.enable_phase_shift_keying:
            if self.enable_psk_delay:
                # note: use ::2 since we only update to non-delay blocks
                _osc_vals_blocks[::2, :, 1] += array([self.phase_osc0_psk_turns, self.phase_osc1_psk_turns,
                                                      self.phase_osc2_psk_turns, self.phase_osc3_psk_turns,
                                                      self.phase_osc4_psk_turns][:self._num_phaser_oscs]).transpose()
            else:
                _osc_vals_blocks[:, :, 1] += array([self.phase_osc0_psk_turns, self.phase_osc1_psk_turns,
                                                    self.phase_osc2_psk_turns, self.phase_osc3_psk_turns,
                                                    self.phase_osc4_psk_turns][:self._num_phaser_oscs]).transpose()

        '''
        COMPILE WAVEFORMS SPECIFIC TO EACH PARAMETER
        '''
        # record phaser waveforms - one for each phase
        for waveform_params in self._waveform_param_list:
            # extract waveform params
            phase_sweep_turns = waveform_params[0]

            # waveform sweep (case - PSK + ramsey): create block_time_list_us with target delay time
            if self.enable_phase_shift_keying and self.enable_psk_delay:
                block_time_list_us = riffle([self.time_heating_us / num_psk_blocks] * num_psk_blocks,
                                            [self.time_psk_delay_us] * (num_psk_blocks - 1))

            # create local copy of _osc_vals_blocks and update with waveform parameters
            # note: no need to deep copy b/c it's filled w/immutables
            _osc_vals_blocks_local = np_copy(_osc_vals_blocks)
            _osc_vals_blocks_local[:, :, 1] += self.phase_osc_sweep_arr * phase_sweep_turns  # apply phase sweep

            # specify sequence as a list of blocks, where each block is a dict
            # note: have to instantiate locally each loop b/c dicts aren't deep copied
            _sequence_blocks_local = [
                {
                    "oscillator_parameters": _osc_vals_blocks_local[_idx_block],
                    "config": {
                        "time_us": block_time_list_us[_idx_block],
                        # don't pulse shape for delay blocks lmao
                        "pulse_shaping": self.enable_pulse_shaping and (block_ampl_scale_list[_idx_block] != 0),
                        "pulse_shaping_config": {
                            "pulse_shape": self.type_pulse_shape,
                            "pulse_shape_rising": self.enable_pulse_shaping,
                            "pulse_shape_falling": self.enable_pulse_shaping,
                            "sample_rate_khz": self.freq_pulse_shape_sample_khz,
                            "rolloff_time_us": self.time_pulse_shape_rolloff_us
                        }
                    }
                } for _idx_block in range(num_blocks)
            ]

            # create QVSA waveform and store data in a holder
            self.waveform_index_to_compiled_wav.append(
                self.spinecho_wizard.compile_waveform(_sequence_blocks_local))

    def _prepare_argument_checks(self) -> TNone:
        """
        Check experiment arguments for validity.
        """
        # todo
        '''
        BICHROMATIC/CAT CHECKS
        '''
        # self.max_ampl_singlepass_pct, self.min_att_singlepass_db = (50., 7.)
        #
        # # ensure single pass values are safe and valid
        # if any((ampl_pct > self.max_ampl_singlepass_pct or ampl_pct < 0.
        #         for ampl_pct in self.ampl_singlepass_default_pct_list)):
        #     raise ValueError(
        #         "Singlepass amplitude outside valid range - [0., {:f}].".format(self.max_ampl_singlepass_pct))
        #
        # if any((att_db > 31.5 or att_db < self.min_att_singlepass_db
        #         for att_db in self.att_singlepass_default_db_list)):
        #     raise ValueError(
        #         "Singlepass attenuation outside valid range - [{:.1f}, 31.5].".format(self.min_att_singlepass_db))

        '''
        CHECK PHASER BASE OSC CONFIG
        '''
        # check that phaser oscillator amplitude config is valid
        if ((not isinstance(self.ampl_osc_frac_list, list)) or
                (len(self.ampl_osc_frac_list) != self._num_phaser_oscs)):
            raise ValueError(
                "Error: phaser oscillator amplitude array must be list of length {:d}.".format(self._num_phaser_oscs))
        elif sum(self.ampl_osc_frac_list) >= 100.:
            raise ValueError("Error: phaser oscillator amplitudes must sum <100.")

        # check that phaser oscillator phase arrays are valid
        if ((not isinstance(self.phase_osc_turns_list, list)) or
                (len(self.phase_osc_turns_list) != self._num_phaser_oscs)):
            raise ValueError(
                "Error: phaser oscillator phase array must be list of length {:d}.".format(self._num_phaser_oscs))

        # check that phaser oscillator frequencies are valid
        if ((not isinstance(self.freq_osc_khz_list, list)) or
                (len(self.freq_osc_khz_list) != self._num_phaser_oscs)):
            raise ValueError(
                "Error: phaser oscillator frequency array must be list of length {:d}.".format(self._num_phaser_oscs))
        max_osc_freq_hz = (max(list(self.freq_osc_sweep_khz_list)) * kHz +
                           max(self.freq_osc_khz_list) * kHz +
                           (self.freq_global_offset_mhz * MHz))
        min_osc_freq_hz = (min(list(self.freq_osc_sweep_khz_list)) * kHz +
                           min(self.freq_osc_khz_list) * kHz +
                           (self.freq_global_offset_mhz * MHz))
        if (max_osc_freq_hz > 12.5 * MHz) or (min_osc_freq_hz < -12.5 * MHz):
            raise ValueError("Error: phaser oscillator frequencies outside valid range of [-12.5, 12.5] MHz.")

        # ensure phaser output frequency falls within valid DUC bandwidth
        phaser_carrier_freq_dev_hz = abs(self.phaser_eggs.freq_center_hz - self.freq_phaser_carrier_mhz * MHz)
        if phaser_carrier_freq_dev_hz >= 300. * MHz:
            raise ValueError("Invalid argument: output frequencies outside +/- 300 MHz phaser DUC bandwidth.")

        '''
        CHECK PHASER WAVEFORM CONFIG
        '''
        # check that PSK schedule is valid
        num_psk_blocks = len(self.phase_osc0_psk_turns)
        psk_schedule_invalid = self.enable_phase_shift_keying and any([
            (not isinstance(psk_schedule, list)) or (len(psk_schedule) != num_psk_blocks)
            for psk_schedule in (
                self.phase_osc0_psk_turns, self.phase_osc1_psk_turns,
                self.phase_osc2_psk_turns, self.phase_osc3_psk_turns
            )
        ])
        if psk_schedule_invalid:
            raise ValueError("Invalid PSK schedule: all PSK schedules must be of same length.")

        # ensure that sweep targets are lists of appropriate length
        if not (isinstance(self.phase_osc_sweep_arr, list) and (
                len(self.phase_osc_sweep_arr) == self._num_phaser_oscs)):
            raise ValueError(
                "Invalid phase_osc_sweep_arr: {:}.\nphase_osc_sweep_arr must be list of length {:d}.".format(
                    self.phase_osc_sweep_arr, self._num_phaser_oscs))
        if not (isinstance(self.freq_osc_sweep_arr, list) and (len(self.freq_osc_sweep_arr) == self._num_phaser_oscs)):
            raise ValueError("Invalid freq_osc_sweep_arr: {:}.\nfreq_osc_sweep_arr must be list of length {:d}.".format(
                self.freq_osc_sweep_arr, self._num_phaser_oscs))

        # check that waveforms are not too many/not sweeping too hard
        num_waveforms_to_record = (len(list(self.phase_osc_sweep_turns_list)))
        if num_waveforms_to_record > PULSESHAPER_MAX_WAVEFORMS:
            raise ValueError("Too many waveforms to record ({:d}) - must be fewer than {:d}.\n"
                             "Reduce length of any of [phase_osc_sweep_turns_list].".format(
                num_waveforms_to_record, PULSESHAPER_MAX_WAVEFORMS))

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_experiment_list),
                13)

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

        # record phaser waveforms
        self.phaser_record()

        # set maximum attenuations for phaser outputs to prevent leakage during configuration
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.channel[0].set_att_mu(0x00)
        delay_mu(self.phaser_eggs.t_sample_mu)
        self.phaser_eggs.channel[1].set_att_mu(0x00)

    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        # predeclare variables ahead of time
        time_start_mu = now_mu() & ~0x7  # store reference time for device synchronization
        ion_state = (-1, 0, int64(0))  # store ion state for adaptive readout
        herald_counter = 0  # store herald attempts
        _loop_iter = 0  # used to check_termination more frequently

        self.pulse_shaper.waveform_load()  # load phaser waveforms from DMA

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
                freq_sweep_hz = config_vals[7]
                waveform_num = int32(config_vals[8])
                freq_dynamical_decoupling_ftw = int32(config_vals[9])
                phase_dynamical_decoupling_cat1_pow = int32(config_vals[10])
                phase_dynamical_decoupling_cat2_pow = int32(config_vals[11])

                herald_counter = 0  # clear herald counter

                # prepare phase arrays for bichromatic
                cat4_phases = [
                    self.phases_cat2_cat_pow[0] + self.phases_cat2_cat_update_dir[0] * phase_cat2_cat_pow,
                    self.phases_cat2_cat_pow[1] + self.phases_cat2_cat_update_dir[1] * phase_cat2_cat_pow,
                ]

                # get corresponding waveform parameters and pulseshaper ID from the index
                waveform_params = self._waveform_param_list[waveform_num]
                phaser_waveform = self.waveform_index_to_pulseshaper_id[waveform_num]

                # configure phaser
                if self.enable_QVSA_pulse:
                    # create frequency update list for phaser oscs and set phaser frequencies
                    freq_update_list = self.freq_osc_base_hz_list + (freq_sweep_hz * self.freq_osc_sweep_arr)
                    # set phaser frequency/phases
                    self.core.break_realtime()
                    self.phaser_eggs.frequency_configure(
                        self.freq_phaser_carrier_hz,  # carrier frequency (via DUC)
                        # oscillator frequencies
                        [freq_update_list[0], freq_update_list[1], freq_update_list[2],
                         freq_update_list[3], freq_update_list[4]],
                        self.phase_global_ch1_turns  # global CH1 phase
                    )

                '''
                BEGIN MAIN SEQUENCE
                '''
                while True:
                    # check heralding OK (otherwise execution is blocked)
                    if herald_counter >= self.max_herald_attempts:
                        print("\t\tWarning: too many heralds. Moving onto next configuration.")
                        self.check_termination()  # check termination b/c we haven't in a while
                        self.core.break_realtime()  # add slack
                        break

                    self.core.break_realtime()  # add slack for execution
                    delay_mu(125000)  # add even more slack lol

                    '''
                    INITIALIZE ION STATE
                    '''
                    # initialize ion in S-1/2 state & SBC to ground state
                    self.initialize_subsequence.run_dma()
                    self.sidebandcool_subsequence.run_dma()

                    # set target profile to ensure we run correctly
                    self.qubit.set_profile(self.profile_729_target)
                    self.qubit.io_update()

                    # synchronize start time to phaser's 320ns frame (which is multiple of coarse RTIO clk)
                    time_start_mu = self.phaser_eggs.get_next_frame_mu()
                    # most important: clear phaser osc HERE & NOW to ensure phase coherent w/ DDSs
                    #   b/c DDSs are phase-tracked wrt time_start_mu
                    self.phaser_eggs.phase_osc_clear()
                    # unimportant: clear DUC phase (b/c why not)
                    self.phaser_eggs.reset_duc_phase()

                    '''
                    CAT #1
                    '''
                    # cat1 - sigma_x (displacement vs cat)
                    if self.enable_cat1_sigmax:
                        self.pulse_sigmax(time_start_mu, 0)
                    # cat1 - bichromatic cat pulse
                    if self.enable_cat1_bichromatic:
                        if self.enable_dynamical_decoupling:
                            self.run_dynamical_decoupling(freq_dynamical_decoupling_ftw,
                                                          phase_dynamical_decoupling_cat1_pow,
                                                          time_start_mu)
                        self.pulse_bichromatic(time_start_mu, self.time_cat1_bichromatic_mu,
                                               self.phases_pulse1_cat_pow,
                                               freq_cat_center_ftw, freq_cat_secular_ftw)
                        self.stop_dynamical_decoupling()
                    # cat1 - anti-sigma_x (return to S-state)
                    if self.enable_cat1_antisigmax:
                        self.pulse_sigmax(time_start_mu, self.phase_cat1_antisigmax_pow)

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
                    if self.enable_dynamical_decoupling:
                        self.run_dynamical_decoupling(freq_dynamical_decoupling_ftw, phase_dynamical_decoupling_cat1_pow,
                                                  time_start_mu)
                    if self.enable_ramsey_delay:
                        delay_mu(time_ramsey_delay_mu)

                    '''
                    QVSA PULSE
                    '''
                    if self.enable_QVSA_pulse:
                        self.phaser_run(phaser_waveform)
                    self.stop_dynamical_decoupling()

                    '''
                    CAT #2
                    '''
                    # cat2 - sigma_x (displacement vs cat)
                    if self.enable_cat2_sigmax:
                        self.pulse_sigmax(time_start_mu, self.phase_cat2_sigmax_pow)
                    # cat2 - bichromatic cat pulse
                    if self.enable_cat2_bichromatic:
                        self.pulse_bichromatic(time_start_mu, time_cat2_cat_mu,
                                               cat4_phases,
                                               freq_cat_center_ftw, freq_cat_secular_ftw)
                        if self.enable_dynamical_decoupling:
                            self.run_dynamical_decoupling(freq_dynamical_decoupling_ftw,
                                                          phase_dynamical_decoupling_cat2_pow,
                                                          time_start_mu)
                        self.stop_dynamical_decoupling()
                    # cat2 - anti-sigma_x (return to S-state)
                    if self.enable_cat2_antisigmax:
                        self.pulse_sigmax(time_start_mu, self.phase_cat2_antisigmax_pow)

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
                    counts_res = self.readout_subsequence.fetch_count()
                else:
                    # return -1 so user knows booboo happened
                    counts_res = -1

                # store results
                self.rescue_subsequence.resuscitate()
                self.update_results(freq_cat_center_ftw,
                                    counts_res,
                                    freq_cat_secular_ftw,
                                    time_cat2_cat_mu,
                                    phase_cat2_cat_pow,
                                    freq_729_readout_ftw,
                                    time_729_readout_mu,
                                    time_ramsey_delay_mu,
                                    freq_sweep_hz,
                                    freq_dynamical_decoupling_ftw,
                                    phase_dynamical_decoupling_cat1_pow,
                                    phase_dynamical_decoupling_cat2_pow,
                                    waveform_params[0])

                # check termination more frequently in case reps are low
                if _loop_iter % 200 == 0:
                    self.check_termination()
                _loop_iter += 1

            # rescue ion as needed & support graceful termination
            self.core.break_realtime()
            self.rescue_subsequence.run(trial_num)
            self.check_termination()

    '''
    HELPER FUNCTIONS
    '''

    @kernel(flags={"fast-math"})
    def pulse_sigmax(self, time_start_mu: TInt64, phas_pow: TInt32) -> TNone:
        """
        Run a phase-coherent sigma_x pulse on the qubit.
        :param time_start_mu: fiducial timestamp for initial start reference (in machine units).
        :param phas_pow: relative phase offset for the beam.
        """
        # set up relevant beam waveforms
        # todo: maybe add like 32ns between these? so we don't booboo?
        self.qubit.set_mu(
            self.freq_sigmax_ftw, asf=self.ampl_sigmax_asf, pow_=phas_pow,
            profile=self.profile_729_target,
            phase_mode=ad9910.PHASE_MODE_TRACKING, ref_time_mu=time_start_mu
        )
        self.qubit.singlepass0.set_mu(
            self.qubit.freq_singlepass0_default_ftw,
            asf=self.qubit.ampl_singlepass0_default_asf, pow_=0,
            profile=self.profile_729_target,
            phase_mode=ad9910.PHASE_MODE_TRACKING, ref_time_mu=time_start_mu
        )
        self.qubit.singlepass1.set_mu(
            self.qubit.freq_singlepass1_default_ftw,
            asf=self.qubit.ampl_singlepass1_default_asf, pow_=0,
            profile=self.profile_729_target,
            phase_mode=ad9910.PHASE_MODE_TRACKING, ref_time_mu=time_start_mu
        )
        self.qubit.cpld.set_all_att_mu(self.att_reg_sigmax)

        # run sigmax pulse
        self.qubit.singlepass0.sw.on()
        self.qubit.singlepass1.sw.off()
        self.qubit.on()
        delay_mu(self.time_sigmax_mu)
        self.qubit.off()
        self.qubit.singlepass1.sw.off()

    @kernel(flags={"fast-math"})
    def pulse_bichromatic(self, time_start_mu: TInt64, time_pulse_mu: TInt64, phas_pow_list: TList(TInt32),
                          freq_carrier_ftw: TInt32, freq_secular_ftw: TInt32) -> TNone:
        """
        Run a phase-coherent bichromatic pulse on the qubit.
        :param time_start_mu: fiducial timestamp for initial start reference (in machine units).
        :param time_pulse_mu: length of pulse (in machine units).
        :param phas_pow_list: relative phase offset for the beams (RSB, BSB) (in pow).
        :param freq_carrier_ftw: carrier frequency (set by the double pass) in FTW.
        :param freq_secular_ftw: bichromatic separation frequency (from central frequency) in FTW.
        """
        # set up relevant beam waveforms
        # todo: maybe add like 32ns between these? so we don't booboo?
        self.qubit.set_mu(
            freq_carrier_ftw, asf=self.ampl_doublepass_default_asf,
            pow_=0, profile=self.profile_729_target,
            phase_mode=ad9910.PHASE_MODE_TRACKING, ref_time_mu=time_start_mu
        )
        self.qubit.singlepass0.set_mu(
            self.qubit.freq_singlepass0_default_ftw - freq_secular_ftw,
            asf=self.ampls_cat_asf[0],
            pow_=phas_pow_list[0],
            profile=self.profile_729_target,
            phase_mode=ad9910.PHASE_MODE_TRACKING, ref_time_mu=time_start_mu
        )
        self.qubit.singlepass1.set_mu(
            self.qubit.freq_singlepass1_default_ftw + freq_secular_ftw,
            asf=self.ampls_cat_asf[1],
            pow_=phas_pow_list[1],
            profile=self.profile_729_target,
            phase_mode=ad9910.PHASE_MODE_TRACKING, ref_time_mu=time_start_mu
        )
        self.qubit.cpld.set_all_att_mu(self.att_reg_bichromatic)

        # run bichromatic pulse
        self.qubit.singlepass0.sw.on()
        self.qubit.singlepass1.sw.on()
        self.qubit.on()
        delay_mu(time_pulse_mu)
        self.qubit.off()
        self.qubit.singlepass1.sw.off()

    @kernel(flags={"fast-math"})
    def pulse_readout_sbr(self, time_pulse_mu: TInt64, freq_readout_ftw: TInt32) -> TNone:
        """
        Run a sideband readout pulse.
        :param time_pulse_mu: length of pulse (in machine units).
        :param freq_readout_ftw: readout frequency (set by the double pass) in FTW.
        """
        # set up relevant beam waveforms
        # todo: maybe add like 32ns between these? so we don't booboo?
        self.qubit.set_mu(freq_readout_ftw, asf=self.ampl_729_readout_asf,
                          pow_=0, profile=self.profile_729_target,
                          phase_mode=ad9910.PHASE_MODE_CONTINUOUS)
        self.qubit.singlepass0.set_mu(self.qubit.freq_singlepass0_default_ftw,
                                      asf=self.qubit.ampl_singlepass0_default_asf, pow_=0,
                                      profile=self.profile_729_target,
                                      phase_mode=ad9910.PHASE_MODE_CONTINUOUS)
        self.qubit.singlepass1.set_mu(self.qubit.freq_singlepass1_default_ftw,
                                      asf=self.qubit.ampl_singlepass1_default_asf, pow_=0,
                                      profile=self.profile_729_target,
                                      phase_mode=ad9910.PHASE_MODE_CONTINUOUS)
        self.qubit.singlepass2.set_mu(self.qubit.freq_singlepass2_default_ftw,
                                      asf=self.qubit.ampl_singlepass2_default_asf, pow_=0,
                                      profile=self.profile_729_RAP,
                                      phase_mode=ad9910.PHASE_MODE_CONTINUOUS)
        self.qubit.cpld.set_all_att_mu(self.att_reg_readout_sbr)

        # run readout pulse
        self.qubit.singlepass0.sw.on()
        # todo: this should be off, right???
        self.qubit.singlepass1.sw.on()
        self.qubit.singlepass2.sw.off()
        self.qubit.on()
        delay_mu(time_pulse_mu)
        self.qubit.off()
        self.qubit.singlepass1.sw.off()

    @kernel(flags={"fast-math"})
    def pulse_readout_rap(self) -> TNone:
        """
        Run a RAP readout pulse.
        """
        # set up relevant beam waveforms
        self.qubit.singlepass0.set_mu(self.qubit.freq_singlepass0_default_ftw,
                                      asf=self.qubit.ampl_singlepass0_default_asf, pow_=0,
                                      profile=self.profile_729_RAP,
                                      phase_mode=ad9910.PHASE_MODE_CONTINUOUS)
        self.qubit.singlepass1.set_mu(self.qubit.freq_singlepass1_default_ftw,
                                      asf=self.qubit.ampl_singlepass1_default_asf, pow_=0,
                                      profile=self.profile_729_RAP,
                                      phase_mode=ad9910.PHASE_MODE_CONTINUOUS)
        self.qubit.singlepass2.set_mu(self.qubit.freq_singlepass2_default_ftw,
                                      asf=self.qubit.ampl_singlepass2_default_asf, pow_=0,
                                      profile=self.profile_729_RAP,
                                      phase_mode=ad9910.PHASE_MODE_CONTINUOUS)
        self.qubit.cpld.set_all_att_mu(self.att_reg_readout_rap)

        # run RAP readout pulse
        self.rap_subsequence.run_rap(self.time_rap_mu)

    @kernel(flags={'fast-math'})
    def run_dynamical_decoupling(self, freq_dynamical_decoupling_ftw: TInt32, phase_pow: TInt32,
                                 time_start_mu: TInt64) -> TNone:
        self.qubit.singlepass2.set_mu(
            freq_dynamical_decoupling_ftw,
            asf=self.ampl_dynamical_decoupling_asf, pow_=phase_pow,
            profile=self.profile_729_target,
            phase_mode=ad9910.PHASE_MODE_TRACKING, ref_time_mu=time_start_mu
        )
        # ensure we are at the right attenuation in case it has been overwritten by other sequences
        self.qubit.singlepass2.set_att_mu(self.att_dynamical_decoupling_mu)
        self.qubit.singlepass2.sw.on()

    @kernel(flags={'fast-math'})
    def stop_dynamical_decoupling(self) -> TNone:
        self.qubit.singlepass2.set_mu(
            self.qubit.freq_singlepass2_default_ftw,
            asf=self.qubit.ampl_singlepass2_default_asf, pow_=0,
            profile=self.profile_729_target,
            phase_mode=ad9910.PHASE_MODE_CONTINUOUS,
        )
        self.qubit.singlepass2.set_att_mu(self.qubit.att_singlepass2_default_mu)
        self.qubit.singlepass2.sw.off()

    '''
    HELPER FUNCTIONS - PHASER
    '''

    @kernel(flags={"fast-math"})
    def phaser_run(self, waveform_id: TInt32) -> TNone:
        """
        Run the main phaser pulse together with supporting functionality.
        :param waveform_id: the ID of the waveform to run.
        """
        # EGGS - START/SETUP
        self.phaser_eggs.phaser_setup(self.att_phaser_mu, self.att_phaser_mu)

        # EGGS - RUN
        # reset DUC phase to start DUC deterministically
        self.phaser_eggs.reset_duc_phase()
        self.pulse_shaper.waveform_playback(waveform_id)

        # EGGS - STOP
        # stop all output & clean up hardware (e.g. eggs amp switches, RF integrator hold)
        # note: DOES unset attenuators (beware turn-on glitch if no filters/switches)
        self.phaser_eggs.phaser_stop()

    @kernel(flags={"fast-math"})
    def phaser_record(self) -> TNone:
        """
        Set up core phaser functionality and record the pulse-shaped waveforms.
        Should be run during initialize_experiment.
        """
        # record phaser sequences onto DMA for each waveform parameter
        for _wav_idx in range(len(self._waveform_param_list)):
            # get waveform for given parameters
            # note: use sync RPC to reduce significant overhead of direct data transfer
            _wav_data_ampl, _wav_data_phas, _wav_data_time = self._get_compiled_waveform(_wav_idx)

            # record phaser pulse sequence and save returned waveform ID
            # note: no need to add slack b/c waveform_record does it for us
            self.waveform_index_to_pulseshaper_id[_wav_idx] = self.pulse_shaper.waveform_record(
                _wav_data_ampl, _wav_data_phas, _wav_data_time
            )

    @rpc
    def _get_compiled_waveform(self, wav_idx: TInt32) -> TTuple([TArray(TFloat, 2),
                                                                 TArray(TFloat, 2),
                                                                 TArray(TInt64, 1)]):
        """
        Return compiled waveform values.
        By returning the large waveform arrays via RPC, we avoid all-at-once large data transfers,
            speeding up experiment compilation and transfer to Kasli.
        :param wav_idx: the index of the compiled waveform to retrieve.
        :return: the compiled waveform corresponding to the waveform index.
        """
        return self.waveform_index_to_compiled_wav[wav_idx]

    '''
    ANALYSIS
    '''

    def analyze_experiment(self):
        pass
