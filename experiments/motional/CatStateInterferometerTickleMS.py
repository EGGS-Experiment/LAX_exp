from artiq.experiment import *
from artiq.coredevice import ad9910

from numpy import array, int32, int64, zeros

from LAX_exp.language import *
from LAX_exp.system.subsequences import (
    InitializeQubit, SidebandCoolContinuousRAM, Readout, ReadoutAdaptive,
    RescueIon, QubitRAP
)

from LAX_exp.system.objects.PulseShaper import available_pulse_shapes
from LAX_exp.system.objects.dds_pulse_shaper import DDSPulseShaper


class CatStateInterferometerTickleMS(LAXExperiment, Experiment):
    """
    Experiment: Cat State Interferometer Tickle MS

    Create and characterize cat states with projective state preparation.
    Uses adaptive readout to reduce timing overheads and extend available coherence times.
    """
    name = 'Cat State Inteferometer Tickle MS'
    kernel_invariants = {
        # subsequences & objects
        'initialize_subsequence', 'sidebandcool_subsequence', 'readout_subsequence',
        'readout_adaptive_subsequence', 'rescue_subsequence', 'rap_subsequence', 'dds_pulse_shaper',

        # ion parameters
        'freq_secular_ftw',
    
        # hardware values - ms - bichromatic
        'enable_ms_gate', 'ampls_ms_asf',
                                                       
        # hardware values - parity
        'enable_parity_pulse', 'time_parity_pulse_mu',

        # hardware values - cat - default
        'time_herald_slack_mu', 'time_adapt_read_slack_mu', 'max_herald_attempts', 'ampls_cat_asf',

        # hardware values - cat1 - bichromatic
        'enable_cat1_bichromatic', 'enable_cat1_herald', 'enable_cat1_quench',
        'time_cat_bichromatic_mu', 'phases_pulse1_cat_pow',

        # hardware values - ramsey
        'enable_ramsey_delay',

        # hardware values - cat2 - bichromatic
        'enable_cat2_bichromatic', 'enable_cat2_herald', 'enable_cat2_quench',
        'phase_cat_update_dir',

        # hardware values - tickle
        'enable_tickle_pulse', 'att_tickle_mu', 'time_tickle_mu',

        # hardware values - readout
        'readout_config', 'freq_rap_center_ftw', 'freq_rap_dev_ftw', 'time_rap_mu',

        #   hardware values - dynamical decoupling
        'enable_dynamical_decoupling', 'ampl_dynamical_decoupling_asf',

        # hardware values - intensity servo
        'enable_servo_relock', 'time_servo_relock_mu',

        # configs
        'profile_729_SBC',
        'profile_729_cat1', 'profile_729_cat2',
        'profile_729_RAP', 'profile_tickle_RAM',
        'profile_729_ms', 'profile_729_parity',
        'att_reg_cat_interferometer', 'att_reg_readout_rap',
        'att_reg_ms_gate', 'att_reg_parity_pulse',
        'config_experiment_list',

        # extras
        'urukul_setup_time_mu', 'profiles', 'urukul_dd_reset_time',
    }

    def build_experiment(self):

        # core arguments
        self.setattr_argument("repetitions", NumberValue(default=50, precision=0, step=1, min=1, max=100000))
        self.setattr_argument("readout_type", EnumerationValue(["None",  "RAP"], default="RAP"),
                              tooltip="None: NO 729nm pulses are applied before state-selective readout.\n"
                                      "RAP (Rapid Adiabatic Passage): Does RAP to measure fock state overlap.\n"
                                      "Note: readout pulses are NOT phase coherent with any bichromatic/sigma_x pulses.")

        # allocate relevant beam profiles
        self.profile_729_RAP = 0
        self.profile_729_SBC = 1
        self.profile_729_cat1 = 2
        self.profile_729_cat2 = 3
        self.profile_729_ms = 4
        self.profile_729_parity = 5

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
        self.setattr_device('dds_dipole')

        # set build arguments
        self._build_arguments_ion_parameters()
        self._build_arguments_dynamical_decoupling()
        self._build_arguments_ms_gate()
        self._build_arguments_parity()
        self._build_arguments_cat_default()
        self._build_arguments_cat1()
        self._build_arguments_ramsey()
        self._build_arguments_cat2()
        self._build_arguments_tickle()
        self._build_arguments_readout()
        self._build_arguments_intensity_servo()

        # # instantiate RAP here since it relies on experiment arguments
        self.rap_subsequence = QubitRAP(
            self, ram_profile=self.profile_729_RAP, ram_addr_start=202, num_samples=250,
            ampl_max_pct=self.ampl_rap_pct, pulse_shape="blackman"
        )

    def _build_arguments_ion_parameters(self):
        """
        Build arguments for ion frequencies
        """
        _argstr = "ion_parameters"

        self.setattr_argument("freq_secular_khz", NumberValue(
            default=710,
            min=500, max=3000, step=0.001,
            unit="kHz", scale=1, precision=6),
                              group=_argstr,
                              tooltip="Secular frequency (in kHz) of the ion")

        self.setattr_argument("freq_carrier_mhz_list", Scannable(
            default=[
                ExplicitScan([101.075]),
                CenterScan(101.075, 0.01, 0.0001, randomize=True),
                RangeScan(101.075, 101.1018, 50, randomize=True),
            ],
            global_min=60., global_max=400, global_step=1,
            unit="MHz", scale=1, precision=6
        ), group=_argstr,
                              tooltip="Carrier frequency of the ion.\n"
                                      "Note: this is applied via the main doublepass DDS.\n")

        self.setattr_argument('freq_cat_carrier_detuning_khz', NumberValue(default=0,
                                                              min=-10, max=10, step=1,
                                                              scale=1, precision=6),
                              group=_argstr, tooltip="Detuning frequency of the ion during resonant catting\n")

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
        self.setattr_argument("target_ms_phase",
                              EnumerationValue(['RSB', 'BSB', 'RSB-BSB', 'RSB+BSB'], default='RSB+BSB'),
                              group=_argstr,
                              tooltip="Phase update array for the singlepass DDSs during the ms gate.\n"
                                      "This configures how phase_ms_turns_list are to be applied to the DDSs.")

        # scanning options
        self.setattr_argument("time_ms_gate_us_list", Scannable(
            default=[
                ExplicitScan([50.]),
                RangeScan(0, 500, 50, randomize=True),
            ],
            global_min=1, global_max=10000, global_step=1,
            unit="us", scale=1, precision=5
        ),
                              group=_argstr,
                              tooltip="Pulse time for the the ms gate")

        self.setattr_argument("freq_ms_secular_detuning_khz_list", Scannable(
            default=[
                ExplicitScan([0]),
                CenterScan(0, 4, 0.1, randomize=True),
                RangeScan(-100, 100, 50, randomize=True),
            ],
            global_min=-100, global_max=10000, global_step=1,
            unit="kHz", scale=1, precision=3
        ), group=_argstr,
                             tooltip="Single-sided detuning from the secular frequency for the ms gate, applied via singlepass DDSs.\n"
                                     "The singlepass1 DDS is treated as the RSB, and will thus have its frequency from the secular DECREASED by this amount.\n"
                                     "Similarly, the singlepass2 DDS is treated as the BSB, and will thus have its frequency from the secular INCREASED by this amount.\n"
                                     "i.e. frequencies for [singlepass1, singlepass2] is set as [beams.freq_mhz.freq_singlepass1_mhz - freq_secular_khz - freq_ms_secular_detuning, "
                                     "beams.freq_mhz.freq_singlepass2_mhz + freq_secular_khz + freq_ms_secular_detuning].")

        self.setattr_argument("phase_ms_turns_list", Scannable(
            default=[
                ExplicitScan([0.]),
                RangeScan(0, 1.0, 11, randomize=True),
            ],
            global_min=-1.0, global_max=1.0, global_step=0.1,
            unit="turns", scale=1, precision=3
        ), group=_argstr,
                              tooltip="Phase sweep values applied to the singlepass DDSs during the ms gate.\n"
                                      "These values are multiplied/scaled by the array specified by target_cat2_cat_phase.")

        self.setattr_argument('phase_ms_dynamical_decoupling_turns_list',
                              Scannable(default=[
                                  ExplicitScan([0.]),
                                  RangeScan(0, 1, 11),
                                  CenterScan(0.5, 1, 0.1)], unit='turns',
                                  global_min=0., global_max=2., global_step=0.1, precision=3, scale=1.0),
                              tooltip='Phase of the third rf tone applied to the 729 single pass AOM during the MS gate\n'
                                      'For dynamical decoupling to work the phase must be configured correctly as we \n'
                                      'need [H_{bi}, H_{dd}]=0 and this only occurs if these Hamiltonians contain \n'
                                      'the same linear combination of c1*sigma_{x}+c2*sigma_{y} which is determined '
                                      'by \n the laser phase.',
                              group=_argstr)

    def _build_arguments_parity(self):

        _name = 'parity_pulse'

        self.setattr_argument('enable_parity_pulse', BooleanValue(False),
                              tooltip='Enable parity analysis pulse after applying the MS gate',
                              group=_name)

        self.setattr_argument('time_parity_pulse_us',
                              NumberValue(default=1.2, min=0.01, max=1e5,
                                          step=0.01, scale=1, precision=4, unit='us'),
                              tooltip='Length of the parity pulse applied after the MS gate. \n'
                                      'This time should correspond to the pi/2 carrier time.', group=_name)

        self.setattr_argument('phase_parity_pulse_turns_list',
                            Scannable(default=[
                                ExplicitScan([0.]),
                                RangeScan(0, 1, 11),
                                CenterScan(0.5, 1, 0.1)], unit='turns',
                                global_min=0., global_max=2., global_step=0.1, precision=3, scale=1.0),
                            tooltip = 'Phase of the parity pulse. \n'
                                      'For a N-ion entangled state a parity pulse the population parity \n'
                                      'should reveal a cos(N \phi) oscilaation',
                            group = _name)

    def _build_arguments_cat_default(self):
        """
        Build arguments for default cat beam parameters.
        """
        _argstr = 'cat.default'
        self.setattr_argument("time_cat_bichromatic_us",
                              NumberValue(default=50, precision=2, step=5, min=0.1, max=10000000, scale=1., unit="us"),
                              group=_argstr,
                              tooltip="Pulse time for the bichromatic pulses.")
        self.setattr_argument("ampls_cat_pct", PYONValue([50., 50.]), group=_argstr,
                              tooltip="DDS amplitudes for the singlepass DDSs during the bichromatic pulses.\n"
                                      "Should be a list of [rsb_ampl_pct, bsb_ampl_pct], which are applied to [singlepass1, singlepass2].")
        self.setattr_argument("atts_cat_db", PYONValue([11., 11.]), group=_argstr,
                              tooltip="DDS attenuations for the singlepass DDSs during the bichromatic pulses.\n"
                                      "Should be a list of [rsb_att_db, bsb_att_db], which are applied to [singlepass1, singlepass2].")

        self.setattr_argument("freq_cat_secular_detuning_khz_list", Scannable(
            default=[
                ExplicitScan([0]),
                CenterScan(0, 4, 0.1, randomize=True),
                RangeScan(100, 100, 50, randomize=True),
            ],
            global_min=-100, global_max=100, global_step=1,
            unit="kHz", scale=1, precision=3
        ), group=_argstr,
                             tooltip="Single-sided detuning frequency for the bichromatic pulses, applied via singlepass DDSs.\n"
                                     "The singlepass1 DDS is treated as the RSB, and will thus have its frequency DECREASED by this amount.\n"
                                     "Similarly, the singlepass2 DDS is treated as the BSB, and will thus have its frequency INCREASED by this amount.\n"
                                     "i.e. frequencies for [singlepass1, singlepass2] is set as [beams.freq_mhz.freq_singlepass1_mhz - freq_secular_khz + freq_cat_carrier_detuning_khz, "
                                     "beams.freq_mhz.freq_singlepass2_mhz + freq_secular_khz + freq_cat_carrier_detuning_khz].")

        self.setattr_argument('phase_cat_dynamical_decoupling_turns_list',
                              Scannable(default=[
                                  ExplicitScan([0.]),
                                  RangeScan(0, 1, 11),
                                  CenterScan(0.5, 1, 0.1)], unit='turns',
                                  global_min=0., global_max=2., global_step=0.1, precision=3, scale=1.0),
                              tooltip='Phase of the third rf tone applied to the 729 single pass AOM during the cats\n'
                                      'For dynamical decoupling to work the phase must be configured correctly as we \n'
                                      'need [H_{bi}, H_{dd}]=0 and this only occurs if these Hamiltonians contain \n'
                                      'the same linear combination of c1*sigma_{x}+c2*sigma_{y} which is determined '
                                      'by \n the laser phase.',
                              group=_argstr)

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

    def _build_arguments_cat1(self):
        """
        Build arguments for bichromatic/cat pulse #1.
        """
        # cat 1 config
        _argstr = "cat1"
        self.setattr_argument("enable_cat1_bichromatic", BooleanValue(default=True), group=_argstr,
                              tooltip="Enables application of the 1st bichromatic pulse.")
        self.setattr_argument("enable_cat1_herald", BooleanValue(default=False), group=_argstr,
                              tooltip="Enables spin-state heralding via state-selective fluorescence. "
                                      "Heralding only progresses if the state is dark, since otherwise, the motional state is destroyed.\n"
                                      "Pulses are applied as [bichromatic, herald, quench].\n"
                                      "Note: uses adaptive readout - ensure adaptive readout arguments are correctly set in the dataset manager.")
        self.setattr_argument("enable_cat1_quench", BooleanValue(default=False), group=_argstr,
                              tooltip="Enables quenching via 854nm to return the spin-state to the S-1/2 state.\n"
                                      "Note: if quench is applied to a superposition state, then the result is a mixed state, not a pure state.\n"
                                      "Pulses are applied as [bichromatic, herald, quench].")
        self.setattr_argument("phases_pulse1_cat_turns", PYONValue([0., 0.]), group=_argstr,
                              tooltip="Relative phases for the singlepass DDSs during the 1st bichromatic pulse.\n"
                                      "Should be a list of [rsb_phase_turns, bsb_phase_turns].\n"
                                      "Note: these phases are applied to the singlepass DDSs, so do not need to be halved, "
                                      "unlike phases applied to the main doublepass.")

    def _build_arguments_ramsey(self):
        _argstr = 'ramsey'
        # ramsey delay between cat1 & cat2
        self.setattr_argument("enable_ramsey_delay", BooleanValue(default=False), group=_argstr,
                              tooltip="Enables a Ramsey delay between the 1st and 2nd bichromatic pulses. "
                                      "Useful for doing motional coherence tests.")
        self.setattr_argument("time_ramsey_delay_us_list", Scannable(
            default=[
                ExplicitScan([100]),
                RangeScan(0, 500, 50, randomize=True),],
            global_min=1, global_max=100000, global_step=1,
            unit="us", scale=1, precision=5),
                              group=_argstr,
                              tooltip="Ramsey delay time between 1st and 2nd bichromatic pulses.")

    def _build_arguments_cat2(self):
        """
        Build arguments for bichromatic/cat pulse #2.
        """
        # cat2 - config
        _argstr = 'cat2'
        self.setattr_argument("enable_cat2_bichromatic", BooleanValue(default=True),
                              group=_argstr,
                              tooltip="Enables application of the 2nd bichromatic pulse.")
        self.setattr_argument("enable_cat2_herald", BooleanValue(default=False), group=_argstr,
                              tooltip="Enables spin-state heralding via state-selective fluorescence. "
                                      "Heralding only progresses if the state is dark, since otherwise, the motional state is destroyed.\n"
                                      "Pulses are applied as [bichromatic, herald, quench].\n"
                                      "Note: uses adaptive readout - ensure adaptive readout arguments are correctly set in the dataset manager.")
        self.setattr_argument("enable_cat2_quench", BooleanValue(default=False), group=_argstr,
                              tooltip="Enables quenching via 854nm to return the spin-state to the S-1/2 state.\n"
                                      "Note: if quench is applied to a superposition state, then the result is a mixed state, not a pure state.\n"
                                      "Pulses are applied as [bichromatic, herald, quench].")

        self.setattr_argument("target_cat2_cat_phase",
                              EnumerationValue(['RSB', 'BSB', 'RSB-BSB', 'RSB+BSB'], default='RSB-BSB'),
                              group=_argstr,
                              tooltip="Phase update array for the singlepass DDSs during the 2nd bichromatic pulse.\n"
                                      "This configures how phase_cat2_turns_list are to be applied to the DDSs.")
        self.setattr_argument("phase_cat2_turns_list", Scannable(
            default=[
                ExplicitScan([0.]),
                RangeScan(0, 1.0, 11, randomize=True),
            ],
            global_min=-1.0, global_max=1.0, global_step=0.1,
            unit="turns", scale=1, precision=3
        ), group=_argstr,tooltip="Phase sweep values applied to the singlepass DDSs during the 2nd bichromatic pulse.\n"
                                "These values are multiplied/scaled by the array specified by target_cat2_cat_phase.")

    def _build_arguments_readout(self):
        """
        Build arguments for readout pulse.
        """
        # RAP-based readout
        self.setattr_argument("att_rap_db",
                              NumberValue(default=8, precision=1, step=0.5, min=8, max=31.5, unit="dB", scale=1.),
                              group="read.RAP")
        self.setattr_argument("ampl_rap_pct",
                              NumberValue(default=50., precision=3, step=5, min=1, max=50, unit="%", scale=1.),
                              group="read.RAP")
        self.setattr_argument("freq_rap_center_mhz",
                              NumberValue(default=100.755, precision=6, step=1e-2, min=1, max=200, unit="MHz",
                                          scale=1.),
                              group='read.RAP')
        self.setattr_argument("freq_rap_dev_khz",
                              NumberValue(default=72., precision=2, step=0.01, min=1, max=1e4, unit="kHz", scale=1.),
                              group='read.RAP')
        self.setattr_argument("time_rap_us",
                              NumberValue(default=400., precision=3, min=1, max=1e7, step=1, unit="us", scale=1.),
                              group="read.RAP")

    def _build_arguments_tickle(self):
        """
        Build core sweep arguments for the tickle pulse.
        """
        _argstr = "tickle"  # string to use for arguments

        # waveform - parameter sweeps
        self.setattr_argument("enable_tickle_pulse", BooleanValue(default=True),
                              group=_argstr,
                              tooltip="Enables the tickle pulse.")

        self.setattr_argument("att_tickle_db",
                              NumberValue(default=25., precision=1, step=0.5, min=0., max=31.5, unit="dB", scale=1.),
                              group=_argstr,
                              tooltip="Attenuation to be used for the urukul channel used for generating the tickle.")
        self.setattr_argument("ampl_tickle_pct",
                              NumberValue(default=50., precision=2, min=0., max=50., unit="%", scale=1.),
                              group=_argstr,
                              tooltip='Amplitude of tickle pulse.')

        self.setattr_argument("time_heating_us",
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

        self.setattr_argument("freq_tickle_detuning_khz_list", Scannable(
            default=[
                ExplicitScan([0]),
                CenterScan(0., 10., 0.001, randomize=True),
                RangeScan(-10, 10, 26, randomize=True),
            ],
            global_min = -1000, global_max=1000, global_step=0.001,
            unit="kHz", scale=1, precision=6),
                              group=_argstr,
                              tooltip="Detuning from secular frequency of tickle pulse (in kHz) applied via the urukul dds.")

        self.setattr_argument("phase_tickle_turns_list", Scannable(
            default=[
                ExplicitScan([0.]),
                RangeScan(0, 1.0, 26, randomize=True),
            ],
            global_min=0.0, global_max=1.0, global_step=1,
            unit="turns", scale=1, precision=5),
                              group=_argstr,
                              tooltip="Phase of tickle pulse (in turns) applied via the urukul dds.")

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

        ### Build Pulse Shaper
        if self.enable_pulse_shaping:
            pulse_shape = self.type_pulse_shape
        else:
            pulse_shape = 'square'
        self.dds_pulse_shaper = DDSPulseShaper(self, dds_target= self.dds_dipole.dds,
                                              ram_profile=self.profile_tickle_RAM,
                                              ram_addr_start=202, num_samples=250,
                                              ampl_max_pct=self.ampl_tickle_pct,
                                               pulse_shape=pulse_shape,
                                               phase_autoclear = 1)

        ### MAGIC NUMBERS ###
        self.time_adapt_read_slack_mu = self.core.seconds_to_mu(20 * us)  # post-heralding slack to fix RTIOUnderflows
        self.time_herald_slack_mu = self.core.seconds_to_mu(150 * us)  # add slack only if herald success
        self.max_herald_attempts = 200  # max herald attempts before config skipped
        self.urukul_setup_time_mu = int64(8)  # extra delay between calls to ad9910.set_mu
        self.urukul_dd_reset_time = int64(650) # half the time for urukul to implement set_mu for DD phase shifting

        # run component preparation
        self._prepare_experiment_readout()

        self._prepare_experiment_dynamical_decoupling()

        (freq_cat_secular_detuning_ftw_list,
         phase_cat2_cat_pow_list, time_ramsey_delay_mu_list, phase_cat_dynamical_decoupling_pow_list) = self._prepare_experiment_cat_general()

        freq_carrier_ftw_list = self._prepare_experiment_ion_parameters()
        freq_tickle_detuning_ftw_list, phase_tickle_pow_list = self._prepare_experiment_tickle()

        (time_ms_gate_mu_list, freq_ms_gate_secular_detuning_khz_list, phase_ms_pow_list,
         phase_ms_dynamical_decoupling_pow_list) = self._prepare_experiment_ms_gate()

        phase_parity_pulse_pow_list = self._prepare_experiment_parity_pulse()
        self.phase_dd_phase_shift_pow = self.qubit.turns_to_pow(0.5)


        # create experiment config
        self.config_experiment_list = create_experiment_config(
            # bichromatic sweeps
            freq_carrier_ftw_list,
            # cat sweeps
            freq_cat_secular_detuning_ftw_list,
            phase_cat2_cat_pow_list,
            time_ramsey_delay_mu_list,
            # tickle sweeps
            freq_tickle_detuning_ftw_list, phase_tickle_pow_list,
            # dynamical decoupling sweeps
            phase_cat_dynamical_decoupling_pow_list, phase_ms_dynamical_decoupling_pow_list,
            # ms gate sweeps
            time_ms_gate_mu_list, freq_ms_gate_secular_detuning_khz_list, phase_ms_pow_list,
            # partity pulse sweeps
            phase_parity_pulse_pow_list,

            config_type=float, shuffle_config=True
        )

        # create profile list for later use
        self.profiles = [self.profile_729_cat1, self.profile_729_cat2, self.profile_729_ms, self.profile_729_parity]

        # placeholder arrays for urukul values (first index is  profile number, second index is channel on urukul 0)
        self.phase_beams_pow_list = zeros((8, 4), dtype=int32)
        self.freq_beams_ftw_list = zeros((8, 4), dtype=int32)
        self.ampl_beams_asf_list = zeros((8, 4), dtype=int32)

        # set timing for intensity servo between shots
        self.time_servo_relock_mu = self.core.seconds_to_mu(self.time_servo_relock_us * us)

        self.set_default_profile_configuration()

    def _prepare_experiment_ion_parameters(self):
        """
        Prepare general ion parameters
        :return: list of carrier frequencies for cat and MS
        """
        self.freq_secular_ftw = self.qubit.frequency_to_ftw(self.freq_secular_khz*kHz)

        freq_carrier_ftw_list = array([self.qubit.frequency_to_ftw(freq_mhz * MHz)
                                          for freq_mhz in self.freq_carrier_mhz_list])

        self.freq_cat_carrier_detuning_ftw = self.qubit.frequency_to_ftw(self.freq_cat_carrier_detuning_khz*kHz)

        return freq_carrier_ftw_list

    def _prepare_experiment_ms_gate(self):
        """
        Prepare experimental values for ms gate

        :return: tuple of ms gate time, detunings from secular frequency, and ms beam phases
        """

        if self.enable_ms_gate:
            time_ms_gate_mu_list = [self.core.seconds_to_mu(time_ms_gate_us/2*us) for time_ms_gate_us in self.time_ms_gate_us_list]
            freq_ms_gate_secular_detuning_khz_list = [self.qubit.frequency_to_ftw(ms_gate_secular_detuning_khz*kHz)
                                                  for ms_gate_secular_detuning_khz in self.freq_ms_secular_detuning_khz_list]
            phase_ms_pow_list = array(
                [self.qubit.turns_to_pow(phase_ms_turns) for phase_ms_turns in self.phase_ms_turns_list])
            if self.enable_dynamical_decoupling:
                    phase_ms_dynamical_decoupling_pow_list = array(
                        [self.qubit.singlepass0.turns_to_pow(phase_dynamical_decoupling_turns) for
                         phase_dynamical_decoupling_turns in self.phase_ms_dynamical_decoupling_turns_list])
            else:
                phase_ms_dynamical_decoupling_pow_list = [0]

        else:
            time_ms_gate_mu_list = [0]
            freq_ms_gate_secular_detuning_khz_list = [0]
            phase_ms_pow_list = [0]
            phase_ms_dynamical_decoupling_pow_list = [0]

        if self.enable_dynamical_decoupling:
            phase_ms_dynamical_decoupling_pow_list = array(
                [self.qubit.singlepass0.turns_to_pow(phase_dynamical_decoupling_turns) for
                 phase_dynamical_decoupling_turns in self.phase_ms_dynamical_decoupling_turns_list])

        self.ampls_ms_asf = array([self.qubit.amplitude_to_asf(ampl_ms_pct/100.) for ampl_ms_pct in self.ampls_ms_pct])

        # specify phase update array based on user arguments
        if self.target_ms_phase == 'RSB':
            self.phase_ms_update_dir = array([1, 0], dtype=int32)
        elif self.target_ms_phase == 'BSB':
            self.phase_ms_update_dir = array([0, 1], dtype=int32)
        elif self.target_ms_phase == 'RSB-BSB':
            self.phase_ms_update_dir = array([1, -1], dtype=int32)
        elif self.target_ms_phase == 'RSB+BSB':
            self.phase_ms_update_dir = array([1, 1], dtype=int32)

        self.att_reg_ms_gate = 0x00000000 | (
                (self.qubit.att_qubit_mu << ((self.qubit.beam.chip_select - 4) * 8)) |
                (self.att_dynamical_decoupling_mu << ((self.qubit.singlepass0.chip_select - 4) * 8)) |
                (att_to_mu(self.atts_ms_db[0] * dB) << ((self.qubit.singlepass1.chip_select - 4) * 8)) |
                (att_to_mu(self.atts_ms_db[1] * dB) << ((self.qubit.singlepass2.chip_select - 4) * 8))
        )
        return (time_ms_gate_mu_list, freq_ms_gate_secular_detuning_khz_list, phase_ms_pow_list,
                phase_ms_dynamical_decoupling_pow_list)

    def _prepare_experiment_parity_pulse(self):
        """
        Prepare experimental values for parity pulse
        :return: list of parity pulses phases (in pow)
        """
        self.time_parity_pulse_mu = self.core.seconds_to_mu(self.time_parity_pulse_us * us)
        if self.enable_parity_pulse:
            phase_parity_pulse_pow_list = [self.qubit.turns_to_pow(phase_parity_pulse)
                                           for phase_parity_pulse in self.phase_parity_pulse_turns_list]
        else:
            phase_parity_pulse_pow_list = [0]

        self.att_reg_parity_pulse = 0x00000000 | (
                (self.qubit.att_qubit_mu << ((self.qubit.beam.chip_select - 4) * 8)) |
                (self.qubit.att_singlepass0_default_mu << ((self.qubit.singlepass0.chip_select - 4) * 8)) |
                (att_to_mu(31.5 * dB) << ((self.qubit.singlepass1.chip_select - 4) * 8)) |
                (att_to_mu(31.5 * dB) << ((self.qubit.singlepass2.chip_select - 4) * 8))
        )

        return phase_parity_pulse_pow_list

    def _prepare_experiment_dynamical_decoupling(self):
        """
        Prepare experiment for dynamical decoupling.
        """
        if self.enable_dynamical_decoupling:
            self.ampl_dynamical_decoupling_asf = self.qubit.singlepass0.amplitude_to_asf(
                self.ampl_dynamical_decoupling_pct / 100.)
            self.att_dynamical_decoupling_mu = self.qubit.singlepass0.cpld.att_to_mu(
                self.att_dynamical_decoupling_dB * dB)
        else:
            self.att_dynamical_decoupling_mu = self.qubit.att_singlepass0_default_mu
            self.ampl_dynamical_decoupling_asf = self.qubit.ampl_singlepass0_default_asf

    def _prepare_experiment_readout(self):
        """
        Prepare experiment values for state readout.
        """
        # configure readout method
        if 'RAP' in self.readout_type:  # 2 is RAP
            self.readout_config = 1

        elif 'None' in self.readout_type:  # -1 is None
            self.readout_config = -1

        # do a final check of the readout_type argument
        if not any(kw in self.readout_type for kw in ('None', 'RAP')):
            raise ValueError("Invalid readout type. Must be one of (None, RAP).")

        # prepare RAP arguments
        self.freq_rap_center_ftw = self.qubit.frequency_to_ftw(self.freq_rap_center_mhz * MHz)
        self.freq_rap_dev_ftw = self.qubit.frequency_to_ftw(self.freq_rap_dev_khz * kHz)
        self.time_rap_mu = self.core.seconds_to_mu(self.time_rap_us * us)

    def _prepare_experiment_cat_general(self):
        """
        Prepare experiment values for cat/bichromatic.
        :return: tuple of (freq_cat_secular_detuning_ftw_list,
            phase_cat2_cat_pow_list, time_ramsey_delay_mu_list, phase_cat_dynamical_decoupling_pow_list)
        """
        '''
        CONVERT VALUES TO MACHINE UNITS - BICHROMATIC/CAT DEFAULTS
        '''
        # defaults - main doublepass (near chamber)
        freq_cat_secular_detuning_ftw_list = array([self.qubit.singlepass0.frequency_to_ftw(freq_khz * kHz)
                                           for freq_khz in self.freq_cat_secular_detuning_khz_list])
        self.ampls_cat_asf = array([self.qubit.singlepass0.amplitude_to_asf(ampl_pct / 100.)
                                    for ampl_pct in self.ampls_cat_pct])
        self.time_cat_bichromatic_mu = self.core.seconds_to_mu(self.time_cat_bichromatic_us/2 * us)

        '''
        CONVERT VALUES TO MACHINE UNITS - BICHROMATIC/CAT PULSES
        '''
        # cat1 values
        self.phases_pulse1_cat_pow = [self.qubit.singlepass0.turns_to_pow(phas_pow)
                                      for phas_pow in self.phases_pulse1_cat_turns]


        '''
        CONVERT VALUES TO MACHINE UNITS - DD PHASE
        '''


        if self.enable_dynamical_decoupling and (self.enable_cat1_bichromatic or self.enable_cat2_bichromatic
                                                or self.enable_ramsey_delay or self.enable_tickle_pulse):
            phase_cat_dynamical_decoupling_pow_list = array(
            [self.qubit.singlepass0.turns_to_pow(phase_dynamical_decoupling_turns) for
             phase_dynamical_decoupling_turns in self.phase_cat_dynamical_decoupling_turns_list])
        else:
            phase_cat_dynamical_decoupling_pow_list = [0]

        # inter-cat ramsey delay
        if self.enable_ramsey_delay:
            time_ramsey_delay_mu_list = [self.core.seconds_to_mu(time_delay_us/2 * us)
                                         for time_delay_us in self.time_ramsey_delay_us_list]
            time_ramsey_dd_phase_flip_mu_list = [self.core.seconds_to_mu(time_delay_us/2 * us)
                                         for time_delay_us in self.time_ramsey_delay_us_list]
        else:
            time_ramsey_delay_mu_list = [0]
            time_ramsey_dd_phase_flip_mu_list = [0]

        # cat2 values
        if self.enable_cat2_bichromatic:
            phase_cat2_cat_pow_list = [self.qubit.singlepass0.turns_to_pow(phas_pow)
                                       for phas_pow in self.phase_cat2_turns_list]
        else:
            phase_cat2_cat_pow_list = [0]

        # specify phase update array based on user arguments
        if self.target_cat2_cat_phase == 'RSB':
            self.phase_cat_update_dir = array([1, 0], dtype=int32)
        elif self.target_cat2_cat_phase == 'BSB':
            self.phase_cat_update_dir = array([0, 1], dtype=int32)
        elif self.target_cat2_cat_phase == 'RSB-BSB':
            self.phase_cat_update_dir = array([1, -1], dtype=int32)
        elif self.target_cat2_cat_phase == 'RSB+BSB':
            self.phase_cat_update_dir = array([1, 1], dtype=int32)

        '''
        CREATE ATTENUATION REGISTERS
        '''
        # attenuation register - bichromatic: main doublepass set to specified experiment argument value
        self.att_reg_cat_interferometer = 0x00000000 | (
                (self.qubit.att_qubit_mu << ((self.qubit.beam.chip_select - 4) * 8)) |
                (self.att_dynamical_decoupling_mu << ((self.qubit.singlepass0.chip_select - 4) * 8)) |
                (att_to_mu(self.atts_cat_db[0] * dB) << ((self.qubit.singlepass1.chip_select - 4) * 8)) |
                (att_to_mu(self.atts_cat_db[1] * dB) << ((self.qubit.singlepass2.chip_select - 4) * 8))
        )

        # attenuation register - readout (RAP): singlepasses set to default
        self.att_reg_readout_rap = 0x00000000 | (
                (att_to_mu(self.att_rap_db * dB) << ((self.qubit.beam.chip_select - 4) * 8)) |
                (self.qubit.att_singlepass0_default_mu << ((self.qubit.singlepass0.chip_select - 4) * 8)) |
                (self.qubit.att_singlepass1_default_mu << ((self.qubit.singlepass1.chip_select - 4) * 8)) |
                (self.qubit.att_singlepass2_default_mu << ((self.qubit.singlepass2.chip_select - 4) * 8))
        )

        # return sweep lists required for
        return (freq_cat_secular_detuning_ftw_list,
                phase_cat2_cat_pow_list, time_ramsey_delay_mu_list, phase_cat_dynamical_decoupling_pow_list)

    def _prepare_experiment_tickle(self):
        """
        Prepare general experiment values for the tickle pulse.
        :return: tuple of (freq_tickle_detuning_hz_list, phase_tickle_list)
        """
        # convert values to convenience units
        self.att_tickle_mu = att_to_mu(self.att_tickle_db * dB)
        freq_tickle_detuning_ftw_list = [self.dds_pulse_shaper.dds_target.frequency_to_ftw(freq_tickle_detuning_khz*kHz)
                                         for freq_tickle_detuning_khz in self.freq_tickle_detuning_khz_list]
        phase_tickle_pow_list = [self.dds_pulse_shaper.dds_target.turns_to_pow(phase_tickle_turns) for phase_tickle_turns in self.phase_tickle_turns_list]
        self.time_tickle_mu = self.core.seconds_to_mu(self.time_heating_us * us)
        self.time_tickle_dd_phase_flip_mu = self.core.seconds_to_mu(self.time_heating_us * us/2)


        # don't apply sweep if tickle is disabled
        if self.enable_tickle_pulse:
            return freq_tickle_detuning_ftw_list, phase_tickle_pow_list
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
                13)

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

        # configure RAP pulse
        if self.readout_config == 1:
            self.rap_subsequence.configure(self.time_rap_mu, self.freq_rap_center_ftw, self.freq_rap_dev_ftw)
            delay_mu(50000)

        # configure tickle
        self.dds_pulse_shaper.sequence_initialize()
        self.dds_pulse_shaper.dds_target.set_att_mu(self.att_tickle_mu)
        self.dds_pulse_shaper.dds_target.sw.off()
        delay_mu(8)

    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        # # predeclare variables ahead of time
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
                freq_carrier_ftw = int32(config_vals[0])
                freq_cat_secular_detuning_ftw = int32(config_vals[1])
                phase_cat2_cat_pow = int32(config_vals[2])
                time_ramsey_delay_mu = int64(config_vals[3])
                freq_tickle_detuning_ftw = int32(config_vals[4])
                phase_tickle_pow = int32(config_vals[5])
                phase_cat_dynamical_decoupling_pow = int32(config_vals[6])
                phase_ms_dynamical_decoupling_pow = int32(config_vals[7])
                time_ms_gate_mu = int64(config_vals[8])
                freq_ms_gate_secular_detuning_ftw = int32(config_vals[9])
                phase_ms_pow = int32(config_vals[10])
                phase_parity_pow = int32(config_vals[11])

                herald_counter = 0  # clear herald counter

                self.update_profile_configuration(freq_carrier_ftw=freq_carrier_ftw,
                                                freq_cat_secular_detuning_ftw=freq_cat_secular_detuning_ftw,
                                                freq_ms_gate_secular_detuning_ftw=freq_ms_gate_secular_detuning_ftw,
                                                phase_cat2_cat_pow=phase_cat2_cat_pow,
                                                phase_cat_dynamical_decoupling_pow=phase_cat_dynamical_decoupling_pow,
                                                phase_ms_dynamical_decoupling_pow=phase_ms_dynamical_decoupling_pow,
                                                phase_ms_pow=phase_ms_pow,
                                                phase_parity_pulse_pow=phase_parity_pow)

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
                    self.setup_beam_profiles()
                    self.qubit.set_cfr1(phase_autoclear=1)
                    self.qubit.singlepass0.set_cfr1(phase_autoclear=1)
                    self.qubit.singlepass1.set_cfr1(phase_autoclear=1)
                    self.qubit.singlepass2.set_cfr1(phase_autoclear=1)

                    # set tickle frequency/phases
                    self.dds_pulse_shaper.dds_target.set_ftw(self.freq_secular_ftw + freq_tickle_detuning_ftw)
                    self.dds_pulse_shaper.dds_target.set_pow(phase_tickle_pow)

                    # set up config of shaped pulses to be fired for tickling
                    # also sets up phase autoclear
                    self.dds_pulse_shaper.configure_train(self.time_tickle_mu)
                    with parallel:
                        self.dds_pulse_shaper.dds_target.cpld.io_update.pulse_mu(8)
                        self.qubit.io_update()
                    # for ururuk channel used for tickling keep RAM enabled but ensure we don't clear phase on io_update
                    self.dds_pulse_shaper.dds_target.set_cfr1(ram_enable=1, phase_autoclear=0,
                                                              ram_destination=ad9910.RAM_DEST_ASF)
                    self.dds_pulse_shaper.dds_target.cpld.io_update.pulse_mu(8)

                    # reset cfr1 so we no longer clear phases on io_update
                    self.qubit.set_cfr1()
                    self.qubit.singlepass0.set_cfr1()
                    self.qubit.singlepass1.set_cfr1()
                    self.qubit.singlepass2.set_cfr1()
                    self.qubit.io_update()

                    '''
                    MS Gate
                    '''
                    if self.enable_ms_gate:
                        self.qubit.cpld.set_all_att_mu(self.att_reg_ms_gate)
                        self.pulse_bichromatic(self.profile_729_ms,
                                               time_ms_gate_mu,
                                               phase_ms_dynamical_decoupling_pow)

                    '''
                    Parity Pulse
                    '''
                    if self.enable_parity_pulse:
                        self.pulse_parity()

                    '''
                    CAT #1
                    '''
                    self.qubit.cpld.set_all_att_mu(self.att_reg_cat_interferometer)
                    # cat1 - bichromatic cat pulse
                    if self.enable_cat1_bichromatic:
                        self.pulse_bichromatic(self.profile_729_cat1,
                                       self.time_cat_bichromatic_mu,
                                       phase_cat_dynamical_decoupling_pow)

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
                        # ensure we do not have a negative ramsey time when DDing
                        # with DD minimuim delay time is ~1.3us, i.e time it takes to switch DD phase
                        if time_ramsey_delay_mu - self.urukul_dd_reset_time < 8:
                            time_ramsey_delay_mu = 8
                        else:
                            time_ramsey_delay_mu = time_ramsey_delay_mu - self.urukul_dd_reset_time
                        if self.enable_dynamical_decoupling:
                            self.set_carrier_phase(phase_cat_dynamical_decoupling_pow - self.phase_dd_phase_shift_pow,
                                                   profile=self.profile_729_cat1)
                            self.qubit.singlepass0_on()
                            self.qubit.on()
                        delay_mu(time_ramsey_delay_mu)
                        if self.enable_dynamical_decoupling:
                            self.qubit.singlepass0_off()
                            self.set_carrier_phase(phase_cat_dynamical_decoupling_pow,
                                                   profile=self.profile_729_cat1)
                            self.qubit.singlepass0_on()
                        delay_mu(time_ramsey_delay_mu)
                        self.qubit.off()

                    '''
                    TICKLE PULSE
                    '''
                    if self.enable_tickle_pulse:
                        # ensure we do not have a negative tickle time when DDing
                        # with DD minimuim tickle time is ~1.3us, i.e time it takes to switch DD phase
                        if self.time_tickle_dd_phase_flip_mu - self.urukul_dd_reset_time < 8:
                            time_dd_phase_flip_mu = 8
                        else:
                            time_dd_phase_flip_mu = self.time_tickle_dd_phase_flip_mu - self.urukul_dd_reset_time
                        if self.enable_dynamical_decoupling:
                            self.set_carrier_phase(phase_cat_dynamical_decoupling_pow - self.phase_dd_phase_shift_pow,
                                                   profile=self.profile_729_cat1)
                            self.qubit.singlepass0_on()
                            self.qubit.on()
                        with parallel:
                            self.dds_pulse_shaper.run_train_single()
                            if self.enable_dynamical_decoupling:
                                # see time to set profile
                                delay_mu(time_dd_phase_flip_mu)
                                self.qubit.singlepass0_off()
                                self.set_carrier_phase(phase_cat_dynamical_decoupling_pow,
                                                       profile=self.profile_729_cat1)
                                self.qubit.singlepass0_on()
                        self.qubit.off()

                    '''
                    CAT #2
                    '''
                    if self.enable_cat2_bichromatic:
                        self.pulse_bichromatic(self.profile_729_cat2,
                                           self.time_cat_bichromatic_mu,
                                            phase_cat_dynamical_decoupling_pow)

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
                    # RAP based readout
                    if self.readout_config == 1:
                        self.pulse_readout_rap()

                    # read out fluorescence & clean up loop
                    self.readout_subsequence.run_dma()
                    counts_res = self.readout_subsequence.fetch_count()
                else:
                    # return -1 so user knows booboo happened
                    counts_res = -1

                # cleanup dds_pulse_shaper
                self.dds_pulse_shaper.sequence_cleanup()

                # store results
                self.update_results(freq_carrier_ftw,
                                    counts_res,
                                    freq_cat_secular_detuning_ftw,
                                    phase_cat2_cat_pow,
                                    time_ramsey_delay_mu << 1, # shift bits over to left to account for halving we did for DD
                                    freq_tickle_detuning_ftw,
                                    phase_tickle_pow,
                                    phase_cat_dynamical_decoupling_pow,
                                    phase_ms_dynamical_decoupling_pow,
                                    time_ms_gate_mu << 1, # shift bits over to left to account for halving we did for DD
                                    freq_ms_gate_secular_detuning_ftw,
                                    phase_ms_pow,
                                    phase_parity_pow)


                # check termination more frequently in case reps are low
                if _loop_iter % 100 == 0:
                    self.check_termination()
                _loop_iter += 1

            # rescue ion as needed & support graceful termination
            self.check_termination()


    @kernel(flags={'fast-math'})
    def pulse_bichromatic(self, profile: TInt32,
                            time_pulse_mu: TInt64,
                            phase_dd_pow: TInt32=0) -> TNone:
        """
        Bichromatic interaction to produce a cat state
        :param profile: urukul profile with the proper settings
        :param time_pulse_mu: how long to apply the bichromatic interaction (which creates the cat state) for
        :param phase_dd_pow: phase (in pow) of dynamical decoupling pulse
        """
        # set everything to correct profile
        self.qubit.off()
        self.qubit.set_profile(profile)
        at_mu((now_mu() + 8) & ~7)
        self.qubit.io_update()
        delay_mu(2000)
        # turn on all beams
        self.qubit.singlepass1_on()
        self.qubit.singlepass2_on()
        if self.enable_dynamical_decoupling:
            self.set_carrier_phase(phase_dd_pow, profile=profile)
            self.qubit.singlepass0_on()
            # correct time to account (and ensure it is non-negative) for call to set_mu to change dd phase
            # minimium CAT time with DD is ~1.3us, i.e time for DD phase shift
            if time_pulse_mu - self.urukul_dd_reset_time < 8:
                time_pulse_mu = 8
            else:
                time_pulse_mu = time_pulse_mu - self.urukul_dd_reset_time
        else:
            self.qubit.singlepass0_off()

        # turn on main doublepass and begin bichromatic
        self.qubit.on()
        delay_mu(time_pulse_mu)
        if self.enable_dynamical_decoupling:
            self.qubit.singlepass0_off()
            self.set_carrier_phase(phase_dd_pow - self.phase_dd_phase_shift_pow, profile=profile)
            self.qubit.singlepass0_on()
        delay_mu(time_pulse_mu)

        # turn off all beams except singlepass 0 to prevent thermal fluctuations on the singlepass
        self.qubit.off()
        self.qubit.singlepass1_off()
        self.qubit.singlepass2_off()
        self.qubit.singlepass0_on()

    @kernel(flags={'fast-math'})
    def pulse_parity(self) -> TNone:
        self.qubit.cpld.set_all_att_mu(self.att_reg_parity_pulse)
        self.qubit.set_profile(self.profile_729_parity)
        at_mu((now_mu() + 8) & ~7)
        self.qubit.io_update()
        delay_mu(2000)

        self.qubit.on()
        self.qubit.singlepass0_on()
        delay_mu(self.time_parity_pulse_mu)
        self.qubit.off()
        self.qubit.singlepass0_off()

    @kernel(flags={"fast-math"})
    def set_carrier_phase(self, phase_dd_pow_: TInt32,
                             profile: TInt32) -> TNone:
        """
        Change phase of singlepass0 for spin refocusing during DD
        :param phase_dd_pow_: phase of DD pulse used
        :param profile: urukul profile with the proper settings
        """
        self.qubit.singlepass0.set_mu(
            self.qubit.freq_singlepass0_default_ftw,
            asf=self.ampl_dynamical_decoupling_asf,
            pow_=phase_dd_pow_,
            profile=profile,
            phase_mode=ad9910.PHASE_MODE_CONTINUOUS
        )


    @kernel(flags={"fast-math"})
    def setup_beam_profiles(self) -> TNone:
        """
        Configure parameters for relevant profiles on urukul
        """
        # ensure all beams are off
        self.qubit.off()
        self.qubit.singlepass0_off()
        self.qubit.singlepass1_off()
        self.qubit.singlepass2_off()
        # set up relevant beam waveforms
        for profile in self.profiles:
            self.qubit.set_mu(
                self.freq_beams_ftw_list[profile][0],
                asf=self.ampl_beams_asf_list[profile][0],
                pow_=self.phase_beams_pow_list[profile][0],
                profile=profile,
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
        """
        Store the profile values common across all profiles for all shots
        """
        for profile in self.profiles:
            self.ampl_beams_asf_list[profile][0] = self.qubit.ampl_qubit_asf

            self.phase_beams_pow_list[profile][0] = 0

            self.freq_beams_ftw_list[profile][1] = self.qubit.freq_singlepass0_default_ftw

        # set up values consistent with bichromatic operations
        for profile in [self.profile_729_cat1, self.profile_729_cat2, self.profile_729_ms]:
            self.ampl_beams_asf_list[profile][1] = self.ampl_dynamical_decoupling_asf


        # set up values consistent across cat profiles
        for profile in [self.profile_729_cat1, self.profile_729_cat2]:
            self.ampl_beams_asf_list[profile][2] = self.ampls_cat_asf[0]
            self.ampl_beams_asf_list[profile][3] = self.ampls_cat_asf[1]

        self.phase_beams_pow_list[self.profile_729_cat1][2] = self.phases_pulse1_cat_pow[0]
        self.phase_beams_pow_list[self.profile_729_cat1][3] = self.phases_pulse1_cat_pow[1]

        # ms gate parameters
        self.ampl_beams_asf_list[self.profile_729_ms][2] = self.ampls_ms_asf[0]
        self.ampl_beams_asf_list[self.profile_729_ms][3] = self.ampls_ms_asf[1]

        # parity pulse parameters
        self.ampl_beams_asf_list[self.profile_729_parity][1] = self.qubit.ampl_singlepass0_default_asf


    @kernel(flags={'fast-math'})
    def update_profile_configuration(self, freq_carrier_ftw: TInt32, freq_cat_secular_detuning_ftw: TInt32,
                                     freq_ms_gate_secular_detuning_ftw: TInt32,
                                     phase_cat2_cat_pow: TInt32,
                                     phase_cat_dynamical_decoupling_pow: TInt32,
                                     phase_ms_dynamical_decoupling_pow: TInt32,
                                     phase_ms_pow: TInt32,
                                     phase_parity_pulse_pow: TInt32) -> TNone:

        """
        Fill out list for easy access of values when setting up the profiles on the urukul
        :param freq_carrier_ftw: center frequency of the red and blue sidebands
        :param freq_cat_secular_detuning_ftw: single-sided detuning of the red and blue sidebands from the center cat frequency
        :param freq_ms_gate_secular_detuning_ftw: detuning of ms gate tones from sidebands
        :param phase_cat2_cat_pow: phase of the second cat pulse
        :param phase_cat_dynamical_decoupling_pow: phase of the dynamical decoupling tone for catting operations
        :param phase_ms_dynamical_decoupling_pow: phase of dynamical decoupling tone for ms gate
        :param phase_ms_pow: phase of the ms pulse
        :param phase_parity_pulse_pow: phase of the parity pulse
        """

        # prepare phase arrays for bichromatic
        cat4_phases = [
            self.phase_cat_update_dir[0] * phase_cat2_cat_pow,
            self.phase_cat_update_dir[1] * phase_cat2_cat_pow,
        ]

        ms_phases = [
            self.phase_ms_update_dir[0] * phase_ms_pow,
            self.phase_ms_update_dir[1] * phase_ms_pow,
        ]

        for profile in self.profiles:
            self.freq_beams_ftw_list[profile][0] = freq_carrier_ftw

        # set up values consistent across cats
        for profile in [self.profile_729_cat1, self.profile_729_cat2]:
            self.freq_beams_ftw_list[profile][2] = self.qubit.freq_singlepass1_default_ftw - self.freq_secular_ftw - freq_cat_secular_detuning_ftw + self.freq_cat_carrier_detuning_ftw
            self.freq_beams_ftw_list[profile][3] = self.qubit.freq_singlepass2_default_ftw + self.freq_secular_ftw + freq_cat_secular_detuning_ftw + self.freq_cat_carrier_detuning_ftw


        # set up values for ms gate
        self.freq_beams_ftw_list[self.profile_729_ms][2] = self.qubit.freq_singlepass1_default_ftw - self.freq_secular_ftw - freq_ms_gate_secular_detuning_ftw
        self.freq_beams_ftw_list[self.profile_729_ms][3] = self.qubit.freq_singlepass2_default_ftw + self.freq_secular_ftw + freq_ms_gate_secular_detuning_ftw

        self.phase_beams_pow_list[self.profile_729_ms][1] = phase_ms_dynamical_decoupling_pow
        self.phase_beams_pow_list[self.profile_729_ms][2] = ms_phases[0]
        self.phase_beams_pow_list[self.profile_729_ms][3] = ms_phases[1]

        # set up values for parity pulse
        self.phase_beams_pow_list[self.profile_729_parity][1] = phase_parity_pulse_pow

        # set up values for cat
        self.phase_beams_pow_list[self.profile_729_cat1][1] = phase_cat_dynamical_decoupling_pow

        self.phase_beams_pow_list[self.profile_729_cat2][1] = phase_cat_dynamical_decoupling_pow
        self.phase_beams_pow_list[self.profile_729_cat2][2] = cat4_phases[0]
        self.phase_beams_pow_list[self.profile_729_cat2][3] = cat4_phases[1]

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
        # run RAP turns on qubit
        self.rap_subsequence.run_rap(self.time_rap_mu)
