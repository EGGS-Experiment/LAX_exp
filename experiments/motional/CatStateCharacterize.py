from artiq.experiment import *
from artiq.coredevice import ad9910
from numpy import array, int32, int64

from LAX_exp.language import *
from LAX_exp.system.subsequences import (
    InitializeQubit, SidebandCoolContinuousRAM, Readout, ReadoutAdaptive, RescueIon
)
# todo: migrate readout to adaptive for speed lol


class CatStateCharacterize(LAXExperiment, Experiment):
    """
    Experiment: Cat State Characterize

    Create and characterize cat states with projective state preparation.
    Uses adaptive readout to reduce timing overheads and extend available coherence times.
    """
    name = 'Cat State Characterize'
    kernel_invariants = {
        # subsequences
        'initialize_subsequence', 'sidebandcool_subsequence', 'readout_subsequence',
        'readout_adaptive_subsequence', 'rescue_subsequence',

        # hardware values - default
        'ampl_doublepass_default_asf',
        'freq_sigmax_ftw', 'ampl_sigmax_asf', 'time_sigmax_mu',
        'time_herald_slack_mu', 'time_adapt_read_slack_mu', 'max_herald_attempts',

        # hardware values - cat & readout
        'ampls_cat_asf', 'time_cat1_bichromatic_mu', 'phases_pulse1_cat_pow', 'phase_cat1_antisigmax_pow',
        'phase_cat2_sigmax_pow', 'phase_cat2_antisigmax_pow', 'phases_cat2_cat_pow', 'phases_cat2_cat_update_dir',
        'ampl_729_readout_asf', 'att_reg_sigmax', 'att_reg_bichromatic', 'att_reg_readout',

        # configs
        'profile_729_SBC', 'profile_729_target', 'config_experiment_list',
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
        self.readout_adaptive_subsequence = ReadoutAdaptive(self, time_bin_us=10, error_threshold=1e-2)
        self.rescue_subsequence =       RescueIon(self)

        # relevant devices
        self.setattr_device('qubit')
        self.setattr_device('pump')
        self.setattr_device('repump_qubit')

        # set build arguments
        self._build_arguments_default()
        self._build_arguments_cat1()
        self._build_arguments_cat2()
        self._build_arguments_readout()

    def _build_arguments_default(self):
        """
        Build arguments for default beam parameters.
        """
        # defaults - sigma_x
        self.setattr_argument("freq_sigmax_mhz",    NumberValue(default=101.0978, precision=6, step=1, min=50., max=400., scale=1., unit="MHz"),
                              group="default.sigmax",
                              tooltip="Frequency for both the sigma_x and anti-sigma_x pulses. "
                                      "Applied to main/chamber doublepass.")
        self.setattr_argument("ampl_sigmax_pct",    NumberValue(default=50., precision=3, step=5, min=0.01, max=50, scale=1., unit="%"),
                              group="default.sigmax",
                              tooltip="DDS amplitude for both the sigma_x and anti-sigma_x pulses. "
                                      "Applied to main/chamber doublepass.")
        self.setattr_argument("att_sigmax_db",      NumberValue(default=8., precision=1, step=0.5, min=8., max=31.5, scale=1., unit="dB"),
                              group="default.sigmax",
                              tooltip="DDS attenuation for both the sigma_x and anti-sigma_x pulses. "
                                      "Applied to main/chamber doublepass.")
        self.setattr_argument("time_sigmax_us",     NumberValue(default=2.24, precision=3, step=0.1, min=0.01, max=10000, scale=1., unit="us"),
                              group="default.sigmax",
                              tooltip="Pulse time for both the sigma_x and anti-sigma_x pulses. "
                                      "Applied to main/chamber doublepass.")

        # defaults - beam values - doublepass (main)
        self.setattr_argument("ampl_doublepass_default_pct",    NumberValue(default=50., precision=3, step=5, min=0.01, max=50, scale=1., unit="%"),
                              group="default.cat",
                              tooltip="DDS amplitude for the main doublepass during the bichromatic pulse.")
        self.setattr_argument("att_doublepass_default_db",      NumberValue(default=8., precision=1, step=0.5, min=8., max=31.5, scale=1., unit="dB"),
                              group="default.cat",
                              tooltip="DDS attenuation for the main doublepass during the bichromatic pulse.")
        self.setattr_argument("ampls_cat_pct",  PYONValue([50., 50.]), group='default.cat',
                              tooltip="DDS amplitudes for the singlepass DDSs during the bichromatic pulses.\n"
                                      "Should be a list of [rsb_ampl_pct, bsb_ampl_pct], which are applied to [singlepass0, singlepass1].")
        self.setattr_argument("atts_cat_db",    PYONValue([13., 13.]), group='default.cat',
                              tooltip="DDS attenuations for the singlepass DDSs during the bichromatic pulses.\n"
                                      "Should be a list of [rsb_att_db, bsb_att_db], which are applied to [singlepass0, singlepass1].")
        self.setattr_argument("freq_cat_center_mhz_list",   Scannable(
                                                                default=[
                                                                    ExplicitScan([101.0918]),
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
        self.setattr_argument("freq_cat_secular_khz_list",  Scannable(
                                                                default=[
                                                                    ExplicitScan([702.01]),
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

    def _build_arguments_cat1(self):
        """
        Build arguments for bichromatic/cat pulse #1.
        """
        # cat #1 config (sigma_x only)
        self.setattr_argument("enable_cat1_sigmax",     BooleanValue(default=True), group='cat1.config',
                              tooltip="Applies a sigma_x pulse BEFORE the 1st bichromatic pulse.\n"
                                      "If sigma_x is applied (i.e. True), the bichromatic pulse creates a pure eigenstate (e.g. coherent state).\n"
                                      "If sigma_x is disabled (i.e. False), the bichromatic pulse creates a superposition state (e.g. cat state).")
        self.setattr_argument("enable_cat1_antisigmax", BooleanValue(default=False), group='cat1.config',
                              tooltip="Applies a sigma_x pulse AFTER the 1st bichromatic pulse.\n"
                                      "If the sigma_x pulse is APPLIED, then this pulse disentangles spin from motion.\n"
                                      "If the sigma_x pulse is DISABLED, then this pulse simply selects whether an "
                                      "odd or even superposition is associated with the dark state.")
        self.setattr_argument("phase_cat1_antisigmax_turns",    NumberValue(default=0.25, precision=3, step=0.1, min=-1.0, max=1.0, scale=1., unit="turns"),
                              group='cat1.config',
                              tooltip="Relative phase applied for the anti-sigma_x pulse.\n"
                                      "Note: this phase is applied via the main doublepass DDS, so values should be halved.")

        # cat #1 config
        self.setattr_argument("enable_cat1_bichromatic",    BooleanValue(default=False), group='cat1.config',
                              tooltip="Enables application of the 1st bichromatic pulse.\n"
                                      "Pulses are applied as [sigma_x, bichromatic, antisigma_x, herald, quench].")
        self.setattr_argument("enable_cat1_herald",         BooleanValue(default=False), group='cat1.config',
                              tooltip="Enables spin-state heralding via state-selective fluorescence. "
                                      "Heralding only progresses if the state is dark, since otherwise, the motional state is destroyed.\n"
                                      "Pulses are applied as [sigma_x, bichromatic, antisigma_x, herald, quench].\n"
                                      "Note: uses adaptive readout - ensure adaptive readout arguments are correctly set in the dataset manager.")
        self.setattr_argument("enable_cat1_quench",         BooleanValue(default=False), group='cat1.config',
                              tooltip="Enables quenching via 854nm to return the spin-state to the S-1/2 state.\n"
                                      "Note: if quench is applied to a superposition state, then the result is a mixed state, not a pure state.\n"
                                      "Pulses are applied as [sigma_x, bichromatic, antisigma_x, herald, quench].")
        self.setattr_argument("time_cat1_bichromatic_us",   NumberValue(default=100, precision=2, step=5, min=0.1, max=10000, scale=1., unit="us"),
                              group="cat1.config",
                              tooltip="Pulse time for the 1st bichromatic pulse.")
        self.setattr_argument("phases_pulse1_cat_turns",    PYONValue([0., 0.]), group='cat1.config',
                              tooltip="Relative phases for the singlepass DDSs during the 1st bichromatic pulse.\n"
                                      "Should be a list of [rsb_phase_turns, bsb_phase_turns].\n"
                                      "Note: these phases are applied to the singlepass DDSs, so do not need to be halved, "
                                      "unlike phases applied to the main doublepass.")

        # ramsey delay between cat1 & cat2
        self.setattr_argument("enable_ramsey_delay",        BooleanValue(default=False), group='cat1.ramsey',
                              tooltip="Enables a Ramsey delay between the 1st and 2nd bichromatic pulses. "
                                      "Useful for doing motional coherence tests.")
        self.setattr_argument("time_ramsey_delay_us_list",  Scannable(
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
        self.setattr_argument("enable_cat2_sigmax",       BooleanValue(default=False),
                              group='cat2.config',
                              tooltip="Applies a sigma_x pulse BEFORE the 2nd bichromatic pulse.\n"
                                      "If sigma_x is applied (i.e. True), the bichromatic pulse creates a pure eigenstate (e.g. coherent state).\n"
                                      "If sigma_x is disabled (i.e. False), the bichromatic pulse creates a superposition state (e.g. cat state).")
        self.setattr_argument("phase_cat2_sigmax_turns",  NumberValue(default=0., precision=3, step=0.1, min=-1.0, max=1.0, scale=1., unit="turns"),
                              group='cat2.config',
                              tooltip="Relative phase applied for the sigma_x pulse.\n"
                                      "Note: this phase is applied via the main doublepass DDS, so values should be halved.")
        self.setattr_argument("enable_cat2_antisigmax", BooleanValue(default=False), group='cat2.config',
                              tooltip="Applies a sigma_x pulse AFTER the 2nd bichromatic pulse.\n"
                                      "If the sigma_x pulse is APPLIED, then this pulse disentangles spin from motion.\n"
                                      "If the sigma_x pulse is DISABLED, then this pulse simply selects whether an "
                                      "odd or even superposition is associated with the dark state.")
        self.setattr_argument("phase_cat2_antisigmax_turns",    NumberValue(default=0., precision=3, step=0.1, min=-1.0, max=1.0, scale=1., unit="turns"),
                              group='cat2.config',
                              tooltip="Relative phase applied for the anti-sigma_x pulse.\n"
                                      "Note: this phase is applied via the main doublepass DDS, so values should be halved.")

        # cat2 - config
        self.setattr_argument("enable_cat2_bichromatic",  BooleanValue(default=False),
                              group='cat2.config',
                              tooltip="Enables application of the 2nd bichromatic pulse.\n"
                                      "Pulses are applied as [sigma_x, bichromatic, antisigma_x, herald, quench].")
        self.setattr_argument("enable_cat2_herald",   BooleanValue(default=False), group='cat2.config',
                              tooltip="Enables spin-state heralding via state-selective fluorescence. "
                                      "Heralding only progresses if the state is dark, since otherwise, the motional state is destroyed.\n"
                                      "Pulses are applied as [sigma_x, bichromatic, antisigma_x, herald, quench].\n"
                                      "Note: uses adaptive readout - ensure adaptive readout arguments are correctly set in the dataset manager.")
        self.setattr_argument("enable_cat2_quench",   BooleanValue(default=False), group='cat2.config',
                              tooltip="Enables quenching via 854nm to return the spin-state to the S-1/2 state.\n"
                                      "Note: if quench is applied to a superposition state, then the result is a mixed state, not a pure state.\n"
                                      "Pulses are applied as [sigma_x, bichromatic, antisigma_x, herald, quench].\n")

        # cat2 - pulse parameters
        self.setattr_argument("time_cat2_cat_us_list",    Scannable(
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
        self.setattr_argument("phases_cat2_cat_turns",      PYONValue([0., 0.]), group='cat2.bichromatic',
                              tooltip="Relative phase OFFSET for the singlepass DDSs during the 2nd bichromatic pulse. "
                                      "These phases are applied IN ADDITION to the values from phase_cat2_cat_turns_list.\n"
                                      "Should be a list of [rsb_phase_turns, bsb_phase_turns].\n"
                                      "Note: these phases are applied to the singlepass DDSs, so do not need to be halved, "
                                      "unlike phases applied to the main doublepass.")
        self.setattr_argument("target_cat2_cat_phase",      EnumerationValue(['RSB', 'BSB', 'RSB-BSB', 'RSB+BSB'], default='RSB-BSB'),
                              group="cat2.bichromatic",
                              tooltip="Phase update array for the singlepass DDSs during the 2nd bichromatic pulse.\n"
                                      "This configures how phase_cat2_cat_turns_list are to be applied to the DDSs.")
        self.setattr_argument("phase_cat2_cat_turns_list",  Scannable(
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
        self.setattr_argument("enable_729_readout",     BooleanValue(default=True), group="readout_729",
                              tooltip="Enables 729nm-based readout (e.g. sideband ratio, BSB rabi).\n"
                                      "If this is disabled, then NO 729nm pulses are applied before state-selective readout.\n"
                                      "Note: readout pulses are NOT phase coherent with any bichromatic/sigma_x pulses.")
        self.setattr_argument("ampl_729_readout_pct",   NumberValue(default=50., precision=3, step=5, min=0.01, max=50, unit="%", scale=1.),
                              group="readout_729",
                              tooltip="729nm DDS amplitude (in percent of full scale) to use for readout.\n"
                                      "This is applied via the main doublepass.")
        self.setattr_argument("att_729_readout_db", NumberValue(default=8., precision=1, step=0.5, min=8., max=31.5, unit="dB", scale=1.),
                              group="readout_729",
                              tooltip="729nm DDS attenuation (in dB) to use for readout. "
                                      "This is applied via the main doublepass.")
        self.setattr_argument("freq_729_readout_mhz_list",   Scannable(
                                                                default=[
                                                                    ExplicitScan([101.4308]),
                                                                    CenterScan(101.9851, 0.01, 0.0002, randomize=True),
                                                                    RangeScan(101.9801, 101.9901, 50, randomize=True),
                                                                ],
                                                                global_min=60., global_max=400, global_step=1,
                                                                unit="MHz", scale=1, precision=6
                                                            ), group="readout_729",
                              tooltip="729nm DDS frequencies to use for readout. "
                                      "This is applied via the main doublepass.")
        self.setattr_argument("time_729_readout_us_list",    Scannable(
                                                                default=[
                                                                    # ExplicitScan([122.9]),
                                                                    RangeScan(0, 200, 200, randomize=True),
                                                                ],
                                                                global_min=1, global_max=100000, global_step=1,
                                                                unit="us", scale=1, precision=5
                                                            ), group="readout_729",
                              tooltip="729nm readout pulse times.")

    def prepare_experiment(self):
        """
        Prepare & precompute experimental values.
        """
        '''
        GENERAL SETUP
        '''
        self._prepare_argument_checks()

        ### MAGIC NUMBERS ###
        # extra slack after heralding to prevent RTIOUnderflow errors
        self.time_adapt_read_slack_mu = self.core.seconds_to_mu(20 * us) # always add slack immediately after adaptive readout
        self.time_herald_slack_mu = self.core.seconds_to_mu(150 * us)   # add slack to RTIOCounter only if heralding succeeds
        self.max_herald_attempts =  200 # max number of herald attempts before config is skipped


        '''
        CONVERT VALUES TO MACHINE UNITS - DEFAULTS
        '''
        # defaults - main doublepass (near chamber)
        self.ampl_doublepass_default_asf = self.qubit.amplitude_to_asf(self.ampl_doublepass_default_pct / 100.)

        # defaults - sigma_x pulses
        self.freq_sigmax_ftw =  self.qubit.frequency_to_ftw(self.freq_sigmax_mhz * MHz)
        self.ampl_sigmax_asf =  self.qubit.amplitude_to_asf(self.ampl_sigmax_pct / 100.)
        self.time_sigmax_mu =   self.core.seconds_to_mu(self.time_sigmax_us * us)

        # defaults - cat
        freq_cat_center_ftw_list =  array([self.qubit.frequency_to_ftw(freq_mhz * MHz)
                                           for freq_mhz in self.freq_cat_center_mhz_list])
        freq_cat_secular_ftw_list = array([self.qubit.singlepass0.frequency_to_ftw(freq_khz * kHz)
                                           for freq_khz in self.freq_cat_secular_khz_list])
        self.ampls_cat_asf =        array([self.qubit.singlepass0.amplitude_to_asf(ampl_pct / 100.)
                                           for ampl_pct in self.ampls_cat_pct])


        '''
        CONVERT VALUES TO MACHINE UNITS - PULSES
        '''
        # cat1 values
        self.phase_cat1_antisigmax_pow =    self.qubit.turns_to_pow(self.phase_cat1_antisigmax_turns)
        self.time_cat1_bichromatic_mu = self.core.seconds_to_mu(self.time_cat1_bichromatic_us * us)
        self.phases_pulse1_cat_pow =    [self.qubit.singlepass0.turns_to_pow(phas_pow)
                                         for phas_pow in self.phases_pulse1_cat_turns]

        # inter-cat ramsey delay
        if self.enable_cat2_bichromatic:
            time_ramsey_delay_mu_list = [self.core.seconds_to_mu(time_delay_us * us)
                                         for time_delay_us in self.time_ramsey_delay_us_list]
        else:
            time_ramsey_delay_mu_list = array([0], dtype=int64)

        # cat2 values
        self.phase_cat2_sigmax_pow =        self.qubit.turns_to_pow(self.phase_cat2_sigmax_turns)
        self.phase_cat2_antisigmax_pow =    self.qubit.turns_to_pow(self.phase_cat2_antisigmax_turns)
        self.phases_cat2_cat_pow =    [self.qubit.singlepass0.turns_to_pow(phas_pow)
                                       for phas_pow in self.phases_cat2_cat_turns]

        if self.enable_cat2_bichromatic:
            time_cat2_cat_mu_list =   array([self.core.seconds_to_mu(time_us * us)
                                             for time_us in self.time_cat2_cat_us_list], dtype=int64)
            phase_cat2_cat_pow_list = array([self.qubit.singlepass0.turns_to_pow(phas_pow)
                                             for phas_pow in self.phase_cat2_cat_turns_list], dtype=int32)
        else:
            time_cat2_cat_mu_list =   array([0], dtype=int64)
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

        # readout pulse
        self.ampl_729_readout_asf =  self.qubit.amplitude_to_asf(self.ampl_729_readout_pct / 100.)

        if self.enable_729_readout:
            freq_729_readout_ftw_list =  array([self.qubit.frequency_to_ftw(freq_mhz * MHz)
                                                for freq_mhz in self.freq_729_readout_mhz_list], dtype=int32)
            time_729_readout_mu_list =   array([self.core.seconds_to_mu(time_us * us)
                                                for time_us in self.time_729_readout_us_list], dtype=int64)
        else:
            freq_729_readout_ftw_list =  array([0], dtype=int32)
            time_729_readout_mu_list =   array([0], dtype=int64)


        '''
        CREATE ATTENUATION REGISTERS
        '''
        # attenuation register - sigma_x: singlepasses set to default
        self.att_reg_sigmax = 0x00000000 | (
                (att_to_mu(self.att_sigmax_db * dB) << ((self.qubit.beam.chip_select - 4) * 8)) |
                (self.qubit.att_singlepass0_default_mu << ((self.qubit.singlepass0.chip_select - 4) * 8)) |
                (self.qubit.att_singlepass1_default_mu << ((self.qubit.singlepass1.chip_select - 4) * 8)) |
                (self.qubit.att_doublepass_inj_default_mu << ((self.qubit.doublepass_inj.chip_select - 4) * 8))
        )

        # attenuation register - bichromatic: main doublepass set to specified experiment argument value
        self.att_reg_bichromatic = 0x00000000 | (
                (att_to_mu(self.att_doublepass_default_db * dB) << ((self.qubit.beam.chip_select - 4) * 8)) |
                (att_to_mu(self.atts_cat_db[0] * dB) << ((self.qubit.singlepass0.chip_select - 4) * 8)) |
                (att_to_mu(self.atts_cat_db[1] * dB) << ((self.qubit.singlepass1.chip_select - 4) * 8)) |
                (self.qubit.att_doublepass_inj_default_mu << ((self.qubit.doublepass_inj.chip_select - 4) * 8))
        )

        # attenuation register - readout: singlepasses set to default
        self.att_reg_readout = 0x00000000 | (
                (att_to_mu(self.att_729_readout_db * dB) << ((self.qubit.beam.chip_select - 4) * 8)) |
                (self.qubit.att_singlepass0_default_mu << ((self.qubit.singlepass0.chip_select - 4) * 8)) |
                (self.qubit.att_singlepass1_default_mu << ((self.qubit.singlepass1.chip_select - 4) * 8)) |
                (self.qubit.att_doublepass_inj_default_mu << ((self.qubit.doublepass_inj.chip_select - 4) * 8))
        )


        '''
        CREATE EXPERIMENT CONFIG
        '''
        self.config_experiment_list = create_experiment_config(
            freq_cat_center_ftw_list, freq_cat_secular_ftw_list,
            time_cat2_cat_mu_list, phase_cat2_cat_pow_list,
            freq_729_readout_ftw_list, time_729_readout_mu_list,
            time_ramsey_delay_mu_list,
            config_type=int64, shuffle_config=True
        )

    def _prepare_argument_checks(self) -> TNone:
        """
        Check experiment arguments for validity.
        """
        # todo
        pass
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

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_experiment_list),
                8)


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

    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        # predeclare variables ahead of time
        time_start_mu = now_mu() & ~0x7     # store reference time for device synchronization
        ion_state = (-1, 0, int64(0))       # store ion state for adaptive readout
        herald_counter = 0                  # store herald attempts

        # MAIN LOOP
        for trial_num in range(self.repetitions):
            for config_vals in self.config_experiment_list:

                '''PREPARE & CONFIGURE'''
                # extract values from config list
                freq_cat_center_ftw =   int32(config_vals[0])
                freq_cat_secular_ftw =  int32(config_vals[1])
                time_cat2_cat_mu =      config_vals[2]
                phase_cat2_cat_pow =    int32(config_vals[3])
                freq_729_readout_ftw =  int32(config_vals[4])
                time_729_readout_mu =   config_vals[5]
                time_ramsey_delay_mu =  config_vals[6]

                # prepare phase arrays for bichromatic
                cat4_phases = [
                    self.phases_cat2_cat_pow[0] + self.phases_cat2_cat_update_dir[0] * phase_cat2_cat_pow,
                    self.phases_cat2_cat_pow[1] + self.phases_cat2_cat_update_dir[1] * phase_cat2_cat_pow,
                ]

                # clear herald counter
                herald_counter = 0

                while True:
                    # check heralding OK (otherwise execution is blocked
                    if herald_counter >= self.max_herald_attempts:
                        print("\t\tWarning: too many heralds. Moving onto next configuration.")
                        self.core.break_realtime()  # add slack
                        break

                    # add slack for execution
                    self.core.break_realtime()


                    '''INITIALIZE ION STATE'''
                    # initialize ion in S-1/2 state & SBC to ground state
                    self.initialize_subsequence.run_dma()
                    self.sidebandcool_subsequence.run_dma()

                    # set target profile to ensure we run correctly
                    self.qubit.set_profile(self.profile_729_target)
                    self.qubit.io_update()

                    # synchronize start time to coarse RTIO clock
                    time_start_mu = now_mu() & ~0x7


                    '''CAT #1'''
                    # cat1 - sigma_x (displacement vs cat)
                    if self.enable_cat1_sigmax:
                        self.pulse_sigmax(time_start_mu, 0)
                    # cat1 - bichromatic cat pulse
                    if self.enable_cat1_bichromatic:
                        self.pulse_bichromatic(time_start_mu, self.time_cat1_bichromatic_mu,
                                               self.phases_pulse1_cat_pow,
                                               freq_cat_center_ftw, freq_cat_secular_ftw)
                    # cat1 - anti-sigma_x (return to S-state)
                    if self.enable_cat1_antisigmax:
                        self.pulse_sigmax(time_start_mu, self.phase_cat1_antisigmax_pow)

                    # cat1 - force herald (to projectively disentangle spin/motion)
                    if self.enable_cat1_herald:
                        ion_state = self.readout_adaptive_subsequence.run()
                        delay_mu(self.time_adapt_read_slack_mu) # add slack following completion
                        self.pump.off()

                        # ensure dark state (flag is 0)
                        if ion_state[0] != 0:
                            herald_counter += 1 # increment herald counter to check for errors
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


                    '''RAMSEY DELAY'''
                    if self.enable_ramsey_delay:
                        delay_mu(time_ramsey_delay_mu)


                    '''CAT #2'''
                    # cat2 - sigma_x (displacement vs cat)
                    if self.enable_cat2_sigmax:
                        self.pulse_sigmax(time_start_mu, self.phase_cat2_sigmax_pow)
                    # cat2 - bichromatic cat pulse
                    if self.enable_cat2_bichromatic:
                        self.pulse_bichromatic(time_start_mu, time_cat2_cat_mu,
                                               cat4_phases,
                                               freq_cat_center_ftw, freq_cat_secular_ftw)
                    # cat2 - anti-sigma_x (return to S-state)
                    if self.enable_cat2_antisigmax:
                        self.pulse_sigmax(time_start_mu, self.phase_cat2_antisigmax_pow)

                    # cat2 - force herald (to projectively disentangle spin/motion)
                    # todo: should we be doing now_mu instead? get_rtio_counter isn't very deterministic ...
                    if self.enable_cat2_herald:
                        ion_state = self.readout_adaptive_subsequence.run()
                        delay_mu(self.time_adapt_read_slack_mu) # add slack following completion
                        self.pump.off()

                        # ensure dark state (flag is 0)
                        if ion_state[0] != 0:
                            herald_counter += 1 # increment herald counter to check for errors
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


                '''READ OUT & STORE RESULTS'''
                # only read out if no boooboo
                if herald_counter < self.max_herald_attempts:
                    # 729nm based readout (sideband ratio, rabi flopping)
                    if self.enable_729_readout:
                        self.pulse_readout(time_729_readout_mu, freq_729_readout_ftw)

                    # read out fluorescence & clean up loop
                    self.readout_subsequence.run_dma()
                    counts_res = self.readout_subsequence.fetch_count()
                else:
                    # return -1 so user knows booboo happened
                    counts_res = -1
                    self.check_termination() # check termination b/c we haven't in a while

                # store results
                self.core.break_realtime()
                self.rescue_subsequence.resuscitate()
                self.update_results(freq_cat_center_ftw,
                                    counts_res,
                                    freq_cat_secular_ftw,
                                    time_cat2_cat_mu,
                                    phase_cat2_cat_pow,
                                    freq_729_readout_ftw,
                                    time_729_readout_mu,
                                    time_ramsey_delay_mu)

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
    def pulse_readout(self, time_pulse_mu: TInt64, freq_readout_ftw: TInt32) -> TNone:
        """
        Run a readout pulse.
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
        self.qubit.cpld.set_all_att_mu(self.att_reg_readout)

        # run readout pulse
        self.qubit.singlepass0.sw.on()
        # todo: this should be off, right???
        self.qubit.singlepass1.sw.on()
        self.qubit.on()
        delay_mu(time_pulse_mu)
        self.qubit.off()
        self.qubit.singlepass1.sw.off()


    '''
    ANALYSIS
    '''
    def analyze_experiment(self):
        pass

