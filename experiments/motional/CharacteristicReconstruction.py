from artiq.experiment import *
from artiq.coredevice import ad9910
from numpy import int32, int64, array, arctan2

from LAX_exp.language import *
from LAX_exp.system.subsequences import InitializeQubit, SidebandCoolContinuousRAM, ReadoutAdaptive, RescueIon

import math
# todo: support 1D phase sweep for speed (a la SDR charread)


class CharacteristicReconstruction(LAXExperiment, Experiment):
    """
    Experiment: Characteristic Reconstruction

    Directly reconstruct the Characteristic function of a given motional state using the Fluhmann technique
    (https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.125.043602).
    """
    name = 'Characteristic Reconstruction'
    kernel_invariants = {
        # subsequences
        'initialize_subsequence', 'sidebandcool_subsequence', 'readout_subsequence', 'rescue_subsequence',

        # hardware parameters
        'ampl_doublepass_default_asf',
        'freq_sigmax_ftw', 'ampl_sigmax_asf', 'time_sigmax_mu', 'phase_antisigmax_pow',
        'time_herald_slack_mu', 'max_herald_attempts', 'time_adapt_read_slack_mu',
        'att_reg_bichromatic', 'att_reg_sigmax',

        # cat state parameters
        'ampls_cat_asf', 'time_motion_cat_mu', 'phases_motion_cat_pow', 'phase_char_axis_pow',
        'phases_char_cat_pow', 'phases_char_cat_update_dir',

        # experiment parameters
        'profile_729_SBC', 'profile_729_target', 'config_experiment_list',
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
            ram_addr_start_729=0, ram_addr_start_854=0, num_samples=200
        )
        self.initialize_subsequence =   InitializeQubit(self)
        self.readout_subsequence =      ReadoutAdaptive(self, time_bin_us=10, error_threshold=1e-2)
        self.rescue_subsequence =       RescueIon(self)

        # get relevant devices
        self.setattr_device('qubit')
        self.setattr_device('pump')
        self.setattr_device('repump_qubit')

        # set build arguments
        self._build_arguments_default()
        self._build_arguments_stateprep()
        self._build_arguments_char()

    def _build_arguments_default(self):
        """
        Build arguments for default beam parameters.
        """
        # defaults - sigma_x
        self.setattr_argument("freq_sigmax_mhz",    NumberValue(default=101.1013, precision=6, step=1, min=50., max=400., scale=1., unit="MHz"),
                              group="default.sigmax")
        self.setattr_argument("ampl_sigmax_pct",    NumberValue(default=50., precision=3, step=5, min=0.01, max=50, scale=1., unit="%"),
                              group="default.sigmax")
        self.setattr_argument("att_sigmax_db",      NumberValue(default=8., precision=1, step=0.5, min=8., max=31.5, scale=1., unit="dB"),
                              group="default.sigmax")
        self.setattr_argument("time_sigmax_us",     NumberValue(default=1.59, precision=2, step=5, min=0.1, max=10000, scale=1., unit="us"),
                              group="default.sigmax")
        self.setattr_argument("phase_antisigmax_turns",    NumberValue(default=0.25, precision=3, step=0.1, min=-1.0, max=1.0, scale=1., unit="turns"),
                              group='default.sigmax',
                              tooltip="Relative phase applied for the anti-sigma_x pulse.\n"
                                      "Note: this phase is applied via the main doublepass DDS, so values should be halved.")

        # defaults - bichromatic
        self.setattr_argument("ampl_doublepass_default_pct",    NumberValue(default=50., precision=3, step=5, min=0.01, max=50, scale=1., unit="%"),
                              group="default.bichromatic",
                              tooltip="DDS amplitude for the main doublepass during the bichromatic pulse.")
        self.setattr_argument("att_doublepass_default_db",      NumberValue(default=8., precision=1, step=0.5, min=8., max=31.5, scale=1., unit="dB"),
                              group="default.bichromatic",
                              tooltip="DDS attenuation for the main doublepass during the bichromatic pulse.")
        self.setattr_argument("ampls_cat_pct",  PYONValue([50., 50.]),
                              group='default.bichromatic',
                              tooltip="DDS amplitudes for the singlepass DDSs during the bichromatic pulses.\n"
                                      "Should be a list of [rsb_ampl_pct, bsb_ampl_pct], which are applied to [singlepass0, singlepass1].")
        self.setattr_argument("atts_cat_db",    PYONValue([13., 13.]),
                              group='default.bichromatic',
                              tooltip="DDS attenuations for the singlepass DDSs during the bichromatic pulses.\n"
                                      "Should be a list of [rsb_att_db, bsb_att_db], which are applied to [singlepass0, singlepass1].")
        self.setattr_argument("freq_cat_center_mhz_list",   Scannable(
                                                                default=[
                                                                    ExplicitScan([101.1013]),
                                                                    CenterScan(101.1013, 0.01, 0.0001, randomize=True),
                                                                    RangeScan(101.1005, 101.1021, 50, randomize=True),
                                                                ],
                                                                global_min=60., global_max=400, global_step=1,
                                                                unit="MHz", scale=1, precision=6
                                                            ),
                              group='default.bichromatic',
                              tooltip="Center frequency for the bichromatic pulses.\n"
                                      "Note: this is applied via the main doublepass DDS.\n"
                                      "The singlepass DDSs center frequencies are set as their default values from the dataset manager "
                                      "(e.g. beams.freq_mhz.freq_singlepass0_mhz).")
        self.setattr_argument("freq_cat_secular_khz_list",  Scannable(
                                                                default=[
                                                                    ExplicitScan([703.101]),
                                                                    CenterScan(703.1, 4, 0.1, randomize=True),
                                                                    RangeScan(701.0, 704.0, 50, randomize=True),
                                                                ],
                                                                global_min=0, global_max=10000, global_step=1,
                                                                unit="kHz", scale=1, precision=3
                                                            ),
                              group='default.bichromatic',
                              tooltip="Single-sided detuning frequency for the bichromatic pulses, applied via singlepass DDSs.\n"
                                      "The singlepass0 DDS is treated as the RSB, and will thus have its frequency DECREASED by this amount.\n"
                                      "Similarly, the singlepass1 DDS is treated as the BSB, and will thus have its frequency INCREASED by this amount.\n"
                                      "i.e. frequencies for [singlepass0, singlepass1] is set as [beams.freq_mhz.freq_singlepass0_mhz - freq_cat_secular_khz, "
                                      "beams.freq_mhz.freq_singlepass1_mhz + freq_cat_secular_khz].")

    def _build_arguments_stateprep(self):
        """
        Build arguments for motional state preparation.
        """
        self.setattr_argument("enable_motion_sigmax",   BooleanValue(default=False), group='motion.config',
                              tooltip="Applies a sigma_x pulse BEFORE the 1st bichromatic pulse.\n"
                                      "If sigma_x is applied (i.e. True), the bichromatic pulse creates a pure eigenstate (e.g. coherent state).\n"
                                      "If sigma_x is disabled (i.e. False), the bichromatic pulse creates a superposition state (e.g. cat state).")
        self.setattr_argument("enable_motion_cat",      BooleanValue(default=False), group='motion.config',
                              tooltip="Enables application of the 1st bichromatic pulse.\n"
                                      "Pulses are applied as [sigma_x, bichromatic, antisigma_x, herald, quench].")
        self.setattr_argument("enable_motion_antisigmax", BooleanValue(default=False), group='motion.config',
                              tooltip="Applies a sigma_x pulse AFTER the bichromatic pulse.\n"
                                      "If the sigma_x pulse is APPLIED, then this pulse disentangles spin from motion.\n"
                                      "If the sigma_x pulse is DISABLED, then this pulse simply selects whether an "
                                      "odd or even superposition is associated with the dark state.")

        self.setattr_argument("enable_motion_herald",   BooleanValue(default=True),
                              group='motion.config',
                              tooltip="Enables spin-state heralding via state-selective fluorescence. "
                                      "Heralding only progresses if the state is dark, since otherwise, the motional state is destroyed.\n"
                                      "Pulses are applied as [sigma_x, bichromatic, antisigma_x, herald, quench].\n"
                                      "Note: uses adaptive readout - ensure adaptive readout arguments are correctly set in the dataset manager.")
        self.setattr_argument("time_motion_cat_us",     NumberValue(default=100, precision=2, step=5, min=0.1, max=10000),
                              group="motion.config",
                              tooltip="Pulse time for the bichromatic pulse.")
        self.setattr_argument("phases_motion_cat_turns",    PYONValue([0., 0.]),
                              group='motion.config',
                              tooltip="Relative phases for the singlepass DDSs during the bichromatic pulse.\n"
                                      "Should be a list of [rsb_phase_turns, bsb_phase_turns].\n"
                                      "Note: these phases are applied to the singlepass DDSs, so do not need to be halved, "
                                      "unlike phases applied to the main doublepass.")

    def _build_arguments_char(self):
        """
        Build arguments for characteristic function reconstruction.
        """
        # sigma_x: select real/imag part of characteristic function
        self.setattr_argument("characteristic_axis",        EnumerationValue(['Both', 'Real', 'Imaginary'], default='Both'),
                              group='char.axis',
                              tooltip="Selects the real/imag component of the characteristic function by "
                                      "either applying a sigma_x operation (Imag), or not (Real). "
                                      "The 'Both' option enables measurement of both real and imag components within a single experiment.")
        self.setattr_argument("phase_char_axis_turns",  NumberValue(default=0.125, precision=3, step=0.1, min=-1.0, max=1.0, scale=1., unit="turns"),
                              group='char.axis',
                              tooltip="Sets the relative phase of the sigma_x operation used "
                                      "to define the real/imag axis of the characteristic function.")

        # bichromatic: characteristic readout protocol
        self.setattr_argument("phases_char_cat_turns",    PYONValue([0., 0.]),
                              group='char.read',
                              tooltip="[rsb_turns, bsb_turns]")
        self.setattr_argument("target_char_cat_phase",    EnumerationValue(['RSB', 'BSB', 'RSB-BSB', 'RSB+BSB'], default='RSB-BSB'),
                              group='char.read',
                              tooltip="todo: document")
        self.setattr_argument("time_char_cat_x_us_list",    Scannable(
                                                                default=[
                                                                    RangeScan(-50, 50, 11, randomize=True),
                                                                    ExplicitScan([100]),
                                                                ],
                                                                global_min=-100000, global_max=100000, global_step=1,
                                                                unit="us", scale=1, precision=5
                                                            ),
                              group='char.read',
                              tooltip="todo: document")
        self.setattr_argument("time_char_cat_y_us_list",    Scannable(
                                                                default=[
                                                                    RangeScan(-50, 50, 11, randomize=True),
                                                                    ExplicitScan([100]),
                                                                ],
                                                                global_min=-100000, global_max=100000, global_step=1,
                                                                unit="us", scale=1, precision=5
                                                            ),
                              group='char.read',
                              tooltip="todo: document")

    def prepare_experiment(self):
        """
        Prepare & precompute experimental values.
        """
        '''
        GENERAL SETUP/CHECK INPUT FOR ERRORS
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

        # defaults - sigma_x waveform
        self.freq_sigmax_ftw =  self.qubit.frequency_to_ftw(self.freq_sigmax_mhz * MHz)
        self.ampl_sigmax_asf =  self.qubit.amplitude_to_asf(self.ampl_sigmax_pct / 100.)
        self.time_sigmax_mu =   self.core.seconds_to_mu(self.time_sigmax_us * us)
        self.phase_antisigmax_pow = self.qubit.turns_to_pow(self.phase_antisigmax_turns)

        # defaults - cat
        freq_cat_center_ftw_list =  array([self.qubit.frequency_to_ftw(freq_mhz * MHz)
                                              for freq_mhz in self.freq_cat_center_mhz_list])
        freq_cat_secular_ftw_list = array([self.qubit.singlepass0.frequency_to_ftw(freq_khz * kHz)
                                              for freq_khz in self.freq_cat_secular_khz_list])
        self.ampls_cat_asf =    array([self.qubit.singlepass0.amplitude_to_asf(ampl_pct / 100.)
                                          for ampl_pct in self.ampls_cat_pct], dtype=int32)


        '''
        CONVERT VALUES TO MACHINE UNITS
        '''
        # motional state - cat/bichromatic
        self.time_motion_cat_mu =       self.core.seconds_to_mu(self.time_motion_cat_us * us)
        self.phases_motion_cat_pow =    [self.qubit.singlepass0.turns_to_pow(phas_pow)
                                         for phas_pow in self.phases_motion_cat_turns]

        # define characteristic axis (via sigma_x)
        self.phase_char_axis_pow =  self.qubit.turns_to_pow(self.phase_char_axis_turns)
        # configure whether real/imag/both axes of the characteristic function are to be measured
        if self.characteristic_axis == "Real":          characteristic_axis_list = [False]
        elif self.characteristic_axis == "Imaginary":   characteristic_axis_list = [True]
        elif self.characteristic_axis == "Both":        characteristic_axis_list = [True, False]

        # characteristic readout - bichromatic
        self.phases_char_cat_pow = array([self.qubit.singlepass0.turns_to_pow(phas_pow)
                                          for phas_pow in self.phases_char_cat_turns], dtype=int32)
        if self.target_char_cat_phase == 'RSB':
            self.phases_char_cat_update_dir = array([1, 0], dtype=int32)
        elif self.target_char_cat_phase == 'BSB':
            self.phases_char_cat_update_dir = array([0, 1], dtype=int32)
        elif self.target_char_cat_phase == 'RSB-BSB':
            self.phases_char_cat_update_dir = array([1, -1], dtype=int32)
        elif self.target_char_cat_phase == 'RSB+BSB':
            self.phases_char_cat_update_dir = array([1, 1], dtype=int32)


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
                (att_to_mu(self.atts_cat_db[0]) << ((self.qubit.singlepass0.chip_select - 4) * 8)) |
                (att_to_mu(self.atts_cat_db[1]) << ((self.qubit.singlepass1.chip_select - 4) * 8)) |
                (self.qubit.att_doublepass_inj_default_mu << ((self.qubit.doublepass_inj.chip_select - 4) * 8))
        )


        '''
        CREATE EXPERIMENT CONFIG
        '''
        # create sampling grid in radial coordinates
        vals_char_mu_pow_list = array([
            [
                self.core.seconds_to_mu(math.sqrt(x_us ** 2. + y_us ** 2.) * us),
                self.qubit.singlepass0.turns_to_pow(arctan2(y_us, x_us) / (2. * math.pi))
            ]
            for x_us in self.time_char_cat_x_us_list
            for y_us in self.time_char_cat_y_us_list
        ], dtype=int64)

        # create an array of values for the experiment to sweep
        self.config_experiment_list = create_experiment_config(
            freq_cat_center_ftw_list, freq_cat_secular_ftw_list,
            vals_char_mu_pow_list, characteristic_axis_list,
            config_type=int64, shuffle_config=True
        )

    def _prepare_argument_checks(self) -> TNone:
        """
        Check experiment arguments for validity.
        """
        pass
        # # ensure singlepass values are safe and valid
        # if any((ampl_pct > self.max_ampl_singlepass_pct or ampl_pct < 0.
        #         for ampl_pct in self.ampl_singlepass_default_pct_list)):
        #     raise ValueError("Singlepass amplitude outside valid range - [0., {:f}].".format(self.max_ampl_singlepass_pct))
        # if any((att_db > 31.5 or att_db < self.min_att_singlepass_db
        #         for att_db in self.att_singlepass_default_db_list)):
        #     raise ValueError("Singlepass attenuation outside valid range - [{:.1f}, 31.5].".format(self.min_att_singlepass_db))

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_experiment_list),
                6)


    '''
    MAIN SEQUENCE
    '''
    @kernel(flags={"fast-math"})
    def initialize_experiment(self) -> TNone:
        # note: no need to remove phase_autoclear from CFR1 b/c PHASE_MODE_TRACKING does it for us
        # record general subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.sidebandcool_subsequence.record_dma()
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        # instantiate relevant variables
        time_start_mu = now_mu() & ~0x7 # store reference time for device synchronization
        ion_state = (-1, 0, int64(0))   # store ion state for adaptive readout
        herald_counter = 0              # store herald attempts

        # main loop
        for trial_num in range(self.repetitions):
            # sweep exp config
            for config_vals in self.config_experiment_list:

                '''PREPARE & CONFIGURE'''
                # extract values from config list
                freq_cat_center_ftw =       int32(config_vals[0])
                freq_cat_secular_ftw =      int32(config_vals[1])
                time_char_cat_mu =          config_vals[2]
                phase_char_cat_pow =        int32(config_vals[3])
                characteristic_axis_bool =  bool(config_vals[4])

                # prepare variables for execution
                char_read_phases = [
                    self.phases_char_cat_pow[0] + self.phases_char_cat_update_dir[0] * phase_char_cat_pow,
                    self.phases_char_cat_pow[1] + self.phases_char_cat_update_dir[1] * phase_char_cat_pow,
                ]

                # clear herald counter
                herald_counter = 0

                while True:
                    # check heralding OK (otherwise execution is blocked)
                    if herald_counter >= self.max_herald_attempts:
                        print("\t\tWarning: too many heralds. Moving onto next configuration.")
                        self.core.break_realtime()  # add slack
                        break

                    # add slack for execution
                    self.core.break_realtime()


                    '''INITIALIZE ION'''
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
                    # anti-sigma_x: return to S-state
                    if self.enable_motion_antisigmax:
                        self.pulse_sigmax(time_start_mu, self.phase_antisigmax_pow, True)

                    # herald ion via state-dependent fluorescence (to projectively disentangle spin/motion)
                    if self.enable_motion_herald:
                        ion_state = self.readout_subsequence.run()
                        delay_mu(self.time_adapt_read_slack_mu)
                        self.pump.off()

                        # ensure dark state (flag is 0)
                        if ion_state[0] != 0:
                            herald_counter += 1
                            continue
                        # otherwise, add minor slack and proceed
                        at_mu(self.core.get_rtio_counter_mu() + self.time_herald_slack_mu)

                    # force break loop by default
                    break


                '''READOUT: DIRECT CHARACTERISTIC MEASUREMENT'''
                # only bother continuing with readout if no boooboo
                if herald_counter < self.max_herald_attempts:
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

                    # read out fluorescence & clean up loop
                    ion_state = self.readout_subsequence.run()
                else:
                    # return -1 so user knows booboo happened
                    ion_state = (-1, -1, int64(0))
                    self.check_termination() # check termination b/c we haven't in a while
                    self.core.break_realtime()

                # save results
                self.rescue_subsequence.resuscitate()
                self.update_results(freq_cat_center_ftw,
                                    ion_state[0],
                                    freq_cat_secular_ftw,
                                    time_char_cat_mu,
                                    phase_char_cat_pow,
                                    characteristic_axis_bool)

            # rescue ion as needed & support graceful termination
            self.core.break_realtime()
            self.rescue_subsequence.run(trial_num)
            self.check_termination()


    '''
    HELPER FUNCTIONS
    '''
    @kernel(flags={"fast-math"})
    def pulse_sigmax(self, time_start_mu: TInt64 = -1, phas_pow: TInt32 = 0x0, is_real: TBool = False) -> TNone:
        """
        Run a phase-coherent sigma_x pulse on the qubit.
        :param time_start_mu: fiducial timestamp for initial start reference (in machine units).
        :param phas_pow: relative phase offset for the beam.
        :param is_real: whether to actually run the pulse (True) or a dummy pulse (False).
        """
        # set up relevant beam waveforms
        self.qubit.set_mu(
            self.freq_sigmax_ftw, asf=self.ampl_sigmax_asf, pow_=phas_pow,
            profile=self.profile_729_target, phase_mode=ad9910.PHASE_MODE_TRACKING, ref_time_mu=time_start_mu
        )
        self.qubit.singlepass0.set_mu(
            self.qubit.freq_singlepass0_default_ftw, asf=self.qubit.ampl_singlepass0_default_asf, pow_=0,
            profile=self.profile_729_target, phase_mode=ad9910.PHASE_MODE_TRACKING, ref_time_mu=time_start_mu
        )
        self.qubit.singlepass1.set_mu(
            self.qubit.freq_singlepass1_default_ftw, asf=self.qubit.ampl_singlepass1_default_asf, pow_=0,
            profile=self.profile_729_target, phase_mode=ad9910.PHASE_MODE_TRACKING, ref_time_mu=time_start_mu
        )
        self.qubit.cpld.set_all_att_mu(self.att_reg_sigmax)

        # run pulse
        self.qubit.singlepass0.sw.on()
        self.qubit.singlepass1.sw.on()

        # only fire sigma_x pulse if measuring Re[\chi(\beta)]
        if is_real: self.qubit.on()
        # otherwise, run dummy pulse to keep sequence as similar as possible
        else:       self.qubit.off()
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
        self.qubit.set_mu(
            freq_carrier_ftw, asf=self.ampl_doublepass_default_asf, pow_=0,
            profile=self.profile_729_target, phase_mode=ad9910.PHASE_MODE_TRACKING, ref_time_mu=time_start_mu
        )
        self.qubit.singlepass0.set_mu(
            self.qubit.freq_singlepass0_default_ftw - freq_secular_ftw, asf=self.ampls_cat_asf[0],
            pow_=phas_pow_list[0], profile=self.profile_729_target,
            phase_mode=ad9910.PHASE_MODE_TRACKING, ref_time_mu=time_start_mu
        )
        self.qubit.singlepass1.set_mu(
            self.qubit.freq_singlepass0_default_ftw + freq_secular_ftw, asf=self.ampls_cat_asf[1],
            pow_=phas_pow_list[1], profile=self.profile_729_target,
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

