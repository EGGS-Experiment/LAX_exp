from artiq.experiment import *
from artiq.coredevice import ad9910

from numpy import array, int32, int64

from LAX_exp.language import *
from LAX_exp.system.subsequences import (
    InitializeQubit, SidebandCoolContinuousRAM, Readout,
    RescueIon,
)
from LAX_exp.extensions import *
from LAX_exp.system.objects.dds_pulse_shaper import DDSPulseShaper
from LAX_exp.system.objects.PulseShaper import available_pulse_shapes


class MolmerSorensen(LAXExperiment, Experiment):
    """
    Experiment: Molmer Sorensen Gate


    """
    name = 'Molmer Sorensen Gate'

    kernel_invariants = {
        # subsequences & objects
        'initialize_subsequence', 'sidebandcool_subsequence', 'readout_subsequence',

        # configs
        'profile_729_SBC', 'profile_729_readout',
        'profile_729_ms', 'profile_729_parity',
        'config_experiment_list',
    }

    def build_experiment(self):

        # allocate relevant beam profiles
        self.profile_729_SBC = 1
        self.profile_729_readout = 2
        self.profile_729_ms = 3
        self.profile_729_parity = 4

        # devices
        self.setattr_device('qubit')
        self.setattr_device('pump')
        self.setattr_device('repump_qubit')
        self.setattr_device('ttl8')

        self.setattr_argument("repetitions", NumberValue(default=50, precision=0, step=1, min=1, max=100000))

        self._build_arguments_ms()
        self._build_arguments_parity()
        self._build_arguments_dynamical_decoupling()
        self._build_arguments_intensity_servo()

        '''
        Build subsequences
        '''
        self.sidebandcool_subsequence = SidebandCoolContinuousRAM(
            self, profile_729=self.profile_729_SBC, profile_854=3,
            ram_addr_start_729=0, ram_addr_start_854=0, num_samples=200
        )
        self.initialize_subsequence = InitializeQubit(self)
        self.readout_subsequence = Readout(self)
        self.rescue_subsequence = RescueIon(self)

    def _build_arguments_ms(self):

        _argstr = 'ms'

        self.setattr_argument('time_gate_us_list', Scannable(default=[ExplicitScan([50]),
                                                                      RangeScan(0.01, 100, 100)],
                                                             global_min=0.01, global_max=1e6,
                                                             global_step=0.01, unit="us",
                                                             precision=4, scale=1),
                              group=_argstr, tooltip='time to apply birchromatic interaction for'
                              )

        self.setattr_argument('freq_carrier_ms_mhz', NumberValue(default=101.3,
                                                                 min=2., max=110,
                                                                 step=0.01, unit="MHz",
                                                                 precision=4, scale=1),
                              group=_argstr, tooltip='carrier frequency in MHz'
                              )

        self.setattr_argument('freq_detuning_ms_khz_list', Scannable(default=[ExplicitScan([40]),
                                                                              RangeScan(0.01, 100, 100)],
                                                                     global_min=0.01, global_max=1e6,
                                                                     global_step=0.01, unit="kHz",
                                                                     precision=4, scale=1),
                              group=_argstr, tooltip='detuning from carrier'
                              )

        self.setattr_argument("ampls_ms_pct", PYONValue([50., 50.]), group=_argstr,
                              tooltip="DDS amplitudes for the singlepass DDSs during the bichromatic pulses.\n"
                                      "Should be a list of [rsb_ampl_pct, bsb_ampl_pct], which are applied to [singlepass1, singlepass2].")
        self.setattr_argument("atts_ms_db", PYONValue([11., 11.]), group=_argstr,
                              tooltip="DDS attenuations for the singlepass DDSs during the bichromatic pulses.\n"
                                      "Should be a list of [rsb_att_db, bsb_att_db], which are applied to [singlepass1, singlepass2].")

        # waveform - pulse shaping
        self.setattr_argument("enable_pulse_shaping", BooleanValue(default=False),
                              group=_argstr,
                              tooltip="Applies pulse shaping to the edges of the ms pulse.")
        self.setattr_argument("type_pulse_shape",
                              EnumerationValue(list(available_pulse_shapes.keys()), default='sine_squared'),
                              group=_argstr,
                              tooltip="Pulse shape type to be used.")

    def _build_arguments_dynamical_decoupling(self):
        _argstr = 'dynamical_decoupling'
        self.setattr_argument("enable_dyanmical_decoupling", BooleanValue(default=False),
                              group=_argstr,
                              tooltip="Enable continuous dynamical decoupling by turning on a carrier tone "
                                      "during the ms gate.")

        self.setattr_argument('dynamical_decoupling_turns_list',
                              Scannable(default=[
                                  ExplicitScan([0.]),
                                  RangeScan(0, 1, 11),
                                  CenterScan(0.5, 1, 0.1)], unit='turns',
                                  global_min=0., global_max=2., global_step=0.1, precision=3, scale=1.0),
                              tooltip='Phase to apply to the dynamical decoupling tone',
                              group=_argstr)

        self.setattr_argument('att_dynamical_decoupling_db', NumberValue(default=5., min=12., max=31.5,
                                                                         step=0.5, precision=1, scale=1, unit='dB'),
                              tooltip='attenuation to apply to the singlepass for ', group=_argstr)

        self.setattr_argument('ampl_dynamical_decoupling_pct', NumberValue(default=50., min=0.01, max=50.,
                                                                           step=0.01, precision=3, scale=1, unit='%'),
                              tooltip='amplitude of rf tone used for continuous dynamical decoupling', group=_argstr)

    def _build_arguments_parity(self):

        _argstr = 'parity_pulse'

        self.setattr_argument('enable_parity_pulse', BooleanValue(False),
                              tooltip='Enable parity analysis pulse after applying the MS gate',
                              group=_argstr)

        self.setattr_argument('time_parity_pulse_us',
                              NumberValue(default=2, min=0.01, max=1e6,
                                          step=0.01, scale=1, precision=4, unit='us'),
                              tooltip='', group=_argstr)

        self.setattr_argument('phase_parity_pulse_turns_list',
                              Scannable(default=[
                                  ExplicitScan([0.]),
                                  RangeScan(0, 1, 11),
                                  CenterScan(0.5, 1, 0.1)], unit='turns',
                                  global_min=0., global_max=2., global_step=0.1, precision=3, scale=1.0),
                              tooltip='',
                              group=_argstr)

        self.setattr_argument('att_parity_pulse_db', NumberValue(default=5., min=5., max=31.5,
                                                                 step=0.5, precision=1, scale=1, unit='dB'),
                              tooltip='', group=_argstr)

        self.setattr_argument('ampl_parity_pulse_pct', NumberValue(default=50., min=0.01, max=50.,
                                                                   step=0.01, precision=3, scale=1, unit='%'),
                              tooltip='', group=_argstr)

    def _build_arguments_intensity_servo(self):
        """
        Build arguments for intensity servo hold between shots
        """
        _argstr = "servo_relock"
        self.setattr_argument('enable_servo_relock', BooleanValue(default=False), group=_argstr)
        self.setattr_argument('time_servo_relock_us', NumberValue(default=2000, precision=3, step=1, min=1,
                                                                  max=10000, scale=1., unit='us'), group=_argstr)

    def prepare_experiment(self):
        '''
        Prepare all necessary functions
        '''

        '''
        Convert MS arguments
        '''
        time_gate_mu_list = [self.core.seconds_to_mu(time_gate_us * us) for time_gate_us in
                             self.time_gate_us_list]
        self.freq_carrier_ms_mu = self.qubit.frequency_to_ftw(self.freq_carrier_ms_mhz * MHz)
        freq_detuning_ms_mu_list = [self.qubit.frequency_to_ftw(detuning_khz * kHz) for detuning_khz in
                                    self.freq_detuning_ms_khz_list]
        self.ampls_ms_asf = array([self.qubit.amplitude_to_asf(ampl_pct / 100.)
                                   for ampl_pct in self.ampls_ms_pct])
        self.time_servo_relock_mu = self.core.seconds_to_mu(self.time_servo_relock_us * us)

        ### Build Pulse Shaper
        if self.enable_pulse_shaping:
            pulse_shape = self.type_pulse_shape
        else:
            pulse_shape = 'square'
        self.dds_pulse_shaper_ms = DDSPulseShaper(self, dds_targets=self.qubit.beam,
                                                  ram_profile=self.profile_729_ms,
                                                  ram_addr_start=502, num_samples=200,
                                                  ampl_max_pcts=self.qubit.asf_to_amplitude(
                                                      self.qubit.ampl_qubit_asf) * 100,
                                                  pulse_shapes='square',
                                                  phase_autoclear=1)
        self.dds_pulse_shaper_ms.add_dds_target(self.qubit.singlepass1,
                                                ampl_max_pct=self.ampls_ms_pct[0],
                                                pulse_shape=pulse_shape,
                                                phase_autoclear=1,
                                                )
        self.dds_pulse_shaper_ms.add_dds_target(self.qubit.singlepass2,
                                                ampl_max_pct=self.ampls_ms_pct[1],
                                                pulse_shape=pulse_shape,
                                                phase_autoclear=1,
                                                )

        '''
        Convert Parity Arguments
        '''

        self.time_parity_pulse_mu = self.core.seconds_to_mu(self.time_parity_pulse_us * us)
        if self.enable_parity_pulse:
            phase_parity_pulse_pow_list = [self.qubit.turns_to_pow(phase_parity_pulse)
                                           for phase_parity_pulse in self.phase_parity_pulse_turns_list]
        else:
            phase_parity_pulse_pow_list = [0]
        self.ampl_parity_pulse_asf = self.qubit.amplitude_to_asf(self.ampl_parity_pulse_pct / 100.)

        '''
        Convert DD Arguments
        '''
        if self.enable_dynamical_decoupling:
            dynamical_decoupling_pow_list = [self.qubit.turns_to_pow(phase_dd)
                                             for phase_dd in self.phase_dynamical_decoupling_turns_list]
        else:
            dynamical_decoupling_pow_list = [0]
        self.ampl_dynamical_decoupling_asf = self.qubit.amplitude_to_asf(self.ampl_dynamical_decoupling_asf / 100.)
        self.phase_dd_phase_shift_pow = self.qubit.turns_to_pow(0.5)
        self.urukul_dd_reset_time = int64(650)  # half the time for urukul to implement set_mu for DD phase shifting

        '''
        Prepare Attenuation Register
        '''

        # attenuation register - bichromatic
        self.att_reg_ms = 0x00000000 | (
                (self.qubit.att_qubit_mu << ((self.qubit.beam.chip_select - 4) * 8)) |
                (att_to_mu(self.att_dynamical_decoupling_db) << ((self.qubit.singlepass0.chip_select - 4) * 8)) |
                (att_to_mu(self.atts_ms_db[0] * dB) << ((self.qubit.singlepass1.chip_select - 4) * 8)) |
                (att_to_mu(self.atts_ms_db[1] * dB) << ((self.qubit.singlepass2.chip_select - 4) * 8))
        )

        # attenuation register - parity
        self.att_reg_parity = 0x00000000 | (
                (self.qubit.att_qubit_mu << ((self.qubit.beam.chip_select - 4) * 8)) |
                (att_to_mu(self.att_parity_pulse_db) << ((self.qubit.singlepass0.chip_select - 4) * 8)) |
                (att_to_mu(31.5 * dB) << ((self.qubit.singlepass1.chip_select - 4) * 8)) |
                (att_to_mu(31.5 * dB) << ((self.qubit.singlepass2.chip_select - 4) * 8))
        )

        self.config_experiment_list = create_experiment_config(
            time_gate_mu_list,
            freq_detuning_ms_mu_list,
            phase_parity_pulse_pow_list,
            dynamical_decoupling_pow_list,

            config_type=float, shuffle_config=True
        )

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_experiment_list),
                5)

    @kernel(flags={'fast_math'})
    def initialize_experiment(self) -> TNone:
        self.initialize_subsequence.record_dma()
        self.sidebandcool_subsequence.record_dma()
        self.readout_subsequence.record_dma()
        self.core.break_realtime()

        # initialize pulse shaper
        self.dds_pulse_shaper_ms.sequence_initialize()
        delay_mu(50000)

    @kernel(flags={'fast_math'})
    def run_main(self) -> TNone:

        for trial_num in range(self.repetitions):
            for config_vals in self.config_experiment_list:
                time_gate_ms_mu = int64(config_vals[0])
                freq_detuning_ms_mu = int32(config_vals[1])
                phase_parity_pulse_pow = int32(config_vals[2])
                dynamical_decoupling_pulse_pow = int32(config_vals[3])

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
                time_acutal_mu = [0]

                # # set cfr1 so we clear phases of all urukul0 channels on next io_update
                self.qubit.off()
                self.qubit.singlepass0_off()
                self.qubit.singlepass1_off()
                self.qubit.singlepass2_off()
                self.qubit.cpld.set_all_att_mu(self.att_reg_ms)

                # setup ms tones
                self.dds_pulse_shaper_ms.dds_targets[0].set_ftw(self.freq_carrier_ms_mu)
                self.dds_pulse_shaper_ms.dds_targets[0].set_pow(0)

                self.dds_pulse_shaper_ms.dds_targets[1].set_ftw(
                    self.qubit.freq_singlepass1_default_ftw - freq_detuning_ms_mu)
                self.dds_pulse_shaper_ms.dds_targets[1].set_pow(0)

                self.dds_pulse_shaper_ms.dds_targets[2].set_ftw(
                    self.qubit.freq_singlepass2_default_ftw + freq_detuning_ms_mu)
                self.dds_pulse_shaper_ms.dds_targets[2].set_pow(0)

                # configure singlepass0 for continuous DD
                self.qubit.singlepass0.set_mu(
                    self.qubit.freq_singlepass0_default_ftw,
                    asf=self.ampl_dynamical_decoupling_asf,
                    pow_=dynamical_decoupling_pulse_pow,
                    profile=self.profile_729_ms,
                    phase_mode=ad9910.PHASE_MODE_CONTINUOUS,
                )

                if self.enable_parity_pulse:
                    self.setup_parity_tones(phase_parity_pulse_pow)

                time_actual_mu = self.dds_pulse_shaper_ms.configure_train_all_dds([time_gate_ms_mu])
                self.qubit.singlepass0.set_cfr1(phase_autoclear=1)

                # # # pulse to clear phase accumulator
                self.qubit.io_update()
                # reset cfr1 so we no longer clear phases on io_update
                self.dds_pulse_shaper_ms.dds_targets[0].set_cfr1(ram_enable=1,
                                                                 phase_autoclear=0,
                                                                 ram_destination=ad9910.RAM_DEST_ASF)
                self.dds_pulse_shaper_ms.dds_targets[1].set_cfr1(ram_enable=1,
                                                                 phase_autoclear=0,
                                                                 ram_destination=ad9910.RAM_DEST_ASF)
                self.dds_pulse_shaper_ms.dds_targets[2].set_cfr1(ram_enable=1,
                                                                 phase_autoclear=0,
                                                                 ram_destination=ad9910.RAM_DEST_ASF)
                self.qubit.singlepass0.set_cfr1()
                self.qubit.io_update()

                # ensure we are on the right profile
                self.qubit.set_profile(self.profile_729_ms)
                delay_mu(4000)

                self.ttl8.on()
                """
                RUN MS GATE
                """
                # make sure external switch is open
                self.qubit.rf_switch.off()
                delay_mu(TIME_ZASWA2_SWITCH_DELAY_MU)

                if time_acutal_mu[0] - self.urukul_dd_reset_time < 8:
                    time_gate_dd_mu = 8
                else:
                    time_gate_dd_mu = time_actual_mu[0] - self.urukul_dd_reset_time

                if self.enable_dynamical_decoupling:
                    self.qubit.singlepass0_on()
                with parallel:
                    self.dds_pulse_shaper_ms.run_train_all_dds()
                    if self.enable_dynamical_decoupling:
                        # see time to set profile
                        delay_mu(time_gate_dd_mu)
                        self.qubit.singlepass0_off()
                        self.qubit.singlepass0.set_mu(
                            self.qubit.freq_singlepass0_default_ftw,
                            asf=self.ampl_dynamical_decoupling_asf,
                            pow_=dynamical_decoupling_pulse_pow + self.phase_dd_phase_shift_pow,
                            profile=self.profile_729_ms,
                            phase_mode=ad9910.PHASE_MODE_CONTINUOUS,
                        )
                        self.qubit.singlepass0_on()
                # make sure external switch is closed
                self.qubit.rf_switch.on()
                delay_mu(TIME_ZASWA2_SWITCH_DELAY_MU)

                # ensure ram enabled
                for dds_target in self.dds_pulse_shaper_ms.dds_targets:
                    dds_target.set_cfr1(ram_enable=0)

                """
                RUN PARITY PULSE
                """
                if self.enable_parity_pulse:
                    self.qubit.cpld.set_all_att_mu(self.att_reg_parity)
                    self.qubit.set_profile(self.profile_729_parity)
                    at_mu((now_mu() + 8) & ~7)
                    self.qubit.io_update()
                    delay_mu(2000)

                    self.qubit.on()
                    self.qubit.singlepass0_on()

                    delay_mu(self.time_parity_pulse_mu)

                    self.qubit.off()
                    self.qubit.singlepass0_off()

                self.ttl8.off()

                # # read out fluorescence & clean up loop
                self.readout_subsequence.run_dma()
                counts_res = self.readout_subsequence.fetch_count()

                self.update_results(
                    time_actual_mu[0],
                    counts_res,
                    phase_parity_pulse_pow,
                    freq_detuning_ms_mu
                )

                self.dds_pulse_shaper_ms.sequence_cleanup()
                self.check_termination()

    @kernel(flags={'fast_math'})
    def setup_ms_tones(self, freq_detuning_mu) -> TNone:
        # ensure all beams are off
        self.qubit.off()
        self.qubit.singlepass0_off()
        self.qubit.singlepass1_off()
        self.qubit.singlepass2_off()

        self.qubit.set_mu(
            self.freq_carrier_ms_mu, asf=self.qubit.ampl_qubit_asf,
            pow_=0, profile=self.profile_729_ms,
            phase_mode=ad9910.PHASE_MODE_CONTINUOUS
        )
        self.qubit.singlepass1.set_mu(
            self.qubit.freq_singlepass1_default_ftw - freq_detuning_mu,
            asf=self.ampls_ms_asf[0],
            pow_=0,
            profile=self.profile_729_ms,
            phase_mode=ad9910.PHASE_MODE_CONTINUOUS,
        )
        self.qubit.singlepass2.set_mu(
            self.qubit.freq_singlepass2_default_ftw + freq_detuning_mu,
            asf=self.ampls_ms_asf[1],
            pow_=0,
            profile=self.profile_729_ms,
            phase_mode=ad9910.PHASE_MODE_CONTINUOUS,
        )

    @kernel(flags={'fast_math'})
    def setup_parity_tones(self, phase_parity_pulse_pow) -> TNone:
        # ensure all beams are off
        self.qubit.off()
        self.qubit.singlepass0_off()
        self.qubit.singlepass1_off()
        self.qubit.singlepass2_off()

        self.qubit.set_mu(
            self.freq_carrier_ms_mu, asf=self.qubit.ampl_qubit_asf,
            pow_=0, profile=self.profile_729_parity,
            phase_mode=ad9910.PHASE_MODE_CONTINUOUS
        )
        self.qubit.singlepass0.set_mu(
            self.qubit.freq_singlepass0_default_ftw,
            asf=self.ampl_parity_pulse_asf,
            pow_=phase_parity_pulse_pow,
            profile=self.profile_729_parity,
            phase_mode=ad9910.PHASE_MODE_CONTINUOUS,
        )
