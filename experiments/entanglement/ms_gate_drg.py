from artiq.experiment import *
from artiq.coredevice import ad9910
from artiq.coredevice import urukul
from common.lib.servers.Pulser.pulse_sequences.pulse_sequences_config import dds729DP

from numpy import array, int32, int64

from LAX_exp.language import *
from LAX_exp.system.subsequences import (
    InitializeQubit, SidebandCoolContinuousRAM, Readout,
    RescueIon,
)
from LAX_exp.extensions import *
from LAX_exp.system.objects.dds_ramper import DDSRamper


class MolmerSorensenDRG(LAXExperiment, Experiment):
    """
    Experiment: Molmer Sorensen Gate DRG


    """
    name = 'Molmer Sorensen Gate Drg'

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

        self.setattr_device('trigger_line')

        self.setattr_argument("repetitions", NumberValue(default=50, precision=0, step=1, min=1, max=100000))

        self.setattr_argument("enable_linetrigger", BooleanValue(default=False),
                              tooltip="Trigger the beginning of each shot from the AC line.")

        self._build_arguemnts_ion_parameters()
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

    def _build_arguemnts_ion_parameters(self):
        _argstr = 'ion_parameters'


        self.setattr_argument('freq_carrier_mhz', NumberValue(default=100.3,
                                                            min=90., max=110,
                                                            step=0.01, unit="MHz",
                                                            precision=4, scale=1),
                      group=_argstr, tooltip='carrier frequency in MHz'
                      )


        self.setattr_argument('freq_secular_khz', NumberValue(default=710,
                                                           min=500., max=2000,
                                                           step=0.01, unit="kHz",
                                                           precision=4, scale=1),
                      group=_argstr, tooltip='secular frequency in kHz'
                      )


    def _build_arguments_ms(self):

        _argstr = 'ms'

        self.setattr_argument('time_gate_us_list', Scannable(default=[ExplicitScan([50]),
                                                                      RangeScan(1, 100, 100)],
                                                             global_min=6.5, global_max=1e6,
                                                             global_step=0.01, unit="us",
                                                             precision=4, scale=1),
                              group=_argstr, tooltip='flat time to apply birchromatic interaction for'
                              )

        self.setattr_argument('freq_carrier_detuning_ms_khz_list', Scannable(default=[ExplicitScan([0.0]),
                                                               RangeScan(-1, 1, 100)],
                                                                 global_min=-200., global_max=200,
                                                                 global_step=0.01, unit="kHz",
                                                                 precision=4, scale=1),
                              group=_argstr, tooltip='detuning from carrier frequency in kHz'
                              )

        self.setattr_argument('freq_secular_detuning_ms_khz_list', Scannable(default=[ExplicitScan([40]),
                                                                              RangeScan(0.01, 100, 100)],
                                                                     global_min=0.01, global_max=1e6,
                                                                     global_step=0.01, unit="kHz",
                                                                     precision=4, scale=1),
                              group=_argstr, tooltip='detuning from secular in kHz'
                              )

        self.setattr_argument("ampls_ms_pct", PYONValue([50., 50.]), group=_argstr,
                              tooltip="DDS amplitudes for the singlepass DDSs during the bichromatic pulses.\n"
                                      "Should be a list of [rsb_ampl_pct, bsb_ampl_pct], which are applied to [singlepass1, singlepass2].")
        self.setattr_argument("atts_ms_db", PYONValue([11., 11.]), group=_argstr,
                              tooltip="DDS attenuations for the singlepass DDSs during the bichromatic pulses.\n"
                                      "Should be a list of [rsb_att_db, bsb_att_db], which are applied to [singlepass1, singlepass2].")

    def _build_arguments_dynamical_decoupling(self):
        _argstr = 'dynamical_decoupling'
        self.setattr_argument("enable_dynamical_decoupling", BooleanValue(default=False),
                              group=_argstr,
                              tooltip="Enable continuous dynamical decoupling by turning on a carrier tone "
                                      "during the ms gate.")

        self.setattr_argument('phase_dynamical_decoupling_turns_list',
                              Scannable(default=[
                                  ExplicitScan([0.]),
                                  RangeScan(0, 1, 11),
                                  CenterScan(0.5, 1, 0.1)], unit='turns',
                                  global_min=0., global_max=2., global_step=0.1, precision=3, scale=1.0),
                              tooltip='Phase to apply to the dynamical decoupling tone',
                              group=_argstr)

        self.setattr_argument('att_dynamical_decoupling_db', NumberValue(default=25., min=12., max=31.5,
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

        self.setattr_argument('time_parity_pulse_us_list', Scannable(default=[ExplicitScan([1.4]),
                                                                      RangeScan(1, 100, 100)],
                                                             global_min=0.01, global_max=1e6,
                                                             global_step=0.01, unit="us",
                                                             precision=4, scale=1),
                              group=_argstr, tooltip='parity pulse time'
                              )

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
        Convert ion parameters
        '''
        self.freq_carrier_ftw = self.qubit.frequency_to_ftw(self.freq_carrier_mhz * MHz)
        self.freq_secular_ftw = self.qubit.frequency_to_ftw(self.freq_secular_khz * kHz)

        '''
        Convert MS arguments
        '''
        time_gate_mu_list = [self.core.seconds_to_mu(time_gate_us * us) for time_gate_us in
                             self.time_gate_us_list]

        # divide by two as it goes to the doublepass
        freq_carrier_detuning_ms_ftw_list = [self.qubit.frequency_to_ftw(freq_carrier_detuning_ms_khz/2 * kHz) for freq_carrier_detuning_ms_khz in
                                   self.freq_carrier_detuning_ms_khz_list]
        freq_secular_detuning_ms_ftw_list = [self.qubit.frequency_to_ftw(freq_secular_detuning_ms_khz * kHz) for freq_secular_detuning_ms_khz in
                                    self.freq_secular_detuning_ms_khz_list]
        self.ampls_ms_asf = array([self.qubit.amplitude_to_asf(ampl_pct / 100.)
                                   for ampl_pct in self.ampls_ms_pct])
        self.time_servo_relock_mu = self.core.seconds_to_mu(self.time_servo_relock_us * us)

        ### Build Pulse Shaper
        self.dds_ramper_ms = DDSRamper(self, self.qubit.singlepass1,
                                        num_samples = 100,
                                       ramp_dest=2,
                                       data_high=self.ampls_ms_pct[0] / 100.,
                                       data_low=0)

        self.dds_ramper_ms.add_dds_target(self.qubit.singlepass2,
                                          ramp_dest=2,
                                          data_high=self.ampls_ms_pct[1] / 100.,
                                          data_low=0)

        '''
        Convert Parity Arguments
        '''
        if self.enable_parity_pulse:
            phase_parity_pulse_pow_list = [self.qubit.turns_to_pow(phase_parity_pulse)
                                           for phase_parity_pulse in self.phase_parity_pulse_turns_list]

            time_parity_pulse_mu_list = [self.core.seconds_to_mu(time_parity_pulse_us * us) for time_parity_pulse_us
                                         in self.time_parity_pulse_us_list]
        else:
            phase_parity_pulse_pow_list = [0]
            time_parity_pulse_mu_list = [0]
        self.ampl_parity_pulse_asf = self.qubit.amplitude_to_asf(self.ampl_parity_pulse_pct / 100.)

        '''
        Convert DD Arguments
        '''
        if self.enable_dynamical_decoupling:
            phase_dynamical_decoupling_pow_list = [self.qubit.turns_to_pow(phase_dd)
                                             for phase_dd in self.phase_dynamical_decoupling_turns_list]
        else:
            phase_dynamical_decoupling_pow_list = [0]
        self.ampl_dynamical_decoupling_asf = self.qubit.amplitude_to_asf(self.ampl_dynamical_decoupling_pct / 100.)
        self.phase_dd_phase_shift_pow = self.qubit.turns_to_pow(0.5)
        self.time_urukul_reset_mu = int64(1250)  # half the time for urukul to implement set_mu for DD phase shifting

        '''
        Prepare Attenuation Register
        '''

        # attenuation register - bichromatic
        if self.enable_dynamical_decoupling:
            self.att_reg_ms = 0x00000000 | (
                    (self.qubit.att_qubit_mu << ((self.qubit.beam.chip_select - 4) * 8)) |
                    (att_to_mu(self.att_dynamical_decoupling_db) << ((self.qubit.singlepass0.chip_select - 4) * 8)) |
                    (att_to_mu(self.atts_ms_db[0] * dB) << ((self.qubit.singlepass1.chip_select - 4) * 8)) |
                    (att_to_mu(self.atts_ms_db[1] * dB) << ((self.qubit.singlepass2.chip_select - 4) * 8))
            )

        else:
            self.att_reg_ms = 0x00000000 | (
                    (self.qubit.att_qubit_mu << ((self.qubit.beam.chip_select - 4) * 8)) |
                    (att_to_mu(31.5*dB) << ((self.qubit.singlepass0.chip_select - 4) * 8)) |
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
            freq_secular_detuning_ms_ftw_list,
            phase_parity_pulse_pow_list,
            phase_dynamical_decoupling_pow_list,
            freq_carrier_detuning_ms_ftw_list,
            time_parity_pulse_mu_list,

            config_type=float, shuffle_config=True
        )

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_experiment_list),
                7)

    @kernel(flags={'fast_math'})
    def initialize_experiment(self) -> TNone:
        self.initialize_subsequence.record_dma()
        self.sidebandcool_subsequence.record_dma()
        self.readout_subsequence.record_dma()
        self.core.break_realtime()

        # initialize pulse shaper
        self.dds_ramper_ms.sequence_initialize()
        delay_mu(50000)


    @kernel(flags={'fast_math'})
    def run_main(self) -> TNone:

        for trial_num in range(self.repetitions):
            for config_vals in self.config_experiment_list:
                time_gate_ms_mu = int64(config_vals[0])
                freq_secular_detuning_ms_ftw = int32(config_vals[1])
                phase_parity_pulse_pow = int32(config_vals[2])
                phase_dynamical_decoupling_pow = int32(config_vals[3])
                freq_carrier_detuning_ms_ftw = int32(config_vals[4])
                time_parity_pulse_mu = int64(config_vals[5])

                self.core.break_realtime()  # add slack for execution
                delay_mu(125000)  # add even more slack lol

                # wait for linetrigger
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

                # # set cfr1 so we clear phases of all urukul0 channels on next io_update
                self.qubit.off()
                self.qubit.singlepass0_off()
                self.qubit.singlepass1_off()
                self.qubit.singlepass2_off()
                self.qubit.cpld.set_all_att_mu(self.att_reg_ms)

                self.setup_ms_tones(freq_carrier_detuning_ms_ftw, freq_secular_detuning_ms_ftw,
                                    phase_dynamical_decoupling_pow)

                # setup ms tones
                if self.enable_parity_pulse:
                    self.setup_parity_tones(phase_parity_pulse_pow)


                # # # pulse to clear phase accumulator
                self.dds_ramper_ms.clear_phase_accumulator_all_dds()
                self.qubit.set_cfr1(phase_autoclear=1)
                self.qubit.singlepass0.set_cfr1(phase_autoclear=1)
                at_mu((now_mu() + 8) & ~7) # synchronize to clock
                self.qubit.io_update()
                # reset cfr1 so we no longer clear phases on io_update
                self.dds_ramper_ms.reset_cfr1_all_dds()
                self.qubit.set_cfr1(phase_autoclear=0)
                self.qubit.singlepass0.set_cfr1(phase_autoclear=0)
                self.qubit.io_update()
                delay_mu(1000)

                # ensure we are on the right profile
                self.qubit.set_profile(self.profile_729_ms)
                self.dds_ramper_ms.configure_ramp_all_dds(1)

                if self.enable_dynamical_decoupling:
                    # scope says this should be 1328 ns - investigate
                    time_gate_dd_mu = (time_gate_ms_mu - self.time_urukul_reset_mu - self.dds_ramper_ms.ramp_firing_delay) >> 1
                    if time_gate_dd_mu < 0:
                        time_gate_dd_mu = 8
                else:
                    time_gate_dd_mu = (time_gate_ms_mu - self.dds_ramper_ms.ramp_firing_delay) >> 1

                """
                RUN MS GATE
                """
                self.ttl8.on()

                '''RISING PORITION OF PULSE'''
                self.qubit.on()
                # start ramp-up when coarse aligned to SYNC_CLK for determinacy
                self.dds_ramper_ms.run_ramp_all_dds()
                delay_mu(self.dds_ramper_ms.drg_time_ramp_mu[0])

                '''FLAT TOP PORITION OF PULSE'''
                if self.enable_dynamical_decoupling:
                    delay_mu(self.dds_ramper_ms.ramp_firing_delay)
                    self.qubit.singlepass0_on()
                delay_mu(time_gate_dd_mu)
                if self.enable_dynamical_decoupling:
                    self.qubit.singlepass0_off()
                    self.qubit.singlepass0.set_mu(
                            self.qubit.freq_singlepass0_default_ftw,
                            asf=self.ampl_parity_pulse_asf,
                            pow_=phase_dynamical_decoupling_pow + self.phase_dd_phase_shift_pow,
                            profile=self.profile_729_ms,
                            phase_mode=ad9910.PHASE_MODE_CONTINUOUS,
                    )
                    self.qubit.singlepass0_on()
                with parallel:
                    # this must go after set_mu because set_mu calls io_update
                    delay_mu(time_gate_dd_mu)
                    self.dds_ramper_ms.configure_falling_ramp_all_dds()

                '''FALLING EDGE OF PULSE'''
                self.qubit.singlepass0_off()
                self.dds_ramper_ms.run_ramp_all_dds()
                delay_mu(self.dds_ramper_ms.drg_time_ramp_mu[0])
                self.dds_ramper_ms.switch_off_all_dds()
                self.qubit.off()
                self.ttl8.off()

                '''ENSURE ALL CFR REGISTERS ARE RESET TO NORMAL VALUES FOR SINGLE TONE OPERATION'''
                self.dds_ramper_ms.reset_cfrs_all_dds()

                """
                RUN PARITY PULSE
                """
                if self.enable_parity_pulse:
                    self.qubit.cpld.set_all_att_mu(self.att_reg_parity)
                    self.qubit.set_profile(self.profile_729_parity)
                    at_mu((now_mu() + 8) & ~7)
                    self.qubit.io_update()
                    # delay_mu(2000)

                    self.qubit.on()
                    self.qubit.singlepass0_on()

                    delay_mu(time_parity_pulse_mu)

                    self.qubit.off()
                    self.qubit.singlepass0_off()

                # # read out fluorescence & clean up loop
                self.readout_subsequence.run_dma()
                counts_res = self.readout_subsequence.fetch_count()

                self.update_results(
                    time_gate_ms_mu,
                    counts_res,
                    phase_parity_pulse_pow,
                    freq_secular_detuning_ms_ftw,
                    phase_dynamical_decoupling_pow,
                    freq_carrier_detuning_ms_ftw << 1,
                    time_parity_pulse_mu
                )

                # self.dds_ramper_ms.sequence_cleanup()
                self.check_termination()

    @kernel(flags={'fast_math'})
    def setup_ms_tones(self, freq_carrier_detuning_ms_ftw: TInt32,
                       freq_secular_detuning_ftw: TInt32,
                       phase_dynamical_decoupling_pow: TInt32) -> TNone:
        # ensure all beams are off
        self.qubit.off()
        self.qubit.singlepass0_off()
        self.qubit.singlepass1_off()
        self.qubit.singlepass2_off()

        self.qubit.set_mu(
            self.freq_carrier_ftw + freq_carrier_detuning_ms_ftw,  asf=self.qubit.ampl_qubit_asf,
            pow_=0, profile=self.profile_729_ms,
            phase_mode=ad9910.PHASE_MODE_CONTINUOUS
        )
        self.qubit.singlepass0.set_mu(
            self.qubit.freq_singlepass0_default_ftw,
            asf=self.ampl_parity_pulse_asf,
            pow_=phase_dynamical_decoupling_pow,
            profile=self.profile_729_ms,
            phase_mode=ad9910.PHASE_MODE_CONTINUOUS,
        )
        self.qubit.singlepass1.set_mu(
            self.qubit.freq_singlepass1_default_ftw - self.freq_secular_ftw - freq_secular_detuning_ftw,
            asf=self.ampls_ms_asf[0],
            pow_=0,
            profile=self.profile_729_ms,
            phase_mode=ad9910.PHASE_MODE_CONTINUOUS,
        )
        self.qubit.singlepass2.set_mu(
            self.qubit.freq_singlepass2_default_ftw + self.freq_secular_ftw + freq_secular_detuning_ftw,
            asf=self.ampls_ms_asf[1],
            pow_=0,
            profile=self.profile_729_ms,
            phase_mode=ad9910.PHASE_MODE_CONTINUOUS,
        )

    @kernel(flags={'fast_math'})
    def setup_parity_tones(self, phase_parity_pulse_pow: TInt32) -> TNone:
        # ensure all beams are off
        self.qubit.off()
        self.qubit.singlepass0_off()
        self.qubit.singlepass1_off()
        self.qubit.singlepass2_off()

        self.qubit.set_mu(
            self.freq_carrier_ftw, asf=self.qubit.ampl_qubit_asf,
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
