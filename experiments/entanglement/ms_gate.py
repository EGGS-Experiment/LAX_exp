from artiq.experiment import *
from artiq.coredevice import ad9910

from numpy import array, int32, int64, zeros

from LAX_exp.language import *
from LAX_exp.system.subsequences import (
    InitializeQubit, SidebandCoolContinuousRAM, Readout, ReadoutAdaptive,
    RescueIon, QubitRAP
)


class MolmerSorensen(LAXExperiment, Experiment):
    """
    Experiment: Molmer Sorensen Gate


    """
    name = 'Molmer Sorensen Gtae'

    kernel_invariants = {}

    def build_experiment(self):

        # allocate relevant beam profiles
        self.profile_729_RAP = 0
        self.profile_729_SBC = 1
        self.profile_729_readout = 2
        self.profile_729_ms = 3
        self.profile_729_parity = 4

        # devices
        self.setattr_device('qubit')
        self.setattr_device('pump')
        self.setattr_device('repump_qubit')

        self.setattr_argument('reps', NumberValue(default=50, min=10, max=1e6,
                                                  step = 1, scale=1, precision=0), tooltip='Number of shots to run')

        self._build_arguements_ms()
        self._build_arguements_parity()


        '''
        Build subsequences
        '''
        self.sidebandcool_subsequence = SidebandCoolContinuousRAM(self, profile_729=self.profile_729_SBC,
        profile_854=3, ram_addr_start_729=0, ram_addr_start_854=0m num_samples=200)
        self.initialize_subsequence = InitializeQubit(self)
        self.readout_subsequence = Readout(self)
        self.rescue_subsequence = RescueIon(self)



    def _build_arguements_ms(self):

        _name = 'ms'

        self.setattr_argument('time_gate_us_list', Scannable(default=[ExplicitScan([50]),
                                                                RangeScan(0.01, 100, 100)],
                                                                global_min=0.01, global_max=1e6,
                                                                global_step = 0.01, unit="us",
                                                                precision=4, scale=1),
                              group= _name, tooltip='time to apply birchromatic interaction for'
                              )


        self.setattr_argument('freq_carrier_ms_mhz',NumberValue(default=101.3,
                                                                min=90, max=101,
                                                                step = 0.01, unit="MHz",
                                                                precision=4, scale=1),
                              group = _name,tooltip='carrier frequency in MHz'
                              )


        self.setattr_argument('detuning_ms_khz_list', Scannable(default=[ExplicitScan([40]),
                                                                RangeScan(0.01, 100, 100)],
                                                                global_min=0.01, global_max=1e6,
                                                                global_step = 0.01, unit="kHz",
                                                                precision=4, scale=1),
                              group = _name,tooltip='detuning from carrier'
                              )

        self.setattr_argument("ampls_ms_pct", PYONValue([50., 50.]), group='default.ms',
                              tooltip="DDS amplitudes for the singlepass DDSs during the bichromatic pulses.\n"
                                      "Should be a list of [rsb_ampl_pct, bsb_ampl_pct], which are applied to [singlepass1, singlepass2].")
        self.setattr_argument("atts_ms_db", PYONValue([11., 11.]), group='default.ms',
                              tooltip="DDS attenuations for the singlepass DDSs during the bichromatic pulses.\n"
                                      "Should be a list of [rsb_att_db, bsb_att_db], which are applied to [singlepass1, singlepass2].")


    def _build_arguements_parity(self):

        _name = '_parity_pulse'

        self.setattr_argument('enable_parity_pulse', BooleanValue(False),
                              tooltip='Enable parity analysis pulse after applying the MS gate',
                              group=_name)

        self.setattr_argument('time_parity_pulse_us',
                              NumberValue(default=2, min=0.01, max=1e6,
                                          step=0.01, scale=1, ndecimals=4, unit='us'),
                              tooltip='', group=_name)

        self.setattr_argument('phase_parity_pulse_turns_list',
                            Scannable(default=[
                                ExplicitScan([0.]),
                                RangeScan(0, 1, 11),
                                CenterScan(0.5, 1, 0.1)], unit='turns',
                                global_min=0., global_max=2., global_step=0.1, precision=3, scale=1.0),
                            tooltip = '',
                            group = _name)

        self.setattr_argument('att_parity_pulse_db', NumberValue(default=5., min=5., max=31.5,
                                                                step=0.5, ndecimals=1, scale=1, unit='dB'),
                                                                tooltip='', group=_name)

        self.setattr_argument('ampl_parity_pulse_pct', NumberValue(default=50., min=0.01, max=50.,
                                                                step=0.01, ndecimals=3, scale=1, unit='%'),
                                                                tooltip='', group=_name)




    def prepare_experiment(self):
        '''
        Prepare all necessary functions
        '''

        '''
        Convert MS arguments
        '''
        time_gate_mu_list = [self.core.seconds_to_mu(time_gate_us*us) for time_gate_us in
                                  self.time_gate_us_list]
        self.freq_carrier_ms_mu = self.qubit.frequency_to_ftw(self.freq_carrier_ms_MHz * MHz)
        freq_detuning_ms_mu_list = [self.qubit.frequency_to_ftw(detuning_khz*kHz) for detuning_kHz in
                                    self.freq_detuning_ms_khz_list]
        self.ampls_ms_asf = array([self.qubit.amplitude_to_asf(ampl_pct / 100.)
                                    for ampl_pct in self.ampls_ms_pct])
        self.ampl_qubit_ms_asf = self.qubit.amplitude_to_asf(8.0)


        '''
        Convert Parity Arguments
        '''

        self.time_parity_pulse_mu = self.core.seconds_to_mu(self.time_parity_pulse_us*us)
        phase_parity_pulse_pow_list = [self.qubit.turns_to_pow(phase_parity_pulse)
                                            for phase_parity_pulse in self.phase_parity_pulse_turns_list]
        self.ampl_parity_pulse_asf = self.qubit.amplitude_to_asf(self.ampl_parity_pulse_pct/100.)



        '''
        Prepare Attenuation Register
        '''

        # attenuation register - bichromatic
        self.att_reg_ms = 0x00000000 | (
                (att_to_mu(self.att_doublepass_default_db * dB) << ((self.qubit.beam.chip_select - 4) * 8)) |
                (att_to_mu(31.5*dB) << ((self.qubit.singlepass0.chip_select - 4) * 8)) |
                (att_to_mu(self.atts_ms_db[0] * dB) << ((self.qubit.singlepass1.chip_select - 4) * 8)) |
                (att_to_mu(self.atts_ms_db[1] * dB) << ((self.qubit.singlepass2.chip_select - 4) * 8))
        )

        # attenuation register - parity
        self.att_reg_parity = 0x00000000 | (
                (att_to_mu(self.att_doublepass_default_db * dB) << ((self.qubit.beam.chip_select - 4) * 8)) |
                (att_to_mu(self.att_parity_pulse_db) << ((self.qubit.singlepass0.chip_select - 4) * 8)) |
                (att_to_mu(31.5*dB) << ((self.qubit.singlepass1.chip_select - 4) * 8)) |
                (att_to_mu(31.5 * dB) << ((self.qubit.singlepass2.chip_select - 4) * 8))
        )

        self.config_experiment_list = create_experiment_config(
           time_gate_mu_list,
            freq_detuning_ms_mu_list,
            phase_parity_pulse_pow_list,

            config_type=float, shuffle_config=True
        )

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_experiment_list),
                3)


    @kernel(flags={'fast_math'})
    def initialize_experiment(self) -> TNone:

        self.initialize_subsequence.record_dma()
        self.sidebandcool_subsequence.record_dma()
        self.readout_subsequence.record_dma()
        self.core.break_realtime()

    @kernel(flags={'fast_math'})
    def run_main(self) -> TNone:
        for trial_num in range(self.repetitions):
            for config_vals in self.config_experiment_list:
                time_gate_ms_mu = int64(config_vals[0])
                freq_detuning_ms_mu = int32(config_vals[1])
                phase_parity_pulse_pow = int32(config_vals[2])

                self.core.break_realtime()  # add slack for execution
                delay_mu(125000)  # add even more slack lol

                self.setup_ms_tones()
                self.setup_parity_tones()

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
                self.qubit.cpld.set_all_att_mu(self.att_reg_ms)
                # self.qubit.singlepass0_off()
                self.setup_beam_profiles()
                self.qubit.set_cfr1(phase_autoclear=1)
                self.qubit.singlepass0.set_cfr1(phase_autoclear=1)
                self.qubit.singlepass1.set_cfr1(phase_autoclear=1)
                self.qubit.singlepass2.set_cfr1(phase_autoclear=1)

                # pulse to clear phase accumulator
                self.qubit.io_update()
                # reset cfr1 so we no longer clear phases on io_update
                self.qubit.set_cfr1()
                self.qubit.singlepass0.set_cfr1()
                self.qubit.singlepass1.set_cfr1()
                self.qubit.singlepass2.set_cfr1()

                self.qubit.io_update()

                """
                RUN MS GATE
                """
                self.qubit.cpld.set_all_att_mu(self.att_reg_ms)
                self.qubit.set_profile(self.profile_729_ms)
                at_mu((now_mu() + 8) & ~7)
                self.qubit.io_update()
                delay_mu(2000)

                with parallel:
                    self.qubit.on()
                    self.qubit.singlepass1_on()
                    self.qubit.singlepass2_on()

                delay_mu(time_gate_ms_mu)

                with parallel:
                    self.qubit.off()
                    self.qubit.singlepass1_off()
                    self.qubit.singlepass2_off()


                """
                RUN PARITY PULSE
                """
                self.qubit.cpld.set_all_att_mu(self.att_reg_parity)
                self.qubit.set_profile(self.profile_729_parity)
                at_mu((now_mu() + 8) & ~7)
                self.qubit.io_update()
                delay_mu(2000)

                with parallel:
                    self.qubit.on()
                    self.qubit.singlepass0_on()

                delay_mu(self.time_parity_pulse_mu)

                with parallel:
                    self.qubit.off()
                    self.qubit.singlepass1_off()
                    self.qubit.singlepass2_off()

                # read out fluorescence & clean up loop
                self.readout_subsequence.run_dma()
                counts_res = self.readout_subsequence.fetch_count()

                self.update_results(
                    time_gate_ms_mu,
                    counts_res,
                    phase_parity_pulse_pow,
                    freq_detuning_ms_mu
                )




    @kernel(flags={'fast_math'})
    def setup_ms_tones(self, freq_detuning_mu):
        # ensure all beams are off
        self.qubit.off()
        self.qubit.singlepass0_off()
        self.qubit.singlepass1_off()
        self.qubit.singlepass2_off()

        self.qubit.set_mu(
            self.freq_carrier_ms_mu, asf=self.ampl_qubit_ms_asf,
            pow_=0, profile=self.profile_729_ms,
            phase_mode=ad9910.PHASE_MODE_CONTINUOUS
        )
        self.qubit.singlepass1.set_mu(
            self.qubit.singlepass1.freq_singlepass1_default_ftw - freq_detuning_mu,
            asf=self.ampl_qubit_ms_asf[0],
            pow_=0,
            profile=self.profile_729_ms,
            phase_mode=ad9910.PHASE_MODE_CONTINUOUS,
        )
        self.qubit.singlepass2.set_mu(
            self.qubit.singlepass2.freq_singlepass2_default_ftw + freq_detuning_mu,
            asf=self.ampl_qubit_ms_asf[1],
            pow_=0,
            profile=self.profile_729_ms,
            phase_mode=ad9910.PHASE_MODE_CONTINUOUS,
        )

    def setup_parity_tones(self, phase_parity_pulse_pow):
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
            self.qubit.singlepass0,
            asf=self.ampl_parity_pulse_asf,
            pow_=phase_parity_pulse_pow,
            profile=self.profile_729_parity,
            phase_mode=ad9910.PHASE_MODE_CONTINUOUS,
        )


    def analyze_experiment(self) -> TNone:
        pass




