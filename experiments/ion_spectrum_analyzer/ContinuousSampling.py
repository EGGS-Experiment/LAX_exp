import numpy as np
from artiq.experiment import *

from LAX_exp.language import *
from LAX_exp.system.subsequences import (
    InitializeQubit, Readout, RescueIon, SidebandCoolContinuousRAM, QubitRAP
)

from LAX_exp.system.objects.SpinEchoWizardRDX import SpinEchoWizardRDX
from LAX_exp.system.objects.PhaserPulseShaper import PhaserPulseShaper

# todo: check kernel invariants
# todo: double check all arg names


class ContinuousSampling(LAXExperiment, Experiment):
    """
    Experiment: Continuous Sampling

    todo: document
    """
    name = 'Continuous Sampling'
    kernel_invariants = {
        # hardware values
        'freq_osc_base_hz_list', 'freq_phaser_duc_hz', 'freq_global_offset_hz', 'att_phaser_mu', 'pulseshaper_vals',
        'sample_period_mu',

        # subsequences
        'initialize_subsequence', 'sidebandcool_subsequence', 'readout_subsequence',
        'rescue_subsequence', 'rap_subsequence',

        # configs
        'profile_729_SBC', 'profile_729_RAP', '_num_phaser_oscs'
    }

    def build_experiment(self):
        # exp-specific variables
        _argstr = "CS"              # string to use for arguments
        self._num_phaser_oscs = 5   # number of phaser oscillators in use

        # core arguments
        self.setattr_argument("num_samples",        NumberValue(default=4, precision=0, step=1, min=1, max=100000))
        self.setattr_argument("sample_period_ms",   NumberValue(default=10., precision=6, min=12, max=1e5, step=1, unit="ms", scale=1.))

        # allocate relevant beam profiles
        self.profile_729_SBC = 1
        self.profile_729_RAP = 2

        # get subsequences
        self.sidebandcool_subsequence =     SidebandCoolContinuousRAM(
            self, profile_729=self.profile_729_SBC, profile_854=3,
            ram_addr_start_729=0, ram_addr_start_854=0,
            num_samples=200
        )
        self.initialize_subsequence =       InitializeQubit(self)
        self.readout_subsequence =          Readout(self)
        self.rescue_subsequence =           RescueIon(self)

        # readout - RAP
        self.setattr_argument("att_rap_db",             NumberValue(default=8, precision=1, step=0.5, min=8, max=31.5, unit="dB", scale=1.), group="RAP")
        self.setattr_argument("ampl_rap_pct",           NumberValue(default=50., precision=3, step=5, min=1, max=50, unit="%", scale=1.), group="RAP")
        self.setattr_argument("freq_rap_center_mhz",    NumberValue(default=101.3318, precision=6, step=1e-2, min=60, max=200, unit="MHz", scale=1.), group='RAP')
        self.setattr_argument("freq_rap_dev_khz",       NumberValue(default=100., precision=2, step=0.01, min=1, max=1e4, unit="kHz", scale=1.), group='RAP')
        self.setattr_argument("time_rap_us",            NumberValue(default=200., precision=3, min=1, max=1e5, step=1, unit="us", scale=1.), group="RAP")

        # waveform - global config
        self.setattr_argument("att_phaser_db",          NumberValue(default=27., precision=1, step=0.5, min=0, max=31.5, unit="dB", scale=1.), group="{}.global".format(_argstr))
        self.setattr_argument("freq_phaser_carrier_mhz", NumberValue(default=86, precision=7, step=1, min=0.001, max=4800, unit="MHz", scale=1.), group="{}.global".format(_argstr))
        self.setattr_argument("freq_global_offset_mhz", NumberValue(default=2., precision=6, step=1., min=-10., max=10., unit="MHz", scale=1.), group="{}.global".format(_argstr))
        self.setattr_argument("phase_phaser_ch1_global_turns",  NumberValue(default=0.1, precision=3, step=0.05, min=-1.1, max=1.1, unit="turns", scale=1.), group="{}.global".format(_argstr))
        self.setattr_argument("osc_num_target_list",    PYONValue([0, 1]), group="{}.global".format(_argstr))

        # waveform - custom specification
        self.setattr_argument("time_osc_pulse_us",      NumberValue(default=500, precision=2, step=500, min=0.04, max=100000000, unit="us", scale=1.),
                              group="{}.waveform".format(_argstr))
        self.setattr_argument("freq_osc_khz_list",      PYONValue([-702., 702., -0.5, 0., 0.]), group="{}.waveform".format(_argstr))
        self.setattr_argument("ampl_osc_frac_list",     PYONValue([40., 40., 1., 1., 0.]), group="{}.waveform".format(_argstr))
        self.setattr_argument("phase_osc_turns_list",   PYONValue([0., 0., 0., 0.5, 0.]), group="{}.waveform".format(_argstr))
        self.setattr_argument("phase_osc_ch1_offset_turns",     PYONValue([0., 0., 0.5, 0.5, 0.]), group="{}.waveform".format(_argstr))

        # waveform - pulse shaping
        self.setattr_argument("enable_pulse_shaping",   BooleanValue(default=False), group='{}.pulse_shaping'.format(_argstr))
        self.setattr_argument("type_pulse_shape",       EnumerationValue(['sine_squared', 'error_function', 'slepian'], default='sine_squared'),
                              group='{}.pulse_shaping'.format(_argstr))
        self.setattr_argument("time_pulse_shape_rolloff_us",    NumberValue(default=100, precision=1, step=100, min=0.2, max=100000, unit="us", scale=1.),
                              group='{}.pulse_shaping'.format(_argstr))
        self.setattr_argument("freq_pulse_shape_sample_khz",    NumberValue(default=1500, precision=0, step=100, min=1, max=5000, unit="kHz", scale=1.),
                              group='{}.pulse_shaping'.format(_argstr))

        # waveform - PSK (Phase-shift Keying)
        self.setattr_argument("enable_phase_shift_keying",  BooleanValue(default=True), group="{}.psk".format(_argstr))
        self.setattr_argument("num_psk_phase_shifts",       NumberValue(default=4, precision=0, step=10, min=1, max=200), group="{}.psk".format(_argstr))
        self.setattr_argument("enable_psk_delay",           BooleanValue(default=True), group="{}.psk".format(_argstr))
        self.setattr_argument("time_psk_delay_us",          NumberValue(default=200., precision=3, min=1, max=1e5, step=1, unit="us", scale=1.), group="{}.psk".format(_argstr))
        self.setattr_argument("phase_osc0_psk_turns",       PYONValue([0., 0.5, 0., 0.5, 0.]), group="{}.psk".format(_argstr))
        self.setattr_argument("phase_osc1_psk_turns",       PYONValue([0., 0.5, 0., 0.5, 0.]), group="{}.psk".format(_argstr))
        self.setattr_argument("phase_osc2_psk_turns",       PYONValue([0., 0., 0., 0., 0.]), group="{}.psk".format(_argstr))
        self.setattr_argument("phase_osc3_psk_turns",       PYONValue([0., 0., 0., 0., 0.]), group="{}.psk".format(_argstr))
        self.setattr_argument("phase_osc4_psk_turns",       PYONValue([0., 0., 0., 0., 0.]), group="{}.psk".format(_argstr))

        # get relevant devices
        self.setattr_device("qubit")
        self.setattr_device('phaser_eggs')

        # instantiate RAP here since it relies on experiment arguments
        self.rap_subsequence = QubitRAP(
            self, ram_profile=self.profile_729_RAP, ram_addr_start=202, num_samples=250,
            ampl_max_pct=self.ampl_rap_pct, pulse_shape="blackman"
        )

        # instantiate helper objects
        self.spinecho_wizard = SpinEchoWizardRDX(self)

    def prepare_experiment(self):
        """
        Prepare experimental values.
        """
        '''SANITIZE & VALIDATE INPUTS'''
        self._prepare_argument_checks()

        '''SUBSEQUENCE PARAMETERS'''
        # set correct phase delays for field geometries (0.5 for osc_2 for dipole)
        # note: sequence blocks are stored as [block_num, osc_num] and hold [ampl_pct, phase_turns]
        # e.g. self.sequence_blocks[2, 5, 0] gives ampl_pct of 5th osc in 2nd block
        # note: create object here instead of build since phase_osc_ch1_offset_turns isn't well-defined until prepare
        self.pulse_shaper = PhaserPulseShaper(self, np.array(self.phase_osc_ch1_offset_turns))

        # prepare RAP arguments
        self.att_rap_mu = att_to_mu(self.att_rap_db * dB)
        self.freq_rap_center_ftw = hz_to_ftw(self.freq_rap_center_mhz * MHz)
        self.freq_rap_dev_ftw = hz_to_ftw(self.freq_rap_dev_khz * kHz)
        self.time_rap_mu = self.core.seconds_to_mu(self.time_rap_us * us)

        '''HARDWARE VALUES - CONFIG'''
        # convert argument units to whatever is most convenient for us
        self.sample_period_mu = self.core.seconds_to_mu(self.sample_period_ms * ms)
        self.att_phaser_mu = att_to_mu(self.att_phaser_db * dB)
        self.freq_global_offset_hz = self.freq_global_offset_mhz * MHz

        # convert build arguments to appropriate values and format as numpy arrays
        self.freq_phaser_duc_hz = self.freq_phaser_carrier_mhz * MHz - self.phaser_eggs.freq_center_hz - self.freq_global_offset_hz
        self.freq_osc_base_hz_list = np.array(self.freq_osc_khz_list) * kHz + self.freq_global_offset_hz

        # configure waveform via pulse shaper & spin echo wizard
        self._prepare_waveform()

    def _prepare_argument_checks(self) -> TNone:
        """
        Check experiment arguments for validity.
        """
        # check that input amplitude/phase arrays are valid
        if isinstance(self.ampl_osc_frac_list, list):
            if len(self.ampl_osc_frac_list) != self._num_phaser_oscs:
                raise ValueError("Error: phaser oscillator amplitude array must have length {:d}.".format(self._num_phaser_oscs))
            elif np.sum(self.ampl_osc_frac_list) >= 100.:
                raise ValueError("Error: phaser oscillator amplitudes must sum <100.")
        else:
            raise ValueError("Error: phaser oscillator amplitude array must be a list.")

        if isinstance(self.phase_osc_turns_list, list):
            if len(self.phase_osc_turns_list) != self._num_phaser_oscs:
                raise ValueError("Error: phaser oscillator phase array must have length {:d}.".format(self._num_phaser_oscs))
        else:
            raise ValueError("Error: phaser oscillator phase array must be a list.")

        # check that phaser oscillator frequencies are valid
        if not isinstance(self.freq_osc_khz_list, list):
            raise ValueError("Error: phaser oscillator frequency array must be a list.")
        elif len(self.freq_osc_khz_list) != self._num_phaser_oscs:
            raise ValueError("Error: phaser oscillator frequency array must have length {:d}.".format(self._num_phaser_oscs))
        max_osc_freq_hz = (
                max(self.freq_osc_khz_list) * kHz +
                (self.freq_global_offset_mhz * MHz)
        )
        min_osc_freq_hz = (
                max(self.freq_osc_khz_list) * kHz +
                (self.freq_global_offset_mhz * MHz)
        )
        if (max_osc_freq_hz > 10. * MHz) or (min_osc_freq_hz < -10. * MHz):
            raise ValueError("Error: phaser oscillator frequencies outside valid range of [-10, 10] MHz.")

        # ensure phaser output frequency falls within valid DUC bandwidth
        phaser_carrier_lower_dev_hz = abs(self.phaser_eggs.freq_center_hz - self.freq_phaser_carrier_mhz * MHz)
        phaser_carrier_upper_dev_hz = abs(self.phaser_eggs.freq_center_hz - self.freq_phaser_carrier_mhz * MHz)
        if (phaser_carrier_upper_dev_hz >= 200. * MHz) or (phaser_carrier_lower_dev_hz >= 200. * MHz):
            raise ValueError("Error: output frequencies outside +/- 300 MHz phaser DUC bandwidth.")

        # check that PSK schedule is valid
        psk_schedule_invalid = self.enable_phase_shift_keying and any([
            (not isinstance(psk_schedule, list)) or (len(psk_schedule) != self.num_psk_phase_shifts + 1)
            for psk_schedule in (
                self.phase_osc0_psk_turns, self.phase_osc1_psk_turns, self.phase_osc2_psk_turns,
                self.phase_osc3_psk_turns, self.phase_osc4_psk_turns
            )
        ])
        if psk_schedule_invalid:
            raise ValueError("Invalid PSK schedule. Must be a list of length num_psk_phase_shifts+1.")

        # ensure that the spinecho-ing makes sense
        if self.enable_psk_delay and not self.enable_phase_shift_keying:
            raise ValueError("Invalid waveform configuration. Cannot have delays enabled without PSKing.")

        # ensure that osc_num_target_list contains a valid selection of oscillators
        if not (
                all(isinstance(val, int) for val in self.osc_num_target_list) and
                (max(self.osc_num_target_list) <= 4 and min(self.osc_num_target_list) >= 0) and
                len(self.osc_num_target_list) <= 4
        ):
            raise ValueError("Invalid target oscillator list. Must be a list of fewer than 4 numbers in [0, 4].")

    def _prepare_waveform(self) -> TNone:
        """
        Calculate waveforms and timings for the EGGS pulse.
        Uses SpinEchoWizardRDX and PhaserPulseShaper objects to simplify waveform compilation.
        """
        '''PREPARE WAVEFORM COMPILATION'''
        # create holding structures for EGGS pulse waveforms
        self.pulseshaper_vals = None        # store compiled waveforms from pulseshaper
        self.pulseshaper_id =   np.int32(0) # store waveform ID for pulseshaper

        # calculate block timings and scale amplitudes for ramsey-ing
        if self.enable_phase_shift_keying:
            if self.enable_psk_delay:
                num_blocks = 2 * self.num_psk_phase_shifts + 1
                block_time_list_us = riffle([self.time_osc_pulse_us] * (self.num_psk_phase_shifts + 1),
                                            [self.time_psk_delay_us] * self.num_psk_phase_shifts)
                block_ampl_scale_list = riffle([1] * (self.num_psk_phase_shifts + 1), [0] * self.num_psk_phase_shifts)
            else:
                num_blocks = self.num_psk_phase_shifts + 1
                block_time_list_us = [self.time_osc_pulse_us] * (self.num_psk_phase_shifts + 1)
                block_ampl_scale_list = [1] * (self.num_psk_phase_shifts + 1)
        else:
            num_blocks = 1
            block_time_list_us = [self.time_osc_pulse_us]
            block_ampl_scale_list = [1] * (self.num_psk_phase_shifts + 1)

        '''PROGRAM & COMPILE WAVEFORM'''
        # create bare waveform block sequence & set amplitudes
        _osc_vals_blocks = np.zeros((num_blocks, self._num_phaser_oscs, 2), dtype=float)
        _osc_vals_blocks[:, :, 0] = np.array(self.ampl_osc_frac_list)
        _osc_vals_blocks[:, :, 0] *= np.array([block_ampl_scale_list]).transpose()

        # set oscillator phases and account for oscillator update delays
        # note: use mean of osc freqs since I don't want to record a waveform for each osc freq
        t_update_delay_s_list = np.array([0, 40e-9, 80e-9, 80e-9, 120e-9])[:self._num_phaser_oscs]
        _osc_vals_blocks[:, :, 1] += (np.array(self.phase_osc_turns_list) +
                                      self.freq_osc_base_hz_list * t_update_delay_s_list)

        # set PSK phase update schedule
        if self.enable_phase_shift_keying:
            if self.enable_psk_delay:
                # note: use ::2 since we only update to non-delay blocks
                _osc_vals_blocks[::2, :, 1] += np.array([
                    self.phase_osc0_psk_turns, self.phase_osc1_psk_turns, self.phase_osc2_psk_turns,
                    self.phase_osc3_psk_turns, self.phase_osc4_psk_turns
                ][:self._num_phaser_oscs]).transpose()
            else:
                _osc_vals_blocks[:, :, 1] += np.array([
                    self.phase_osc0_psk_turns, self.phase_osc1_psk_turns, self.phase_osc2_psk_turns,
                    self.phase_osc3_psk_turns, self.phase_osc4_psk_turns
                ][:self._num_phaser_oscs]).transpose()

        # specify sequence as a dict of blocks, where each block is a dict
        _sequence_blocks = [
            {
                "oscillator_parameters": _osc_vals_blocks[i],
                "config": {
                    "time_us": block_time_list_us[i],
                    "pulse_shaping": self.enable_pulse_shaping,
                    "pulse_shaping_config": {
                        "pulse_shape":          self.type_pulse_shape,
                        "pulse_shape_rising":   self.enable_pulse_shaping,
                        "pulse_shape_falling":  self.enable_pulse_shaping,
                        "sample_rate_khz":      self.freq_pulse_shape_sample_khz,
                        "rolloff_time_us":      self.time_pulse_shape_rolloff_us
                    }
                }
            } for i in range(num_blocks)
        ]
        # create QVSA waveform and store data in a holder
        self.pulseshaper_vals = self.spinecho_wizard.compile_waveform(_sequence_blocks)

    @property
    def results_shape(self):
        return (self.num_samples, 3)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self) -> TNone:
        # record general subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.sidebandcool_subsequence.record_dma()
        self.readout_subsequence.record_dma()
        self.core.break_realtime()

        # configure RAP pulse
        self.rap_subsequence.configure(self.time_rap_mu, self.freq_rap_center_ftw, self.freq_rap_dev_ftw)
        delay_mu(25000)

        ### PHASER INITIALIZATION ###
        self.phaser_record()

        # set maximum attenuations for phaser outputs to prevent leakage
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.channel[0].set_att_mu(0x00)
        delay_mu(self.phaser_eggs.t_sample_mu)
        self.phaser_eggs.channel[1].set_att_mu(0x00)

        # configure global phaser configs (e.g. DUC)
        self.phaser_eggs.frequency_configure(
            self.freq_phaser_duc_hz,    # carrier frequency (via DUC)
            # oscillator frequencies
            [self.freq_osc_base_hz_list[0], self.freq_osc_base_hz_list[1], self.freq_osc_base_hz_list[2],
             self.freq_osc_base_hz_list[3], self.freq_osc_base_hz_list[4]],
            self.phase_phaser_ch1_global_turns  # global CH1 phase
        )

    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        # load waveform DMA handles
        self.pulse_shaper.waveform_load()
        delay_mu(500000)

        # other setup
        _loop_iter = 0  # used to check_termination more frequently
        _time_start_mu = now_mu()   # ensure samples are evenly spaced

        # MAIN LOOP
        for sample_num in range(self.num_samples):
            # ensure results are samples are evenly/deterministically spaced
            at_mu(_time_start_mu + (sample_num + 1) * self.sample_period_mu)

            '''STATE PREPARATION'''
            self.initialize_subsequence.run_dma()
            self.sidebandcool_subsequence.run_dma()

            '''WAVEFORM SEQUENCE'''
            self.phaser_run(self.pulseshaper_id)

            '''READOUT (RAP-BASED)'''
            self.qubit.set_att_mu(self.att_rap_mu)
            self.rap_subsequence.run_rap(self.time_rap_mu)
            self.readout_subsequence.run_dma()
            _time_stop_mu = now_mu()    # record stop time
            counts = self.readout_subsequence.fetch_count()

            # update dataset
            self.update_results(
                _time_start_mu + (sample_num + 1) * self.sample_period_mu,
                counts,
                _time_stop_mu
            )
            self.core.break_realtime()

            '''LOOP CLEANUP'''
            # resuscitate ion & detect deaths
            self.rescue_subsequence.resuscitate()
            self.rescue_subsequence.detect_death(counts)

            # check termination more frequently in case reps are low
            if _loop_iter % 50 == 0:
                self.check_termination()
            _loop_iter += 1


    '''
    HELPER FUNCTIONS - PHASER
    '''
    @kernel(flags={"fast-math"})
    def phaser_run(self, waveform_id: TInt32) -> TNone:
        """
        Run the main EGGS pulse together with supporting functionality.
        Arguments:
            waveform_id     (TInt32)    : the ID of the waveform to run.
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
        self.phaser_stop_rdx()

    @kernel(flags={"fast-math"})
    def phaser_stop_rdx(self) -> TNone:
        """
        Stop oscillators on phaser without clearing certain phase accumulators.
        Set maximum attenuation to prevent output leakage.
        Can be used following an EGGS pulse to stop the phaser.
        """
        # set amplitudes of all oscillators to 0
        for i in range(5):
            with parallel:
                self.phaser_eggs.channel[0].oscillator[i].set_amplitude_phase(amplitude=0., phase=0., clr=0)
                self.phaser_eggs.channel[1].oscillator[i].set_amplitude_phase(amplitude=0., phase=0., clr=0)
                delay_mu(self.phaser_eggs.t_sample_mu)

        # clear only phase accumulators of target oscs
        for osc_num in self.osc_num_target_list:
            with parallel:
                self.phaser_eggs.channel[0].oscillator[osc_num].set_amplitude_phase(amplitude=0., phase=0., clr=1)
                self.phaser_eggs.channel[1].oscillator[osc_num].set_amplitude_phase(amplitude=0., phase=0., clr=1)
                delay_mu(self.phaser_eggs.t_sample_mu)

        # add delay for oscillator updates to account for pipeline latency
        delay_mu(2560)  # 8 frame periods
        # stop phaser amp switches & deactivate integrator hold
        with parallel:
            self.phaser_eggs.ch0_amp_sw.off()
            self.phaser_eggs.ch1_amp_sw.off()
            self.phaser_eggs.int_hold.off()

        # add delay time after EGGS pulse to allow RF servo to re-lock
        delay_mu(self.phaser_eggs.time_phaser_holdoff_mu)

        # switch off EGGS attenuators to prevent phaser leakage
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.channel[0].set_att_mu(0x00)
        delay_mu(self.phaser_eggs.t_frame_mu)
        self.phaser_eggs.channel[1].set_att_mu(0x00)

    @kernel(flags={"fast-math"})
    def phaser_record(self) -> TNone:
        """
        Set up core phaser functionality and record the pulse-shaped waveforms.
        Should be run during initialize_experiment.
        """
        # get waveform for given sweep phase
        _wav_data_ampl, _wav_data_phas, _wav_data_time = self.pulseshaper_vals

        # record phaser pulse sequence and save returned waveform ID
        delay_mu(1000000)  # 1ms
        self.pulseshaper_id = self.pulse_shaper.waveform_record(_wav_data_ampl, _wav_data_phas, _wav_data_time)
        self.core.break_realtime()

