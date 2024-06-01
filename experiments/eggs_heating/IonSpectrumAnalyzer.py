import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.system.subsequences import SqueezeConfigurable
import LAX_exp.experiments.eggs_heating.EGGSHeating as EGGSHeating


class IonSpectrumAnalyzer(EGGSHeating.EGGSHeating):
    """
    Experiment: Ion Spectrum Analyzer

    ***todo: redocument***
    Cool the ions to the ground state of motion via sideband cooling,
    then apply bichromatic heating tones, and try to read out the fluorescence.
    """
    name = 'Ion Spectrum Analyzer'


    def build_experiment(self):
        # run build_experiment for EGGS Heating
        super().build_experiment()

        # ISA - frequency offset
        self.setattr_argument("freq_ISA_sideband_offset_khz_list",          Scannable(
                                                                                default=[
                                                                                    ExplicitScan([0]),
                                                                                    CenterScan(0, 5, 0.5, randomize=True)
                                                                                ],
                                                                                global_min=-8000, global_max=8000, global_step=10,
                                                                                unit="kHz", scale=1, ndecimals=3
                                                                            ), group='EGGS_Heating.frequencies')
        # ISA - CH1 turns
        self.setattr_argument("phase_eggs_heating_ch1_turns_list",          Scannable(
                                                                                default=[
                                                                                    ExplicitScan([0.]),
                                                                                    RangeScan(0, 1.0, 9, randomize=True),
                                                                                ],
                                                                                global_min=0.0, global_max=1.0, global_step=1,
                                                                                unit="turns", scale=1, ndecimals=3
                                                                            ), group='EGGS_Heating.waveform.time_phase')

        # ISA - antisqueezing
        self.setattr_argument("enable_ISA_antisqueezing",                   BooleanValue(default=False), group='ISA.antisqueezing')
        self.setattr_argument("ampl_ISA_antisqueezing_rsb_pct",             NumberValue(default=40., ndecimals=2, step=10, min=0.0, max=99), group='ISA.antisqueezing')
        self.setattr_argument("ampl_ISA_antisqueezing_bsb_pct",             NumberValue(default=40., ndecimals=2, step=10, min=0.0, max=99), group='ISA.antisqueezing')
        self.setattr_argument("phase_ISA_antisqueezing_rsb_turns",          NumberValue(default=0.5, ndecimals=3, step=0.1, min=-1., max=1.), group='ISA.antisqueezing')
        self.setattr_argument("phase_ISA_antisqueezing_bsb_turns",          NumberValue(default=0.5, ndecimals=3, step=0.1, min=-1., max=1.), group='ISA.antisqueezing')

        self.setattr_argument("enable_ISA_antisqueezing_dipole",            BooleanValue(default=False), group='ISA.antisqueezing')
        self.setattr_argument("ampl_ISA_antisqueezing_dipole_pct",          NumberValue(default=20., ndecimals=2, step=10, min=0.0, max=99), group='ISA.antisqueezing')
        self.setattr_argument("phase_ISA_antisqueezing_dipole_turns",       NumberValue(default=0.5, ndecimals=3, step=0.1, min=-1., max=1.), group='ISA.antisqueezing')

        # ISA - extrinsic squeezing (parametric)
        self.setattr_argument("freq_squeeze_khz",                           NumberValue(default=1542.2, ndecimals=3, step=10, min=1, max=400000), group='squeeze_configurable')
        self.setattr_argument("phase_antisqueeze_turns_list",               Scannable(
                                                                                    default=[
                                                                                        ExplicitScan([0.]),
                                                                                        RangeScan(0, 1.0, 6, randomize=True)
                                                                                    ],
                                                                                    global_min=0.0, global_max=1.0, global_step=1,
                                                                                    unit="turns", scale=1, ndecimals=3
                                                                                ), group='squeeze_configurable')
        self.setattr_argument("time_squeeze_us",                            NumberValue(default=50., ndecimals=3, step=100, min=1, max=1000000), group='squeeze_configurable')
        self.squeeze_subsequence =                                          SqueezeConfigurable(self)

        # get relevant devices
        self.setattr_device('dds_parametric')
        self.setattr_device('ttl10')


    def prepare_experiment(self):
        """
        Prepare experimental values.
        """
        # run prepare_experiment from EGGS Heating
        super().prepare_experiment()


        '''ISA - PARAMETERS/UNIT CONVERSION'''
        self.freq_ISA_sideband_offset_hz_list =         np.array(list(self.freq_ISA_sideband_offset_khz_list)) * kHz
        self.phase_antisqueeze_pow_list =               np.array([self.dds_parametric.turns_to_pow(phase_turns)
                                                                          for phase_turns in self.phase_antisqueeze_turns_list])
        self.phase_eggs_heating_ch1_turns_list =        np.array(list(self.phase_eggs_heating_ch1_turns_list))


        '''ISA - MODIFY EXPERIMENTAL CONFIG'''
        # copy config_eggs_heating_list and insert ISA sweep columns
        num_columns =                                   np.shape(self.config_eggs_heating_list)[1] + 3
        self.config_ISA_list =                          np.stack(np.meshgrid(*self.config_eggs_heating_list.transpose(),
                                                                             self.freq_ISA_sideband_offset_hz_list,
                                                                             self.phase_antisqueeze_pow_list,
                                                                             self.phase_eggs_heating_ch1_turns_list),
                                                                 -1).reshape(-1, num_columns)
        # todo: ensure ISA antisqueezing works properly


        '''ISA - INTRINSIC & EXTRINSIC SQUEEZING'''
        # note: this needs to happen last in the prepare_experiment stage to prevent its configs from being overwritten
        self._prepare_squeezing()


    def _prepare_squeezing(self):
        """
        Calculate values for intrinsic & extrinsic squeezing.
        """
        '''
        EXTRINSIC/PARAMETRIC SQUEEZING
        '''
        # prepare extrinsic/parametric squeezing
        self.freq_squeeze_ftw = self.dds_parametric.frequency_to_ftw(self.freq_squeeze_khz * kHz)
        self.time_squeeze_mu =  self.core.seconds_to_mu(self.time_squeeze_us * us)

        '''
        INTRINSIC SQUEEZING
        '''
        # calculate ISA antisqueezing time as half of the eggs heating time
        self.time_ISA_antisqueeze_mu =              np.int64(round(self.time_eggs_heating_mu / 2))
        # ensure halfway time (for ISA antisqueezing) is a multiple of the phaser sample period
        if self.time_ISA_antisqueeze_mu % self.phaser_eggs.t_sample_mu:
            # round ISA antisqueezing time to a multiple of the phaser sample period
            t_sample_multiples =                    round(self.time_ISA_antisqueeze_mu / self.phaser_eggs.t_sample_mu + 0.5)
            self.time_ISA_antisqueeze_mu =          np.int64(t_sample_multiples * self.phaser_eggs.t_sample_mu)

        # prepare internal antisqueezing
        if self.enable_ISA_antisqueezing:           self.phaser_run = self.phaser_run_ISA_antisqueezing
        else:                                       self.phaser_run = self.phaser_run_nopsk

        # convert intrinsic ISA squeezing values to pct
        self.ampl_ISA_antisqueezing_rsb_frac =      self.ampl_ISA_antisqueezing_rsb_pct / 100.
        self.ampl_ISA_antisqueezing_bsb_frac =      self.ampl_ISA_antisqueezing_bsb_pct / 100.
        self.ampl_ISA_antisqueezing_dipole_frac =   self.ampl_ISA_antisqueezing_dipole_pct / 100. * float(self.enable_ISA_antisqueezing_dipole)

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_ISA_list),
                9)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def run_main(self):
        self.core.reset()

        # get custom sequence handles
        _handle_eggs_pulseshape_rise =      self.core_dma.get_handle('_PHASER_PULSESHAPE_RISE')
        _handle_eggs_pulseshape_fall =      self.core_dma.get_handle('_PHASER_PULSESHAPE_FALL')
        self.core.break_realtime()

        # used to check_termination more frequently
        _loop_iter = 0


        # MAIN LOOP
        for trial_num in range(self.repetitions):

            # sweep experiment configurations
            for config_vals in self.config_ISA_list:

                '''CONFIGURE'''
                # extract values from config list
                freq_readout_ftw =          np.int32(config_vals[0])
                carrier_freq_hz =           config_vals[1]
                sideband_freq_hz =          config_vals[2]
                ampl_rsb_frac =             config_vals[3]
                ampl_bsb_frac =             config_vals[4]
                ampl_dd_frac =              config_vals[5]
                phase_rsb_turns =           config_vals[6]
                time_readout_mu =           np.int64(config_vals[7])
                offset_freq_hz =            config_vals[8]
                phase_antisqueeze_pow =     np.int32(config_vals[9])
                phase_ch1_turns =           config_vals[10]
                self.core.break_realtime()

                # configure EGGS tones and set readout
                self.phaser_configure(carrier_freq_hz, sideband_freq_hz, offset_freq_hz, phase_rsb_turns, phase_ch1_turns)
                self.core.break_realtime()
                self.qubit.set_mu(freq_readout_ftw, asf=self.sidebandreadout_subsequence.ampl_sideband_readout_asf, profile=0)
                self.core.break_realtime()

                # configure squeezing/antisqueezing
                self.squeeze_subsequence.configure(self.freq_squeeze_ftw, phase_antisqueeze_pow, self.time_squeeze_mu)
                self.core.break_realtime()

                '''STATE PREPARATION'''
                # initialize ion in S-1/2 state
                self.initialize_subsequence.run_dma()
                # sideband cool
                self.sidebandcool_subsequence.run_dma()
                # squeeze ion
                self.squeeze_subsequence.squeeze()

                '''ION SPECTRUM ANALYZER'''
                # EGGS - START/SETUP
                # activate integrator hold
                self.ttl10.on()
                # set phaser attenuators
                at_mu(self.phaser_eggs.get_next_frame_mu())
                self.phaser_eggs.channel[0].set_att(self.att_eggs_heating_db * dB)
                delay_mu(self.phaser_eggs.t_sample_mu)
                self.phaser_eggs.channel[1].set_att(self.att_eggs_heating_db * dB)

                # reset DUC phase to start DUC deterministically
                self.phaser_eggs.reset_duc_phase()
                self.core_dma.playback_handle(_handle_eggs_pulseshape_rise)

                # EGGS - RUN
                self.phaser_run(ampl_rsb_frac, ampl_bsb_frac, ampl_dd_frac)

                # EGGS - STOP
                self.core_dma.playback_handle(_handle_eggs_pulseshape_fall)
                self.phaser_eggs.phaser_stop()
                # deactivate integrator hold
                self.ttl10.off()
                # add delay time after EGGS pulse to allow RF servo to re-lock
                delay_mu(self.time_rf_servo_holdoff_mu)

                '''READOUT'''
                # antisqueeze ion
                self.squeeze_subsequence.antisqueeze()
                # shelve ion in the D-5/2 state
                self.sidebandreadout_subsequence.run_time(time_readout_mu)
                self.readout_subsequence.run_dma()
                counts = self.readout_subsequence.fetch_count()

                # update dataset
                with parallel:
                    self.update_results(
                        freq_readout_ftw,
                        counts,
                        carrier_freq_hz,
                        sideband_freq_hz,
                        offset_freq_hz,
                        time_readout_mu,
                        phase_antisqueeze_pow,
                        phase_rsb_turns,
                        phase_ch1_turns
                    )
                    self.core.break_realtime()

                # resuscitate ion
                self.rescue_subsequence.resuscitate()

                # death detection
                self.rescue_subsequence.detect_death(counts)
                self.core.break_realtime()

                # check termination more frequently in case reps are low
                if (_loop_iter % 50) == 0:
                    self.check_termination()
                    self.core.break_realtime()
                _loop_iter += 1

            # rescue ion as needed
            self.rescue_subsequence.run(trial_num)

            # support graceful termination
            with parallel:
                self.check_termination()
                self.core.break_realtime()

        # CLEANUP
        self.core.break_realtime()
        self.phaser_eggs.reset_oscillators()


    '''
    HELPER FUNCTIONS - PHASER
    '''
    @kernel(flags={"fast-math"})
    def phaser_setup(self) -> TNone:
        """
        Set up core phaser functionality and record the pulse-shaped waveforms.
        Should be run during initialize_experiment.

        Note: need to redefined for ISA since our phaser_configure is different.
        """
        # todo: document better
        # get starting phase values for pulse shaping ### todo: document better
        carrier_freq_hz =   self.config_ISA_list[0, 1]
        sideband_freq_hz =  self.config_ISA_list[0, 2]
        offset_freq_hz =    self.config_ISA_list[0, 3]
        phase_rsb_turns =   self.config_ISA_list[0, 9]
        phase_ch1_turns =   self.config_ISA_list[0, 10]
        self.core.break_realtime()

        # configure EGGS tones and set readout frequency; also necessary to ensure phase delays are correctly set
        self.phaser_configure(carrier_freq_hz, sideband_freq_hz, offset_freq_hz, phase_rsb_turns, phase_ch1_turns)

        # record phaser rising pulse shape DMA sequence
        self.core.break_realtime()
        if self.enable_pulse_shaping:
            with self.core_dma.record('_PHASER_PULSESHAPE_RISE'):
                # set amplitude values at given time
                for ampl_val_list in self.ampl_pulse_shape_frac_list:
                    self.phaser_pulseshape_point(ampl_val_list[0], ampl_val_list[1], ampl_val_list[2])
                    delay_mu(self.time_pulse_shape_delay_mu)
        else:
            with self.core_dma.record('_PHASER_PULSESHAPE_RISE'):
                pass

        # record phaser falling pulse shape DMA sequence
        self.core.break_realtime()
        if self.enable_pulse_shaping:
            with self.core_dma.record('_PHASER_PULSESHAPE_FALL'):
                # set amplitude values at given time
                for ampl_val_list in self.ampl_pulse_shape_reverse_frac_list:
                    self.phaser_pulseshape_point(ampl_val_list[0], ampl_val_list[1], ampl_val_list[2])
                    delay_mu(self.time_pulse_shape_delay_mu)
        else:
            with self.core_dma.record('_PHASER_PULSESHAPE_FALL'):
                pass

    @kernel(flags={"fast-math"})
    def phaser_configure(self, carrier_freq_hz: TFloat, sideband_freq_hz: TFloat, offset_freq_hz: TFloat,
                         phase_rsb_turns: TFloat, phase_ch1_turns: TFloat) -> TNone:
        """
        Configure the tones on phaser for EGGS.
        Puts the same RSB and BSB on both channels, and sets a third oscillator to 0 Hz in case dynamical decoupling is used.

        Arguments:
            carrier_freq_hz         (float)     : the carrier frequency (in Hz).
            sideband_freq_hz        (float)     : the sideband frequency (in Hz).
            offset_freq_hz          (float)     : the offset frequency between the carrier and the sideband center.
            phase_rsb_turns         (float)     : the phase for the rsb tone (in turns)
            phase_ch1_turns         (float)     : the phase for phaser Channel 1 (relative to Channel 0).
        """
        '''
        CALCULATE PHASE DELAYS
        '''
        # calculate phase delays between CH0 and CH1, accounting for the relative CH1 latency
        self.phase_ch1_turns = phase_ch1_turns + (carrier_freq_hz * self.phaser_eggs.time_latency_ch1_system_ns * ns)

        # calculate phase delays for each oscillator to account for inherent update latencies and system latencies
        # oscillator 0 (RSB)
        self.phase_phaser_turns_arr[0, 0] = phase_rsb_turns
        self.phase_phaser_turns_arr[1, 0] = phase_rsb_turns
        # oscillator 1 (BSB)
        self.phase_phaser_turns_arr[0, 1] = ((sideband_freq_hz + offset_freq_hz) * self.phaser_eggs.t_sample_mu * ns) + self.phase_eggs_heating_bsb_turns
        self.phase_phaser_turns_arr[1, 1] = ((sideband_freq_hz + offset_freq_hz) * self.phaser_eggs.t_sample_mu * ns) + self.phase_eggs_heating_bsb_turns
        # oscillator 2 (carrier) (note: ch1 has 0.5 turns to put carrier in dipole config)
        self.phase_phaser_turns_arr[0, 2] = 0.
        self.phase_phaser_turns_arr[1, 2] = 0.5
        self.core.break_realtime()

        '''
        SET CARRIER FREQUENCY
        '''
        # set carrier offset frequency via the DUC
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.channel[0].set_duc_frequency(carrier_freq_hz - self.phaser_eggs.freq_center_hz)
        delay_mu(self.phaser_eggs.t_frame_mu)
        self.phaser_eggs.channel[1].set_duc_frequency(carrier_freq_hz - self.phaser_eggs.freq_center_hz)
        delay_mu(self.phaser_eggs.t_frame_mu)
        # strobe updates for both channels
        self.phaser_eggs.duc_stb()

        # set DUC phase delay to compensate for inter-channel latency
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.channel[1].set_duc_phase(self.phase_ch1_turns)
        # todo: do I need to add another get_next_frame_mu?
        self.phaser_eggs.duc_stb()

        '''
        SET OSCILLATOR (i.e. sideband) FREQUENCIES
        '''
        # synchronize to frame
        at_mu(self.phaser_eggs.get_next_frame_mu())
        # set oscillator 0 (RSB)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[0].set_frequency(-sideband_freq_hz + offset_freq_hz)
            self.phaser_eggs.channel[1].oscillator[0].set_frequency(-sideband_freq_hz + offset_freq_hz)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # set oscillator 1 (BSB)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[1].set_frequency(sideband_freq_hz + offset_freq_hz)
            self.phaser_eggs.channel[1].oscillator[1].set_frequency(sideband_freq_hz + offset_freq_hz)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # set oscillator 2 (carrier)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[2].set_frequency(0.)
            self.phaser_eggs.channel[1].oscillator[2].set_frequency(0.)
            delay_mu(self.phaser_eggs.t_sample_mu)

    @kernel(flags={"fast-math"})
    def phaser_run_ISA_antisqueezing(self, ampl_rsb_frac: TFloat, ampl_bsb_frac: TFloat, ampl_dd_frac: TFloat):
        """
        Activate phaser channel outputs for EGGS heating.
        Sets the same RSB, BSB, and dynamical decoupling amplitudes for both channels.
        Adjusts the BSB phase by a set value after a set time (should be calculated in _prepare_squeezing).
        Arguments:
            ampl_rsb_frac   (float) : the red sideband amplitude (as a decimal fraction).
            ampl_bsb_frac   (float) : the blue sideband amplitude (as a decimal fraction).
            ampl_dd_frac    (float) : the dynamical decoupling amplitude (as a decimal fraction).
        """
        '''
        ISA - INTRINSIC SQUEEZING
        '''
        # set oscillator 0 (RSB)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[0].set_amplitude_phase(amplitude=ampl_rsb_frac, phase=self.phase_phaser_turns_arr[0, 0], clr=0)
            self.phaser_eggs.channel[1].oscillator[0].set_amplitude_phase(amplitude=ampl_rsb_frac, phase=self.phase_phaser_turns_arr[1, 0], clr=0)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # set oscillator 1 (BSB)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[1].set_amplitude_phase(amplitude=ampl_bsb_frac, phase=self.phase_phaser_turns_arr[0, 1], clr=0)
            self.phaser_eggs.channel[1].oscillator[1].set_amplitude_phase(amplitude=ampl_bsb_frac, phase=self.phase_phaser_turns_arr[1, 1], clr=0)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # set oscillator 2 (carrier)
        with parallel:
            self.phaser_eggs.channel[0].oscillator[2].set_amplitude_phase(amplitude=ampl_dd_frac, phase=self.phase_phaser_turns_arr[0, 2], clr=0)
            self.phaser_eggs.channel[1].oscillator[2].set_amplitude_phase(amplitude=ampl_dd_frac, phase=self.phase_phaser_turns_arr[1, 2], clr=0)

        # heat for first half (before antisqueezing)
        delay_mu(self.time_ISA_antisqueeze_mu)

        '''
        ISA - INTRINSIC ANTISQUEEZING
        '''
        # adjust oscillator 0 (BSB) phase for antisqueezing
        with parallel:
            self.ttl8.on()
            self.phaser_eggs.channel[0].oscillator[0].set_amplitude_phase(amplitude=self.ampl_ISA_antisqueezing_rsb_frac,
                                                                          phase=self.phase_phaser_turns_arr[0, 0] + self.phase_ISA_antisqueezing_rsb_turns, clr=0)
            self.phaser_eggs.channel[1].oscillator[0].set_amplitude_phase(amplitude=self.ampl_ISA_antisqueezing_rsb_frac,
                                                                          phase=self.phase_phaser_turns_arr[1, 0] + self.phase_ISA_antisqueezing_rsb_turns, clr=0)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # adjust oscillator 1 (BSB) phase for antisqueezing
        with parallel:
            self.phaser_eggs.channel[0].oscillator[1].set_amplitude_phase(amplitude=self.ampl_ISA_antisqueezing_bsb_frac,
                                                                          phase=self.phase_phaser_turns_arr[0, 1] + self.phase_ISA_antisqueezing_bsb_turns, clr=0)
            self.phaser_eggs.channel[1].oscillator[1].set_amplitude_phase(amplitude=self.ampl_ISA_antisqueezing_bsb_frac,
                                                                          phase=self.phase_phaser_turns_arr[1, 1] + self.phase_ISA_antisqueezing_bsb_turns, clr=0)
            delay_mu(self.phaser_eggs.t_sample_mu)
        # turn off oscillator 2 (carrier) during antisqueezing
        with parallel:
            self.phaser_eggs.channel[0].oscillator[2].set_amplitude_phase(amplitude=self.ampl_ISA_antisqueezing_dipole_frac,
                                                                          phase=self.phase_phaser_turns_arr[0, 2] + self.phase_ISA_antisqueezing_dipole_turns, clr=0)
            self.phaser_eggs.channel[1].oscillator[2].set_amplitude_phase(amplitude=self.ampl_ISA_antisqueezing_dipole_frac,
                                                                          phase=self.phase_phaser_turns_arr[1, 2] + self.phase_ISA_antisqueezing_dipole_turns, clr=0)

        # heat for second half
        delay_mu(self.time_ISA_antisqueeze_mu)
        self.ttl8.off()

