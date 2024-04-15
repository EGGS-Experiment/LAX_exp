import numpy as np
from artiq.experiment import *


class PhaserSet(EnvExperiment):
    """
    Utility: Phaser Set

    Set the phaser output for fast user testing.
    """
    name = 'Phaser Set'

    def build(self):
        # global arguments
        self.setattr_argument("repetitions",                    NumberValue(default=5, ndecimals=0, step=1, min=1, max=10000))
        self.setattr_argument("cleanup",                        BooleanValue(default=False), group='global')
        self.setattr_argument("freq_carrier_mhz",               NumberValue(default=82, ndecimals=5, step=1, min=1, max=10000), group='global')
        self.setattr_argument("time_pulse_ms",                  NumberValue(default=0.25, ndecimals=5, step=1, min=0.000001, max=100000), group='global')
        # self.setattr_argument("clear_dac_phase_accumulator", BooleanValue(default=False), group='global')
        self.setattr_argument("synchronize_oscillator_updates", BooleanValue(default=True), group='global')

        # channel-specific arguments
        self.setattr_argument("att_ch0_db",                     NumberValue(default=8., ndecimals=1, step=0.5, min=0, max=31.5), group='channels')
        self.setattr_argument("att_ch1_db",                     NumberValue(default=8., ndecimals=1, step=0.5, min=0, max=31.5), group='channels')
        self.setattr_argument("time_latency_ch1_system_ns",     NumberValue(default=0.17, ndecimals=3, step=0.1, min=-1.0, max=1.0), group='channels')
        self.setattr_argument("apply_ch1_latency_on_oscillators", BooleanValue(default=False), group='channels')
        self.setattr_argument("phase_inherent_ch1_turns",       NumberValue(default=0., ndecimals=3, step=0.1, min=-1.0, max=1.0), group='channels')

        # oscillators - channel 0
        # note: oscillator waveforms are in the format [freq_khz, ampl_pct, phas_turns]
        self.setattr_argument("waveform_ch0_osc0",  PYONValue(default=[-767.5, 40., 0.475]), group='oscillators_ch0')
        self.setattr_argument("waveform_ch0_osc1",  PYONValue(default=[767.5, 40., 0.25]), group='oscillators_ch0')
        self.setattr_argument("waveform_ch0_osc2",  PYONValue(default=[0., 10., 0.]), group='oscillators_ch0')
        self.setattr_argument("waveform_ch0_osc3",  PYONValue(default=[0., 0., 0.]), group='oscillators_ch0')
        self.setattr_argument("waveform_ch0_osc4",  PYONValue(default=[0., 0., 0.]), group='oscillators_ch0')

        # oscillators - channel 1
        # note: oscillator waveforms are in the format [freq_khz, ampl_pct, phas_turns]
        self.setattr_argument("waveform_ch1_osc0",  PYONValue(default=[767.5, 40., 0.475]), group='oscillators_ch1')
        self.setattr_argument("waveform_ch1_osc1",  PYONValue(default=[-767.5, 40., 0.25]), group='oscillators_ch1')
        self.setattr_argument("waveform_ch1_osc2",  PYONValue(default=[0., 10., 0.5]), group='oscillators_ch1')
        self.setattr_argument("waveform_ch1_osc3",  PYONValue(default=[0., 0., 0.]), group='oscillators_ch1')
        self.setattr_argument("waveform_ch1_osc4",  PYONValue(default=[0., 0., 0.]), group='oscillators_ch1')

        # build constants
        self._build_constants()

    def _build_constants(self):
        """
        Instantiate fixed/constant objects/variables.
        """
        # get relevant devices
        self.setattr_device('core')
        self.setattr_device('core_dma')
        self.setattr_device('phaser0')
        self.setattr_device('ttl8')
        self.setattr_device('ttl9')

        # hardware values
        self.t_sample_mu =          np.int64(40)
        self.t_frame_mu =           np.int64(320)
        self.time_output_delay_mu = np.int64(1953)

        # preallocate DAC34H84 0x1F register config to speed up
        # writes when clearing the phase
        self.dac_register_1f =      np.int32(0)

        # preallocate phase holders
        self.phase_ch1_turns =      np.float(0)

        self.phase_ch0_osc0 =       np.float(0)
        self.phase_ch0_osc1 =       np.float(0)
        self.phase_ch0_osc2 =       np.float(0)
        self.phase_ch0_osc3 =       np.float(0)
        self.phase_ch0_osc4 =       np.float(0)

        self.phase_ch1_osc0 =       np.float(0)
        self.phase_ch1_osc1 =       np.float(0)
        self.phase_ch1_osc2 =       np.float(0)
        self.phase_ch1_osc3 =       np.float(0)
        self.phase_ch1_osc4 =       np.float(0)


    def prepare(self):
        # check oscillator arguments and prepare their values
        self._prepare_oscillators()

        ### EXPERIMENT ###
        self._iter_repetitions =    np.arange(self.repetitions)

        ### TIMING ###
        # todo: double check this is OK and reproduces normal behavior
        self.time_pulse_mu =                        self.core.seconds_to_mu(self.time_pulse_ms * ms)
        if self.time_pulse_mu % self.t_frame_mu:
            t_frame_multiples =                     round(self.time_pulse_mu / self.t_frame_mu + 0.5)
            self.time_pulse_mu =                    np.int64(self.t_frame_mu * t_frame_multiples)

        ### FREQUENCIES ###
        self.freq_carrier_hz =                      self.freq_carrier_mhz * MHz
        self.freq_center_hz =       self.get_dataset('eggs.freq_center_mhz') * MHz

    def _prepare_oscillators(self):
        """
        Extract oscillator parameters from argument PYONValues and
        check their validity.
        """
        # collate oscillator waveform values
        self.waveform_ch0_osc_list =    np.array([self.waveform_ch0_osc0,
                                                  self.waveform_ch0_osc1,
                                                  self.waveform_ch0_osc2,
                                                  self.waveform_ch0_osc3,
                                                  self.waveform_ch0_osc4])
        self.waveform_ch1_osc_list =    np.array([self.waveform_ch1_osc0,
                                                  self.waveform_ch1_osc1,
                                                  self.waveform_ch1_osc2,
                                                  self.waveform_ch1_osc3,
                                                  self.waveform_ch1_osc4])

        # check that oscillator waveform parameter lists all have correct dimensionality
        for i in range(5):
            assert len(self.waveform_ch0_osc_list[i]) == 3, "Error: CH0 Osc.{:d} - invalid waveform".format(i)
            assert len(self.waveform_ch1_osc_list[i]) == 3, "Error: CH1 Osc.{:d} - invalid waveform".format(i)
        # check that oscillator waveform parameters are all valid
        for i in range(5):
            assert abs(self.waveform_ch0_osc_list[i, 0]) < 1e4, "Error: CH0 Osc.{:d} - invalid frequency".format(i)
            assert abs(self.waveform_ch1_osc_list[i, 0]) < 1e4, "Error: CH1 Osc.{:d} - invalid frequency".format(i)

            assert 0. <= self.waveform_ch0_osc_list[i, 1] <= 100., "Error: CH0 Osc.{:d} - invalid amplitude".format(i)
            assert 0. <= self.waveform_ch1_osc_list[i, 1] <= 100., "Error: CH1 Osc.{:d} - invalid amplitude".format(i)

            assert -1. <= self.waveform_ch0_osc_list[i, 2] <= 1., "Error: CH0 Osc.{:d} - invalid phase".format(i)
            assert -1. <= self.waveform_ch1_osc_list[i, 2] <= 1., "Error: CH1 Osc.{:d} - invalid phase".format(i)

        # extract oscillator waveform values
        self.freq_ch0_osc_hz_list =     self.waveform_ch0_osc_list[:, 0] * kHz
        self.ampl_ch0_osc_frac_list =   self.waveform_ch0_osc_list[:, 1] / 100.
        self.phas_ch0_osc_turns_list =  self.waveform_ch0_osc_list[:, 2]

        self.freq_ch1_osc_hz_list =     self.waveform_ch1_osc_list[:, 0] * kHz
        self.ampl_ch1_osc_frac_list =   self.waveform_ch1_osc_list[:, 1] / 100.
        self.phas_ch1_osc_turns_list =  self.waveform_ch1_osc_list[:, 2]

        # check total amplitude is valid
        assert 0. <= np.sum(self.ampl_ch0_osc_frac_list) < 1., "Error: CH0 - amplitude sum is > 100%"
        assert 0. <= np.sum(self.ampl_ch1_osc_frac_list) < 1., "Error: CH1 - amplitude sum is > 100%"


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def _run_prepare(self):
        """
        Prepare phaser for output.
        """
        # wait for FIFO to clear before reset
        self.core.wait_until_mu(now_mu())
        self.core.break_realtime()
        self.core.reset()

        # set phaser attenuators
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[0].set_att(self.att_ch0_db * dB)
        delay_mu(self.t_sample_mu)
        self.phaser0.channel[1].set_att(self.att_ch1_db * dB)

        # store current status of DAC34H84's 0x1F (synchronization) register
        at_mu(self.phaser0.get_next_frame_mu())
        self.dac_register_1f = self.phaser0.dac_read(0x1F)
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def run(self):
        # PREPARE (HARDWARE)
        self._run_prepare()

        # MAIN LOOP
        for i in self._iter_repetitions:
            self.core.break_realtime()

            # configure oscillator frequencies
            self.phaser_configure()
            # clear old hardware config and resync
            self.phaser_reset()
            # set oscillator waveforms
            self.phaser_run()

        # CLEAN UP
        self.core.break_realtime()
        if self.cleanup:
            # clear outputs
            self.phaser_reset()

            # reset attenuators
            at_mu(self.phaser0.get_next_frame_mu())
            self.phaser0.channel[0].set_att(31.5 * dB)
            delay_mu(self.t_sample_mu)
            self.phaser0.channel[1].set_att(31.5 * dB)


    # HELPER FUNCTIONS
    @kernel(flags={"fast-math"})
    def phaser_configure(self):
        """
        Configure the tones on phaser.
        Puts the same RSB and BSB on both channels, and sets a third oscillator to 0 Hz in case dynamical decoupling is used.
        """
        # calculate phase delays between CH0 and CH1
        at_mu(now_mu() + 1000000)
        self.phase_ch1_turns =          (self.phase_inherent_ch1_turns +
                                         (self.freq_carrier_hz * self.time_latency_ch1_system_ns * ns))

        # calculate phase delays for each oscillator to account for inherent update latencies and system latencies
        # oscillator 0
        self.phase_ch0_osc0 =           self.phas_ch0_osc_turns_list[0]
        self.phase_ch1_osc0 =           self.phas_ch1_osc_turns_list[0]

        # oscillator 1
        self.phase_ch0_osc1 =           self.freq_ch0_osc_hz_list[1] * self.t_sample_mu * ns + self.phas_ch0_osc_turns_list[1]
        self.phase_ch1_osc1 =           self.freq_ch1_osc_hz_list[1] * self.t_sample_mu * ns + self.phas_ch1_osc_turns_list[1]

        # oscillator 2
        self.phase_ch0_osc2 =           self.freq_ch0_osc_hz_list[2] * 2 * self.t_sample_mu * ns + self.phas_ch0_osc_turns_list[2]
        self.phase_ch1_osc2 =           self.freq_ch1_osc_hz_list[2] * 2 * self.t_sample_mu * ns + self.phas_ch1_osc_turns_list[2]

        # oscillator 3
        self.phase_ch0_osc2 =           self.freq_ch0_osc_hz_list[3] * 3 * self.t_sample_mu * ns + self.phas_ch0_osc_turns_list[3]
        self.phase_ch1_osc2 =           self.freq_ch1_osc_hz_list[3] * 3 * self.t_sample_mu * ns + self.phas_ch1_osc_turns_list[3]

        # oscillator 4
        self.phase_ch0_osc2 =           self.freq_ch0_osc_hz_list[4] * 4 * self.t_sample_mu * ns + self.phas_ch0_osc_turns_list[4]
        self.phase_ch1_osc2 =           self.freq_ch1_osc_hz_list[4] * 4 * self.t_sample_mu * ns + self.phas_ch1_osc_turns_list[4]
        at_mu(now_mu() + 1000000)


        # set carrier offset frequency via the DUC
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[0].set_duc_frequency(self.freq_carrier_hz - self.freq_center_hz)
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[1].set_duc_frequency(self.freq_carrier_hz - self.freq_center_hz)
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.duc_stb()

        # set DUC phase delay to compensate for inter-channel latency
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[1].set_duc_phase(self.phase_ch1_turns)
        self.phaser0.duc_stb()


        # set oscillator frequencies
        time_oscillators_start_mu = self.phaser0.get_next_frame_mu()
        # set oscillator 0
        at_mu(time_oscillators_start_mu)
        with parallel:
            self.phaser0.channel[0].oscillator[0].set_frequency(self.phas_ch0_osc_turns_list[0])
            self.phaser0.channel[1].oscillator[0].set_frequency(self.phas_ch1_osc_turns_list[0])

        # set oscillator 1
        at_mu(time_oscillators_start_mu + self.t_sample_mu)
        with parallel:
            self.phaser0.channel[0].oscillator[1].set_frequency(self.phas_ch0_osc_turns_list[1])
            self.phaser0.channel[1].oscillator[1].set_frequency(self.phas_ch1_osc_turns_list[1])

        # set oscillator 2
        at_mu(time_oscillators_start_mu + 2 * self.t_sample_mu)
        with parallel:
            self.phaser0.channel[0].oscillator[2].set_frequency(self.phas_ch0_osc_turns_list[2])
            self.phaser0.channel[1].oscillator[2].set_frequency(self.phas_ch1_osc_turns_list[2])

        # set oscillator 3
        at_mu(time_oscillators_start_mu + 3 * self.t_sample_mu)
        with parallel:
            self.phaser0.channel[0].oscillator[3].set_frequency(self.phas_ch0_osc_turns_list[3])
            self.phaser0.channel[1].oscillator[3].set_frequency(self.phas_ch1_osc_turns_list[3])

        # set oscillator 4
        at_mu(time_oscillators_start_mu + 4 * self.t_sample_mu)
        with parallel:
            self.phaser0.channel[0].oscillator[4].set_frequency(self.phas_ch0_osc_turns_list[4])
            self.phaser0.channel[1].oscillator[4].set_frequency(self.phas_ch1_osc_turns_list[4])

        # add slack
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def phaser_run(self):
        # get synchronized start time
        time_start_mu = self.phaser0.get_next_frame_mu()

        # set oscillator 0
        at_mu(time_start_mu)
        with parallel:
            self.phaser0.channel[0].oscillator[0].set_amplitude_phase(amplitude=self.ampl_ch0_osc_frac_list[0], phase=self.phase_ch0_osc0, clr=0)
            self.phaser0.channel[1].oscillator[0].set_amplitude_phase(amplitude=self.ampl_ch1_osc_frac_list[0], phase=self.phase_ch1_osc0, clr=0)

            # set debug TTL
            with sequential:
                delay_mu(self.time_output_delay_mu)
                self.ttl8.on()

        # set oscillator 1
        at_mu(time_start_mu + self.t_sample_mu)
        with parallel:
            self.phaser0.channel[0].oscillator[1].set_amplitude_phase(amplitude=self.ampl_ch0_osc_frac_list[1], phase=self.phase_ch0_osc1, clr=0)
            self.phaser0.channel[1].oscillator[1].set_amplitude_phase(amplitude=self.ampl_ch1_osc_frac_list[1], phase=self.phase_ch1_osc1, clr=0)

            # set debug TTL
            with sequential:
                delay_mu(self.time_output_delay_mu)
                self.ttl9.on()

        # set oscillator 2
        at_mu(time_start_mu + 2 * self.t_sample_mu)
        with parallel:
            self.phaser0.channel[0].oscillator[2].set_amplitude_phase(amplitude=self.ampl_ch0_osc_frac_list[2], phase=self.phase_ch0_osc2, clr=0)
            self.phaser0.channel[1].oscillator[2].set_amplitude_phase(amplitude=self.ampl_ch1_osc_frac_list[2], phase=self.phase_ch1_osc2, clr=0)

        # set oscillator 3
        at_mu(time_start_mu + 3 * self.t_sample_mu)
        with parallel:
            self.phaser0.channel[0].oscillator[3].set_amplitude_phase(amplitude=self.ampl_ch0_osc_frac_list[3], phase=self.phase_ch0_osc3, clr=0)
            self.phaser0.channel[1].oscillator[3].set_amplitude_phase(amplitude=self.ampl_ch1_osc_frac_list[3], phase=self.phase_ch1_osc3, clr=0)

        # set oscillator 4
        at_mu(time_start_mu + 4 * self.t_sample_mu)
        with parallel:
            self.phaser0.channel[0].oscillator[4].set_amplitude_phase(amplitude=self.ampl_ch0_osc_frac_list[4], phase=self.phase_ch0_osc4, clr=0)
            self.phaser0.channel[1].oscillator[4].set_amplitude_phase(amplitude=self.ampl_ch1_osc_frac_list[4], phase=self.phase_ch1_osc4, clr=0)

        # wait for heating time
        delay_mu(self.time_pulse_mu)

    @kernel(flags={"fast-math"})
    def phaser_stop(self):
        """
        Stop the oscillators on phaser.
        """
        # set oscillator 0
        at_mu(time_start_mu)
        with parallel:
            self.phaser0.channel[0].oscillator[0].set_amplitude_phase(amplitude=0., phase=self.phase_ch0_osc0, clr=1)
            self.phaser0.channel[1].oscillator[0].set_amplitude_phase(amplitude=0., phase=self.phase_ch1_osc0, clr=1)

        # set oscillator 1
        at_mu(time_start_mu + self.t_sample_mu)
        with parallel:
            self.phaser0.channel[0].oscillator[1].set_amplitude_phase(amplitude=0., phase=self.phase_ch0_osc1, clr=1)
            self.phaser0.channel[1].oscillator[1].set_amplitude_phase(amplitude=0., phase=self.phase_ch1_osc1, clr=1)

        # set oscillator 2
        at_mu(time_start_mu + 2 * self.t_sample_mu)
        with parallel:
            self.phaser0.channel[0].oscillator[2].set_amplitude_phase(amplitude=0., phase=self.phase_ch0_osc2, clr=1)
            self.phaser0.channel[1].oscillator[2].set_amplitude_phase(amplitude=0., phase=self.phase_ch1_osc2, clr=1)

        # set oscillator 3
        at_mu(time_start_mu + 3 * self.t_sample_mu)
        with parallel:
            self.phaser0.channel[0].oscillator[3].set_amplitude_phase(amplitude=0., phase=self.phase_ch0_osc3, clr=1)
            self.phaser0.channel[1].oscillator[3].set_amplitude_phase(amplitude=0., phase=self.phase_ch1_osc3, clr=1)

        # set oscillator 4
        at_mu(time_start_mu + 4 * self.t_sample_mu)
        with parallel:
            self.phaser0.channel[0].oscillator[4].set_amplitude_phase(amplitude=0., phase=self.phase_ch0_osc4, clr=1)
            self.phaser0.channel[1].oscillator[4].set_amplitude_phase(amplitude=0., phase=self.phase_ch1_osc4, clr=1)

        # send debug triggers
        with parallel:
            self.ttl8.off()
            self.ttl9.off()

    @kernel(flags={"fast-math"})
    def phaser_reset(self):
        """
        Clear and reset the signal sources of the phaser.
        """
        # clear oscillators
        at_mu(self.phaser0.get_next_frame_mu())
        with parallel:
            self.phaser0.channel[0].oscillator[0].set_amplitude_phase(amplitude=0., phase=0., clr=1)
            self.phaser0.channel[1].oscillator[0].set_amplitude_phase(amplitude=0., phase=0., clr=1)
            delay_mu(self.t_sample_mu)
        with parallel:
            self.phaser0.channel[0].oscillator[1].set_amplitude_phase(amplitude=0., phase=0., clr=1)
            self.phaser0.channel[1].oscillator[1].set_amplitude_phase(amplitude=0., phase=0., clr=1)
            delay_mu(self.t_sample_mu)
        with parallel:
            self.phaser0.channel[0].oscillator[2].set_amplitude_phase(amplitude=0., phase=0., clr=1)
            self.phaser0.channel[1].oscillator[2].set_amplitude_phase(amplitude=0., phase=0., clr=1)
            delay_mu(self.t_sample_mu)
        with parallel:
            self.phaser0.channel[0].oscillator[3].set_amplitude_phase(amplitude=0., phase=0., clr=1)
            self.phaser0.channel[1].oscillator[3].set_amplitude_phase(amplitude=0., phase=0., clr=1)
            delay_mu(self.t_sample_mu)
        with parallel:
            self.phaser0.channel[0].oscillator[4].set_amplitude_phase(amplitude=0., phase=0., clr=1)
            self.phaser0.channel[1].oscillator[4].set_amplitude_phase(amplitude=0., phase=0., clr=1)
            delay_mu(self.t_sample_mu)

        # clear DAC phase accumulator
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.dac_write(0x1F, self.dac_register_1f & ~np.int32(1 << 1))
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.dac_write(0x1F, self.dac_register_1f | (1 << 1))

        # align DUCs of both channels by clearing their phase accumulators
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[0].set_duc_cfg(clr_once=1)
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[1].set_duc_cfg(clr_once=1)
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.duc_stb()

        self.core.break_realtime()

    def analyze(self):
        print("\tconfig:")
        print("\tch1 global latency: {:.3f}\n".format(self.phase_ch1_turns))

        print("\tosc0:")
        print("\t\tphase ch0 osc0: {:.3f}\n".format(self.phase_ch0_osc0))
        print("\t\tphase ch1 osc0: {:.3f}\n".format(self.phase_ch1_osc0))

        print("\tosc1:")
        print("\t\tphase ch0 osc1: {:.3f}".format(self.phase_ch0_osc1))
        print("\t\tphase ch1 osc1: {:.3f}\n".format(self.phase_ch1_osc1))

        print("\tosc2:")
        print("\t\tphase ch0 osc2: {:.3f}".format(self.phase_ch0_osc2))
        print("\t\tphase ch1 osc2: {:.3f}\n".format(self.phase_ch1_osc2))

        print("\tosc3:")
        print("\t\tphase ch0 osc3: {:.3f}".format(self.phase_ch0_osc3))
        print("\t\tphase ch1 osc3: {:.3f}\n".format(self.phase_ch1_osc3))

        print("\tosc4:")
        print("\t\tphase ch0 osc4: {:.3f}".format(self.phase_ch0_osc4))
        print("\t\tphase ch1 osc4: {:.3f}\n".format(self.phase_ch1_osc4))
