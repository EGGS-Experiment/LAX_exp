import numpy as np
from artiq.experiment import *


class PhaserTesting(EnvExperiment):

    def build(self):
        self.freq_carrier_hz_list =                 np.array([213.7]) * MHz
        self.freq_sideband_hz_list =                np.array([771.1]) * kHz
        self.time_pulse_ms =                        1.0


    def prepare(self):
        # get relevant devices
        self.setattr_device('core')
        self.setattr_device('core_dma')
        self.setattr_device('phaser0')
        self.setattr_device('ttl8')
        self.setattr_device('ttl9')

        # hardware values
        self.t_sample_mu =                          np.int64(40)
        self.t_frame_mu =                           np.int64(320)
        self.freq_center_hz =                       85.0 * MHz

        # system values
        self.time_output_delay_mu =                 np.int64(1953)
        self.phase_inherent_ch1_turns =             -0.425
        self.time_latency_ch1_system_ns =           2.33

        # preallocate DAC34H84 0x1F register config to speed up
        # writes when clearing the phase
        self.dac_register_1f =                      np.int32(0)

        # preallocate phase holders
        self.phase_ch1_turns =                      np.float(0)

        self.phase_ch0_osc0 =                       np.float(0)
        self.phase_ch0_osc1 =                       np.float(0)
        self.phase_ch0_osc2 =                       np.float(0)

        self.phase_ch1_osc0 =                       np.float(0)
        self.phase_ch1_osc1 =                       np.float(0)
        self.phase_ch1_osc2 =                       np.float(0)

        ### TIMING ###
        self.time_pulse_mu =                        self.core.seconds_to_mu(self.time_pulse_ms * ms)
        if self.time_pulse_mu % self.t_frame_mu:
            t_frame_multiples =                     round(self.time_pulse_mu / self.t_frame_mu + 0.5)
            self.time_pulse_mu =                    np.int64(self.t_frame_mu * t_frame_multiples)

        ### FREQUENCIES ###
        self.config_frequencies_list =              np.zeros((len(self.freq_carrier_hz_list) *
                                                              len(self.freq_sideband_hz_list),
                                                              2), dtype=float)
        self.config_frequencies_list[:, :2] =       np.stack(np.meshgrid(self.freq_carrier_hz_list,
                                                                         self.freq_sideband_hz_list),
                                                             -1).reshape(-1, 2)

    @kernel(flags={"fast-math"})
    def run(self):
        self.core.reset()

        # initialize
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[0].set_att(8 * dB)
        delay_mu(self.t_sample_mu)
        self.phaser0.channel[1].set_att(8 * dB)

        # todo: store in a method
        # store current status of DAC34H84's 0x1F (synchronization) register
        at_mu(self.phaser0.get_next_frame_mu())
        self.dac_register_1f = self.phaser0.dac_read(0x1F)
        self.core.break_realtime()


        # run
        for config_vals in self.config_frequencies_list:
            self.core.break_realtime()

            # configure oscillator frequencies
            self.phaser_configure(config_vals[0], config_vals[1])
            self.core.break_realtime()
            # clear old hardware config and resync
            self.phaser_reset()
            # set oscillator waveforms
            self.phaser_run()

        # self.tmp()


    @kernel(flags={"fast-math"})
    def tmp(self):
        self.phaser_reset()
        # reset attenuators
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[0].set_att(31.5 * dB)
        delay_mu(self.t_sample_mu)
        self.phaser0.channel[1].set_att(31.5 * dB)


    # HELPER FUNCTIONS
    @kernel(flags={"fast-math"})
    def phaser_configure(self, carrier_freq_hz: TFloat, sideband_freq_hz: TFloat):
        """
        Configure the tones on phaser for EGGS.
        Puts the same RSB and BSB on both channels, and sets a third oscillator to 0 Hz in case dynamical decoupling is used.

        Arguments:
            carrier_freq_hz         (float)     : the maximum waiting time (in machine units) for the trigger signal.
            sideband_freq_hz        (float)     : the holdoff time (in machine units)
        """
        # calculate phase delays between CH0 and CH1
        self.phase_ch1_turns =          (self.phase_inherent_ch1_turns +
                                         (carrier_freq_hz * self.time_latency_ch1_system_ns * ns))

        # calculate phase delays for each oscillator to account for inherent update latencies and system latencies
        # oscillator 0 (RSB)
        self.phase_ch0_osc0 =           0.
        self.phase_ch1_osc0 =           0.

        # oscillator 1 (BSB)
        self.phase_ch0_osc1 =           sideband_freq_hz * self.t_sample_mu * ns
        self.phase_ch1_osc1 =           sideband_freq_hz * self.t_sample_mu * ns

        # oscillator 2 (carrier) (note: ch1 has 0.5 turns to put carrier in dipole config)
        self.phase_ch0_osc2 =           0.
        self.phase_ch1_osc2 =           0.
        self.core.break_realtime()


        # set carrier offset frequency via the DUC
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[0].set_duc_frequency(carrier_freq_hz - self.freq_center_hz)
        # delay_mu(self.t_frame_mu)
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[1].set_duc_frequency(carrier_freq_hz - self.freq_center_hz)
        # delay_mu(self.t_frame_mu)
        at_mu(self.phaser0.get_next_frame_mu())
        # strobe updates for both channels
        self.phaser0.duc_stb()

        # set DUC phase delay to compensate for inter-channel latency
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[1].set_duc_phase(self.phase_ch1_turns)
        self.phaser0.duc_stb()


        # set sideband frequencies
        at_mu(self.phaser0.get_next_frame_mu())
        # set oscillator 0 (RSB)
        with parallel:
            self.phaser0.channel[0].oscillator[0].set_frequency(-sideband_freq_hz)
            self.phaser0.channel[1].oscillator[0].set_frequency(-sideband_freq_hz)
            delay_mu(self.t_sample_mu)
        # set oscillator 1 (BSB)
        with parallel:
            self.phaser0.channel[0].oscillator[1].set_frequency(sideband_freq_hz)
            self.phaser0.channel[1].oscillator[1].set_frequency(sideband_freq_hz)
            delay_mu(self.t_sample_mu)
        # set oscillator 2 (carrier)
        with parallel:
            self.phaser0.channel[0].oscillator[2].set_frequency(0.)
            self.phaser0.channel[1].oscillator[2].set_frequency(0.)
            delay_mu(self.t_sample_mu)

    @kernel(flags={"fast-math"})
    def phaser_run(self):
        # set oscillator 0
        at_mu(self.phaser0.get_next_frame_mu())
        time_start_mu = now_mu()
        with parallel:
            self.phaser0.channel[0].oscillator[0].set_amplitude_phase(amplitude=0.0, phase=self.phase_ch0_osc0, clr=0)
            self.phaser0.channel[1].oscillator[0].set_amplitude_phase(amplitude=0.0, phase=self.phase_ch1_osc0, clr=0)

            with sequential:
                delay_mu(self.time_output_delay_mu)
                self.ttl8.on()

        # set oscillator 1
        at_mu(time_start_mu + self.t_sample_mu)
        with parallel:
            self.phaser0.channel[0].oscillator[1].set_amplitude_phase(amplitude=0.0, phase=self.phase_ch0_osc1, clr=0)
            self.phaser0.channel[1].oscillator[1].set_amplitude_phase(amplitude=0.0, phase=self.phase_ch1_osc1, clr=0)

            with sequential:
                delay_mu(self.time_output_delay_mu)
                self.ttl9.on()

        # set oscillator 2
        at_mu(time_start_mu + 2 * self.t_sample_mu)
        with parallel:
            self.phaser0.channel[0].oscillator[2].set_amplitude_phase(amplitude=0.9, phase=self.phase_ch0_osc2, clr=0)
            self.phaser0.channel[1].oscillator[2].set_amplitude_phase(amplitude=0.9, phase=self.phase_ch1_osc2, clr=0)

            with sequential:
                delay_mu(self.time_output_delay_mu)

        # wait for heating time
        with parallel:
            self.ttl8.off()
            self.ttl9.off()
            delay_mu(self.time_pulse_mu)

    @kernel(flags={"fast-math"})
    def phaser_reset(self):
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


        # clear DAC phase accumulator
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.dac_write(0x1F, self.dac_register_1f & ~np.int32(1 << 1))
        # delay_mu(self.t_frame_mu)
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.dac_write(0x1F, self.dac_register_1f | (1 << 1))


        # align DUCs of both channels by clearing their phase accumulators
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[0].set_duc_cfg(clr_once=1)
        # delay_mu(self.t_frame_mu)
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[1].set_duc_cfg(clr_once=1)
        # delay_mu(self.t_frame_mu)
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.duc_stb()

    def analyze(self):
        print("\tconfig:")
        print("\t\t{}".format(self.config_frequencies_list / MHz))

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
