import numpy as np
from artiq.experiment import *

class pt1(EnvExperiment):

    def build(self):
        # 82, 83, 86.7 carriers
        # 1.202, 1 kHz, 200 Hz step
        # RSB/BSB/DD 40%/40%/20%
        self.freq_carrier_hz_list =                 np.array([82]) * MHz
        self.freq_sideband_hz_list =                np.array([1200]) * kHz
        self.time_pulse_ms =                        0.1


    def _prepare(self):
        # get relevant devices
        self.setattr_device('core')
        self.setattr_device('core_dma')
        self.setattr_device('phaser0')
        self.setattr_device('ttl13')
        self.setattr_device('ttl16')
        self.setattr_device('ttl17')

        # hardcoded values
        self.t_sample_mu =                          np.int64(40)
        self.t_frame_mu =                           np.int64(320)
        self.freq_center_hz =                       85. * MHz

        # preallocate phase holders
        self.phase_ch1_system_latency_turns =       np.float(0)
        self.phase_ch0_osc1 =                       np.float(0)
        self.phase_ch0_osc2 =                       np.float(0)
        self.phase_ch1_osc0 =                       np.float(0)
        self.phase_ch1_osc1 =                       np.float(0)
        self.phase_ch1_osc2 =                       np.float(0)

        ### TIMING ###
        self.time_pulse_mu =                        self.core.seconds_to_mu(self.time_pulse_ms * ms)
        if self.time_pulse_mu % self.t_frame_mu:
            t_frame_multiples = round(self.time_pulse_mu / self.t_frame_mu + 0.5)
            self.time_pulse_mu = np.int64(self.t_frame_mu * t_frame_multiples)

        ### FREQUENCIES ###
        self.config_frequencies_list =              np.zeros((len(self.freq_carrier_hz_list) * len(self.freq_sideband_hz_list), 2), dtype=float)
        self.config_frequencies_list[:, :2] =       np.stack(np.meshgrid(self.freq_carrier_hz_list, self.freq_sideband_hz_list), -1).reshape(-1, 2)

    def prepare(self):
        self._prepare()

        # global latency values
        self.time_output_delay_mu =                 np.int64(1933)
        self.phas_ch1_global_latency_turns =        0.16
        self.time_ch1_system_latency_ns =           3.712

        # store ch0 latency values
        self.time_ch0_osc1_latency_ns =             40.
        self.time_ch0_osc2_latency_ns =             80.

        # store ch1 latency values
        self.time_ch1_osc0_latency_ns =             0.
        self.time_ch1_osc1_latency_ns =             40.
        self.time_ch1_osc2_latency_ns =             80.


    @kernel(flags={"fast-math"})
    def run(self):
        # initialize
        self.core.reset()
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[0].set_att(30. * dB)
        delay_mu(self.t_sample_mu)
        self.phaser0.channel[1].set_att(30. * dB)

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
            # 2 ms delay
            delay_mu(10)

        # tmp remove
        self.phaser_reset()

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

        # align DUCs of both channels by clearing their phase accumulators
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[0].set_duc_cfg(clr_once=1)
        delay_mu(self.t_sample_mu)
        # delay_mu(self.t_frame_mu)
        self.phaser0.channel[1].set_duc_cfg(clr_once=1)
        delay_mu(self.t_sample_mu)
        # delay_mu(self.t_frame_mu)
        self.phaser0.duc_stb()


    # HELPER FUNCTIONS
    @kernel(flags={"fast-math"})
    def phaser_configure(self, carrier_freq_hz: TFloat, sideband_freq_hz: TFloat):
        # calculate ch1 latency through system
        self.phase_ch1_system_latency_turns = (carrier_freq_hz - sideband_freq_hz) * (self.time_ch1_system_latency_ns * ns)

        # osc 0
        self.phase_ch1_osc0 = (self.phas_ch1_global_latency_turns + self.phase_ch1_system_latency_turns)

        # osc 1
        self.phase_ch0_osc1 = (sideband_freq_hz) * ((self.time_ch0_osc1_latency_ns) * ns)
        self.phase_ch1_osc1 = (sideband_freq_hz) * ((self.time_ch1_osc1_latency_ns) * ns) + (self.phas_ch1_global_latency_turns + self.phase_ch1_system_latency_turns)

        # osc 2
        self.phase_ch0_osc2 = 0.
        self.phase_ch1_osc2 = (self.phas_ch1_global_latency_turns + self.phase_ch1_system_latency_turns)

        # set carrier offset frequency via the DUC
        at_mu(self.phaser0.get_next_frame_mu())
        self.phaser0.channel[0].set_duc_frequency(carrier_freq_hz - self.freq_center_hz)
        delay_mu(self.t_sample_mu)
        # delay_mu(self.t_frame_mu)
        self.phaser0.channel[1].set_duc_frequency(carrier_freq_hz - self.freq_center_hz)
        delay_mu(self.t_sample_mu)
        # delay_mu(self.t_frame_mu)
        self.phaser0.duc_stb()

        at_mu(self.phaser0.get_next_frame_mu())
        # set osc 0: rsb
        with parallel:
            self.phaser0.channel[0].oscillator[0].set_frequency(-sideband_freq_hz)
            self.phaser0.channel[1].oscillator[0].set_frequency(-sideband_freq_hz)
            delay_mu(self.t_sample_mu)
        # set osc 1: bsb
        with parallel:
            self.phaser0.channel[0].oscillator[1].set_frequency(sideband_freq_hz)
            self.phaser0.channel[1].oscillator[1].set_frequency(sideband_freq_hz)
            delay_mu(self.t_sample_mu)
        # set osc 2: carrier
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
            self.phaser0.channel[0].oscillator[0].set_amplitude_phase(amplitude=0., phase=0., clr=0)
            self.phaser0.channel[1].oscillator[0].set_amplitude_phase(amplitude=0., phase=self.phase_ch1_osc0, clr=0)
            with sequential:
                delay_mu(self.time_output_delay_mu)
                self.ttl16.on()

        # set oscillator 1
        at_mu(time_start_mu + self.t_sample_mu)
        with parallel:
            self.phaser0.channel[0].oscillator[1].set_amplitude_phase(amplitude=0., phase=self.phase_ch0_osc1, clr=0)
            self.phaser0.channel[1].oscillator[1].set_amplitude_phase(amplitude=0., phase=self.phase_ch1_osc1, clr=0)

            with sequential:
                delay_mu(self.time_output_delay_mu)
                self.ttl13.on()

        # set oscillator 2
        at_mu(time_start_mu + 2 * self.t_sample_mu)
        with parallel:
            self.phaser0.channel[0].oscillator[2].set_amplitude_phase(amplitude=0., phase=self.phase_ch0_osc2, clr=0)
            self.phaser0.channel[1].oscillator[2].set_amplitude_phase(amplitude=0., phase=self.phase_ch1_osc2 + 0.5, clr=0)

            with sequential:
                delay_mu(self.time_output_delay_mu)
                self.ttl17.on()

        # wait for heating time
        with parallel:
            self.ttl13.off()
            self.ttl16.off()
            self.ttl17.off()
            delay_mu(self.time_pulse_mu)


    def analyze(self):
        print("\tconfig:")
        print("\t\t{}".format(self.config_frequencies_list / MHz))

        print("\tosc0:")
        print("\t\tphase ch1 osc0: {:.3f}\n".format(self.phase_ch1_osc0))

        print("\tosc1:")
        print("\t\tphase ch0 osc1: {:.3f}".format(self.phase_ch0_osc1))
        print("\t\tphase ch1 osc1: {:.3f}\n".format(self.phase_ch1_osc1))

        print("\tosc2:")
        print("\t\tphase ch0 osc2: {:.3f}".format(self.phase_ch0_osc2))
        print("\t\tphase ch1 osc2: {:.3f}\n".format(self.phase_ch1_osc2))
