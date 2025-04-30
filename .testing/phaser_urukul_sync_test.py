import numpy as np
from artiq.experiment import *
from artiq.coredevice import ad9910


class phaserurukulsynctest(EnvExperiment):
    """
    phaserurukulsynctest tmp
    """
    kernel_invariants = {
        "dds", "phaser", "ttl",
        "phaser_t_sample_mu", "phaser_t_frame_mu",

        "freq_hz", "time_delays_mu",

        "time_duc_ns", "duc_extra_phase",
        "time_osc_ns", "osc_extra_phase", "freq_osc_hz_base",

        "phas_dds_turns", "ampl_dds_pct", "att_dds_db", "profile_target",
        "freq_osc_hz", "phas_osc_turns", "ampl_osc_frac", "att_phaser_db",
        "freq_dds_ftw", "phas_dds_pow", "ampl_dds_asf", "att_dds_mu"
    }

    def build(self):
        self.setattr_device("core")
        self.setattr_device('ttl15')

        self.setattr_device('phaser1')

        self.setattr_device('urukul1_ch3')
        self.setattr_device('urukul1_cpld')

        # basic values
        self.phaser_t_sample_mu =   np.int64(40)
        self.phaser_t_frame_mu =    np.int64(320)

        # central values
        self.freq_hz =          75.00 * MHz
        self.time_delays_mu =   [0, 0]

        self.time_duc_ns =      61.086857
        self.duc_extra_phase =  0.3361

        self.freq_osc_hz_base = -9 * MHz
        # self.time_osc_ns =      -898.6479
        # self.time_osc_ns =      -14345.34
        self.time_osc_ns =      663.134
        self.osc_extra_phase =  0.

        # urukul, then phaser
        self.phas_dds_turns =   0.
        self.ampl_dds_pct =     45.
        self.att_dds_db =       8.
        self.profile_target =   6

        # phaser-specific
        self.freq_osc_hz =      np.array([0., 0., 0., 0., 0.]) * MHz
        self.phas_osc_turns =   np.array([0., 0., 0., 0., 0.])
        self.ampl_osc_frac =    np.array([50., 0., 00., 0., 0.]) / 100.
        self.att_phaser_db =    1.

    def prepare(self):
        self.ttl =      self.get_device('ttl15')
        self.dds =      self.get_device('urukul1_ch3')
        self.phaser =   self.get_device('phaser1')

        # convert units for DDS
        self.freq_dds_ftw = self.dds.frequency_to_ftw(self.freq_hz)
        self.phas_dds_pow = self.dds.turns_to_pow(self.phas_dds_turns)
        self.ampl_dds_asf = self.dds.amplitude_to_asf(self.ampl_dds_pct / 100.)
        self.att_dds_mu =   self.dds.cpld.att_to_mu(self.att_dds_db * dB)

        self.time_delays_mu = [
            self.core.seconds_to_mu(val_ns * ns)
            for val_ns in self.time_delays_mu
        ]

        self.time_start_mu_phaser = np.int64(0)
        self.time_start_mu_urukul = np.int64(0)

        self.phas_duc_turns_dj = 0.
        self.time_delay_default = np.int64(0)

    @kernel(flags={'fast-math'})
    def run(self) -> TNone:
        """
        Main run kernel.
        """
        '''
        PREPARE
        '''
        # set up hardware
        self.device_setup()
        self.core.break_realtime()
        delay_mu(10000000)

        # get starting reference times
        self.time_start_mu_phaser = self.phaser.get_next_frame_mu()
        self.time_start_mu_urukul = self.time_start_mu_phaser & ~(0x7)

        '''
        SET CARRIER FREQUENCY
        '''
        with ((parallel)):

            # send debug TTL
            self.ttl.on()

            # run urukul
            with sequential:
                self.dds.set_mu(
                    self.freq_dds_ftw, asf=self.ampl_dds_asf,
                    pow_=self.phas_dds_pow, profile=self.profile_target,
                    phase_mode=ad9910.PHASE_MODE_TRACKING, ref_time_mu=self.time_start_mu_urukul
                )
                self.dds.sw.on()

            # run phaser
            with sequential:
                # clear DUC phase accumulator for both channels
                at_mu(self.time_start_mu_phaser)
                self.phaser.channel[0].set_duc_cfg(clr_once=1)
                self.phaser.channel[1].set_duc_cfg(clr_once=1)
                self.phaser.duc_stb()

                # load phase
                at_mu(self.time_start_mu_phaser + 3 * self.phaser_t_frame_mu)
                self.phas_duc_turns_dj = (
                                                 self.core.mu_to_seconds(now_mu() - self.time_start_mu_phaser) + (self.time_duc_ns * ns)
                                         ) * (self.freq_hz - self.freq_osc_hz_base) + self.duc_extra_phase
                self.time_delay_default = now_mu() - self.time_start_mu_phaser
                self.phaser.channel[0].set_duc_phase(self.phas_duc_turns_dj)
                self.phaser.duc_stb()

                # at_mu(self.time_start_mu_phaser + 4 * self.phaser_t_frame_mu)
                # self.phaser.channel[1].set_duc_phase(self.phas_duc_turns_dj)
                # self.phaser.duc_stb()

                # # set osc freqs
                # # at_mu(self.phaser.get_next_frame_mu())
                # at_mu(self.time_start_mu_phaser + 5 * self.phaser_t_frame_mu)
                # for i in range(5):
                #     with parallel:
                #         self.phaser.channel[0].oscillator[i].set_frequency(self.freq_osc_hz[i])
                #         self.phaser.channel[1].oscillator[i].set_frequency(self.freq_osc_hz[i])
                #         delay_mu(self.phaser_t_sample_mu)

                # set osc ampl/phases
                at_mu(self.time_start_mu_phaser + 6 * self.phaser_t_frame_mu)
                for i in range(5):
                    at_mu(self.time_start_mu_phaser + 6 * self.phaser_t_frame_mu + i * self.phaser_t_sample_mu)
                    osc_turns_dj = (self.core.mu_to_seconds(now_mu() - self.time_start_mu_phaser) +
                                    (self.time_osc_ns * ns)) * self.freq_osc_hz_base + self.phas_osc_turns[i]
                    with parallel:
                        # osc_turns_dj = (self.core.mu_to_seconds(6 * self.phaser_t_frame_mu) +
                        #                 (self.time_osc_ns * ns)) * self.freq_osc_hz_base + self.phas_osc_turns[i]
                        self.phaser.channel[0].oscillator[i].set_amplitude_phase(
                            amplitude=self.ampl_osc_frac[i], phase=osc_turns_dj, clr=0
                        )
                        self.phaser.channel[1].oscillator[i].set_amplitude_phase(
                            amplitude=self.ampl_osc_frac[i], phase=osc_turns_dj, clr=0
                        )
                        # delay_mu(self.phaser_t_sample_mu)

        # stop TTL
        at_mu(self.time_start_mu_phaser + 8000)
        self.ttl.off()

    @kernel(flags={'fast-math'})
    def reset_duc_phase(self, time_run_mu: TInt64) -> TNone:
        # synchronize to frame
        at_mu(time_run_mu)

        # clear DUC phase accumulator for both channels
        self.phaser.channel[0].set_duc_cfg(clr_once=1)
        delay_mu(self.phaser_t_frame_mu)
        self.phaser.channel[1].set_duc_cfg(clr_once=1)
        delay_mu(self.phaser_t_frame_mu)

        # strobe update register for both DUCs
        self.phaser.duc_stb()

    @kernel(flags={'fast-math'})
    def device_setup(self):
        """
        Set up basic functionality for devices.
        """
        # clean slate
        self.core.reset()
        self.core.break_realtime()
        delay_mu(1000000)

        # debug TTL
        self.ttl.off()

        # DDS
        self.dds.sw.off()
        self.dds.cpld.get_att_mu()
        self.core.break_realtime()

        self.dds.set_mu(
            self.freq_dds_ftw, asf=self.ampl_dds_asf,
            pow_=0, profile=self.profile_target
        )
        self.dds.cpld.io_update.pulse_mu(8)
        self.dds.set_att_mu(self.att_dds_mu)
        self.dds.cpld.set_profile(self.profile_target)
        self.dds.cpld.io_update.pulse_mu(8)
        self.core.break_realtime()

        # phaser - clear oscs
        self.phaser.channel[0].set_att(self.att_phaser_db)
        for i in range(5):
            # synchronize to frame
            at_mu(self.phaser.get_next_frame_mu())
            with parallel:
                # self.phaser.channel[0].oscillator[i].set_frequency(self.freq_osc_hz_base)
                # self.phaser.channel[1].oscillator[i].set_frequency(self.freq_osc_hz_base)
                self.phaser.channel[0].oscillator[i].set_frequency(0.)
                self.phaser.channel[1].oscillator[i].set_frequency(0.)
            self.core.break_realtime()

        for i in range(5):
            # synchronize to frame
            at_mu(self.phaser.get_next_frame_mu())
            with parallel:
                self.phaser.channel[0].oscillator[i].set_amplitude_phase(amplitude=0., phase=0., clr=1)
                self.phaser.channel[1].oscillator[i].set_amplitude_phase(amplitude=0., phase=0., clr=1)
            self.core.break_realtime()

        # phaser - set carrier values (on DUC)
        at_mu(self.phaser.get_next_frame_mu())
        self.phaser.channel[0].set_duc_frequency(self.freq_hz - self.freq_osc_hz_base)
        delay_mu(self.phaser_t_frame_mu)
        self.phaser.channel[1].set_duc_frequency(self.freq_hz - self.freq_osc_hz_base)
        delay_mu(self.phaser_t_frame_mu)
        # strobe updates for both channels
        self.phaser.duc_stb()

        # set DUC phase delay to compensate for inter-channel latency
        at_mu(self.phaser.get_next_frame_mu())
        self.phaser.channel[1].set_duc_phase(0.)
        self.phaser.duc_stb()
        self.core.break_realtime()

        # at_mu(self.phaser.get_next_frame_mu())
        # for i in range(5):
        #     with parallel:
        #         self.phaser.channel[0].oscillator[i].set_amplitude_phase(
        #             amplitude=self.ampl_osc_frac[i],
        #             phase=self.phas_osc_turns[i], clr=0
        #         )
        #         self.phaser.channel[1].oscillator[i].set_amplitude_phase(
        #             amplitude=self.ampl_osc_frac[i],
        #             phase=self.phas_osc_turns[i], clr=0
        #         )
        #         delay_mu(self.phaser_t_sample_mu)

        # phaser - clear oscs
        self.phaser.channel[0].set_att(self.att_phaser_db)
        for i in range(5):
            # synchronize to frame
            # at_mu(self.phaser.get_next_frame_mu())
            with parallel:
                self.phaser.channel[0].oscillator[i].set_frequency(self.freq_osc_hz_base)
                self.phaser.channel[1].oscillator[i].set_frequency(self.freq_osc_hz_base)
                delay_mu(self.phaser_t_sample_mu)
                # self.phaser.channel[0].oscillator[i].set_frequency(0.)
                # self.phaser.channel[1].oscillator[i].set_frequency(0.)
            self.core.break_realtime()

    def analyze(self):
        # print(self.time_start_mu_phaser)
        # print(self.time_start_mu_urukul)
        print(self.phas_duc_turns_dj)
        print(self.time_delay_default)
        pass

