import numpy as np
from artiq.experiment import *
from artiq.coredevice.ad9910 import RAM_MODE_BIDIR_RAMP, RAM_DEST_ASF, RAM_DEST_FTW

_DMA_HANDLE_URUKUL_TEST = "urukul_test_run"



class UrukulWindowingTest(EnvExperiment):
    """
    Urukul Windowing Test

    Test windowing on urukul0.
    """

    def build(self):
        """
        Ensure that the necessary devices & arguments are set.
        """
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # general
        self.setattr_argument("repetitions",                            NumberValue(default=5000, ndecimals=0, step=1, min=1, max=1000000))

        # sequence timing
        self.setattr_argument("time_pulse_ms",                          NumberValue(default=5, ndecimals=5, step=1, min=0.00001, max=10000))
        self.setattr_argument("time_sequence_delay_ms",                 NumberValue(default=1, ndecimals=5, step=1, min=1, max=10000))

        # signal parameters
        self.setattr_argument("freq_signal_mhz",                        NumberValue(default=110, ndecimals=5, step=1, min=1, max=10000))


    def prepare(self):
        # alias devices
        self.dds_board =                                                self.get_device("urukul0_cpld")
        self.dds_channel =                                              self.get_device("urukul0_ch0")

        # conversion values
        self.frequency_to_ftw =                                         (1 << 30) / 6.25e6
        self.pct_to_asf =                                               0x7FFF

        # timing values
        self.time_sample_mu =                                           int64(40)
        self.time_frame_mu =                                            self.phaser0.t_frame
        self.time_pulse_mu =                                            self.core.seconds_to_mu(self.time_pulse_ms * ms)
        self.time_sequence_delay_mu =                                   self.core.seconds_to_mu(self.time_sequence_delay_ms * ms)

        # signal values
        self.freq_sideband_mu =                                         np.int32(self.freq_sideband_mhz * MHz * self.frequency_to_ftw)

        # ensure phaser pulse time is a multiple of the phaser frame period
        # (4 ns/clock * 8 clock cycles * 10 words = 320ns)
        if self.time_pulse_mu % self.time_frame_mu:
            t_frame_multiples = round(self.time_pulse_mu / self.time_frame_mu + 0.5)
            self.time_pulse_mu = int64(self.time_frame_mu * t_frame_multiples)

        # set up windowing variables
        # double the sample period since we use two oscillators
        self.time_window_sample_mu =                                    7 * self.time_sample_mu
        num_samples =                                                   np.int32(self.time_pulse_mu / (2 * self.time_window_sample_mu))
        max_amplitude =                                                 int32(0x7FFF / 2 - 1)

        # calculate windowing values - hann window
        self.ampl_window_mu_list =                                      np.power(np.sin(np.arange(num_samples) / num_samples * np.pi), 2)
        self.ampl_window_mu_list =                                      int32(self.ampl_window_mu_list * max_amplitude)
        print('\tnum samples: {}'.format(num_samples))
        # get ram data
        self.ram_values =                                               np.linspace(0x1C28F5C2, 0x1D70A3D6, 1024, dtype=np.int32)
        pass

    @kernel
    def run(self):
        self.core.reset()

        # set ram profile
        self.dds.set_profile_ram(
            start=0,
            end=len(self.ram_values)-1,
            step=1,
            mode=RAM_MODE_BIDIR_RAMP
        )

        # set profile & update
        self.dds.cpld.set_profile(0)
        self.dds.cpld.io_update.pulse_mu(8)
        delay(10 * ms)

        # write ram
        self.dds.write_ram(self.ram_values)
        delay(100 * ms)

        # write to CFR1 to enable RAM modulation
        self.dds.set_cfr1(
            ram_enable=1,
            ram_destination=RAM_DEST_FTW
        )
        self.dds.cpld.io_update.pulse_mu(8)
        self.core.break_realtime()

        # set waveform
        self.dds.set_frequency(110 * MHz)
        self.dds.cpld.io_update.pulse_mu(8)
        self.dds.set_att(8 * dB)
        self.dds.cfg_sw(True)
        self.core.break_realtime()

        while True:
            delay(1 * ms)
            self.dds.cpld.set_profile(0)
            delay(1 * ms)
            self.dds.cpld.set_profile(1)
            self.core.break_realtime()

    @kernel(flags={'fast-math'})
    def run(self):
        # initialize
        self.core.reset()
        self.initialize_phaser()

        # record window
        self.window_record()
        handle_test = self.core_dma.get_handle(_DMA_HANDLE_PHASER_TEST)
        self.core.break_realtime()

        # continually output window
        for i in range(self.repetitions):
            self.core.break_realtime()

            # align to next phaser frame
            at_mu(self.phaser0.get_next_frame_mu())

            # output signal
            self.core_dma.playback_handle(handle_test)

            # wait a holdoff period
            delay_mu(self.time_sequence_delay_mu)

    @kernel(flags={'fast-math'})
    def initialize_phaser(self):
        # initialize phaser board
        self.core.break_realtime()
        self.phaser0.init(debug=True)
        self.core.break_realtime()

        # set nco to center eggs rf around 85 MHz exactly
        self.ph_channel.set_nco_frequency(-217.083495 * MHz)
        self.ph_channel.set_nco_phase(0.)
        self.phaser0.dac_sync()
        self.core.break_realtime()

        # trf setup, and disable rf output while we set things up
        self.ph_channel.set_att(0 * dB)
        self.ph_channel.en_trf_out(rf=0, lo=0)
        self.core.break_realtime()

        # duc
        self.ph_channel.set_duc_frequency(0 * MHz)
        self.ph_channel.set_duc_cfg()
        self.phaser0.duc_stb()
        self.core.break_realtime()

        # set oscillators (i.e. sidebands)
        self.ph_channel.oscillator[0].set_frequency_mu(-self.freq_sideband_mu)
        delay_mu(self.time_sample_mu)
        self.ph_channel.oscillator[0].set_amplitude_phase_mu(0, clr=0)
        delay_mu(self.time_sample_mu)
        self.ph_channel.oscillator[1].set_frequency_mu(self.freq_sideband_mu)
        delay_mu(self.time_sample_mu)
        self.ph_channel.oscillator[1].set_amplitude_phase_mu(0, clr=0)
        self.core.break_realtime()

        # re-enable rf output
        self.ph_channel.en_trf_out(rf=1, lo=0)
        self.core.break_realtime()

    @kernel(flags={'fast-math'})
    def window_record(self):
        # record windowing sequence
        with self.core_dma.record(_DMA_HANDLE_PHASER_TEST):
            # iterate across window amplitude values in mu
            for amp_val_mu in self.ampl_window_mu_list:
                # set amplitude for oscillator 0
                self.ph_channel.oscillator[0].set_amplitude_phase_mu(amp_val_mu)
                delay_mu(self.time_window_sample_mu)

                # set amplitude for oscillator 1
                self.ph_channel.oscillator[1].set_amplitude_phase_mu(amp_val_mu)
                delay_mu(self.time_window_sample_mu)

        self.core.break_realtime()

    def analyze(self):
        print('\t\turukul test done')
