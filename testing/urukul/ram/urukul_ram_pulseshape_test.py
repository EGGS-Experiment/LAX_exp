from artiq.experiment import *
from artiq.coredevice.ad9910 import *
from artiq.coredevice.ad9910 import _AD9910_REG_CFR1
from artiq.coredevice.urukul import DEFAULT_PROFILE

import numpy as np


class UrukulRAMAmplitude(EnvExperiment):
    """
    Urukul RAM Amplitude Test
    Test amplitude modulation via RAM on the AD9910.
    """


    def build(self):
        self.setattr_device("core")
        self.setattr_device("core_dma")
        self.setattr_device("scheduler")

        # experiment arguments
        self.setattr_argument("repetitions",            NumberValue(default=100000, precision=0, step=1, min=1, max=10000))

        # DDS parameters
        self.setattr_argument("dds_name",               StringValue(default='urukul0_ch3'), group='dds')
        self.setattr_argument("att_dds_db",             NumberValue(default=3., precision=1, step=0.5, min=0., max=31.5), group='dds')
        self.setattr_argument("freq_dds_mhz",           NumberValue(default=100., precision=6, step=0.5, min=0., max=400), group='dds')
        self.setattr_argument("ampl_dds_max_pct",       NumberValue(default=50., precision=3, step=5., min=0., max=100.), group='dds')
        self.setattr_argument("phas_dds_rev_turns",     NumberValue(default=0.0, precision=3, step=0.1, min=-1., max=1.), group='dds')

        # modulation parameters
        self.setattr_argument("sample_rate_khz",        NumberValue(default=10000, precision=1, step=1000, min=1., max=150000), group='modulation')
        self.setattr_argument("time_pulse_us",          NumberValue(default=2, precision=1, step=1000, min=1., max=150000), group='modulation')
        self.setattr_argument("time_body_us",           NumberValue(default=2, precision=1, step=1000, min=1., max=150000), group='modulation')

        # debug triggers
        self.setattr_device("ttl8")
        self.setattr_device("ttl9")


    def prepare(self):
        """
        Prepare values for speedy evaluation.
        """
        '''GET DEVICES'''
        try:
            self.dds = self.get_device(self.dds_name)
            # add DDS name to kernel invariants
            kernel_invariants = getattr(self, "kernel_invariants", set())
            self.kernel_invariants = kernel_invariants | {'dds'}
        except Exception as e:
            print("Error: invalid DDS channel.")
            raise e

        '''CONVERT VALUES TO MACHINE UNITS'''
        # DDS values
        self.att_dds_mu =           self.dds.cpld.att_to_mu(self.att_dds_db * dB)
        self.freq_dds_ftw =         self.dds.frequency_to_ftw(self.freq_dds_mhz * MHz)
        self.phas_dds_rev_pow =     self.dds.turns_to_pow(self.phas_dds_rev_turns)

        # timing
        self.time_pulse_mu =    self.core.seconds_to_mu(self.time_pulse_us * us)
        self.time_body_mu =     self.core.seconds_to_mu(self.time_body_us * us)
        self.time_holdoff_mu =  self.core.seconds_to_mu(2000. * us)


        '''SPECFIY RAM PARAMETERS'''
        # specify RAM profile values
        # note: MUST USE PROFILE0 FOR BIDIRECTIONAL RAMP
        _MAX_RAM_LENGTH =                   1000    # actually 1024 words long, but add some
        self.ram_profile =                  0
        self.ram_addr_start =               0x00

        # create modulation array
        self.data_modulation_length =           round((self.sample_rate_khz * kHz) * (self.time_pulse_us * us))

        # ensure that waveform array stays within the RAM
        if self.data_modulation_length >= _MAX_RAM_LENGTH:
            raise Exception("Error: Waveform is too long ({:d} words).".format(self.data_modulation_length))
        elif (self.ram_addr_start + self.data_modulation_length) >= _MAX_RAM_LENGTH:
            raise Exception("Error: Max RAM address exceeds allowable range (address: {:d}).".format(self.data_modulation_length + self.ram_addr_start))

        # convert specified waveform sample rate to multiples of the SYNC_CLK (i.e. waveform update clock) period
        self.freq_dds_sync_clk_hz =             1e9 / 4.    # SYNC_CLK HAS 4ns PERIOD
        self.data_modulation_step_num_clks =    round(self.freq_dds_sync_clk_hz / (self.sample_rate_khz * kHz))

        # calculate time for ram sequence
        self.time_ramp_mu =                     self.core.seconds_to_mu((1. / self.freq_dds_sync_clk_hz) *
                                                                        self.data_modulation_step_num_clks *
                                                                        self.data_modulation_length)


        '''CALCULATE WAVEFORM'''
        # scale waveform
        _wav_x_scale =                      np.pi / 2.
        _wav_y_scale =                      self.ampl_dds_max_pct / 100.

        # create scaled x-axis
        _wav_x_vals =                       np.linspace(0., 1., self.data_modulation_length) * _wav_x_scale
        # calculate waveform y-values, then normalize and rescale
        _wav_y_vals =                       np.array([self._waveform_calc(x_val) for x_val in _wav_x_vals])
        _wav_y_vals =                       _wav_y_vals / np.max(_wav_y_vals) * _wav_y_scale

        # create empty array to store values
        self.data_modulation_arr =          [np.int32(0)] * self.data_modulation_length
        # convert amplitude data to RAM in ampl. mod. mode (i.e. 64-bit word) and store in data_modulation_arr
        self.dds.amplitude_to_ram(_wav_y_vals, self.data_modulation_arr)
        # pre-reverse data_modulation_arr since write_ram makes a booboo and reverses the array
        self.data_modulation_arr =          self.data_modulation_arr[::-1]

    def _waveform_calc(self, x: TFloat) -> TFloat:
        """
        User function that returns the waveform shape.
        """
        return np.sin(x) ** 2.
        # return 1.


    """
    MAIN FUNCTIONS
    """
    @kernel(flags={"fast-math"})
    def run_initialize(self) -> TNone:
        # reset core
        self.core.reset()
        self.core.break_realtime()

        # initialize urukul (idk why but everyone does it each time)
        self.dds.cpld.init()
        self.core.break_realtime()
        delay_mu(1000000)
        # also initialize DDS for fun? i guess?
        self.dds.init()
        self.core.break_realtime()
        delay_mu(1000000)

        # disable RAM mode
        self.dds.set_cfr1()
        self.dds.cpld.io_update.pulse_mu(8)

        # set up essential DDS parameters
        self.dds.sw.off()
        self.dds.set_att_mu(self.att_dds_mu)
        self.dds.set_cfr2(matched_latency_enable=1)

        # set up non-modulated DDS waveform parameters
        self.dds.set_ftw(self.freq_dds_ftw)
        self.dds.cpld.io_update.pulse_mu(8)
        self.dds.set_pow(0x00)
        self.dds.cpld.io_update.pulse_mu(8)

        # set up debug TTLs
        self.ttl8.off()
        self.ttl9.off()
        self.core.break_realtime()


    @kernel(flags={"fast-math"})
    def run(self) -> TNone:
        # prepare devices for experiment
        self.run_initialize()


        """
        PREPARE RAM
        """
        # configure parameters for RAM profile
        self.dds.set_profile_ram(
            start=self.ram_addr_start, end=self.ram_addr_start + (self.data_modulation_length - 1),
            step=self.data_modulation_step_num_clks,
            profile=self.ram_profile, mode=RAM_MODE_BIDIR_RAMP
        )
        self.core.break_realtime()

        # set RAM profile
        self.dds.cpld.set_profile(self.ram_profile)
        self.dds.cpld.io_update.pulse_mu(8)
        # write waveform to RAM profile
        delay_mu(1000000)   # 1ms
        self.dds.write_ram(self.data_modulation_arr)
        delay_mu(1000000)   # 1ms
        self.core.break_realtime()


        """
        MAIN LOOP
        """
        for i in range(self.repetitions):

            '''PREPARE LOOP'''
            # add slack
            self.core.break_realtime()
            delay_mu(self.time_holdoff_mu)

            # tmp remove
            self.ttl9.on()
            # tmp remove

            # initialize as profile 0 (necessary for bidirectional ramp mode)
            self.dds.cpld.set_profile(0)
            self.dds.cpld.io_update.pulse_mu(8)

            # enable RAM mode and clear DDS phase accumulator
            self.dds.write32(_AD9910_REG_CFR1,
                               (1 << 31) |              # ram_enable
                               (RAM_DEST_ASF << 29) |   # ram_destination
                               (1 << 16) |              # select_sine_output
                               (1 << 13)                # phase_autoclear
                               )

            # prime RAM sequence by pulsing IO_UPDATE
            self.dds.cpld.io_update.pulse_mu(8)


            '''RAMP-UP'''
            time_start_mu = now_mu() & ~7

            # start ramp-up (coarse align to SYNC_CLK)
            at_mu(time_start_mu)
            self.dds.cpld.set_profile(1)

            # disable phase autoclear
            self.dds.write32(_AD9910_REG_CFR1,
                               (1 << 31) |              # ram_enable
                               (RAM_DEST_ASF << 29) |   # ram_destination
                               (1 << 16)                # select_sine_output
                               )
            self.dds.cpld.io_update.pulse_mu(8)

            # open DDS switch at appropriate time
            # todo: actually do this
            # at_mu(time_start_mu + 416 + 63 - 140)
            at_mu(time_start_mu + 416 + 63 - 140 - 244)
            self.ttl8.on()
            self.dds.sw.on()

            # wait for ramp-up to finish
            delay_mu(self.time_ramp_mu)
            self.ttl8.off()


            '''PULSE DELAY'''
            # todo: set pow reg
            # self.dds.set_pow(self.phas_dds_rev_pow)
            # wait for main pulse
            delay_mu(self.time_body_mu)


            '''RAMP-DOWN'''
            time_stop_mu = now_mu() & ~7

            # start ramp-down (coarse align to SYNC_CLK)
            at_mu(time_stop_mu)
            self.dds.cpld.set_profile(0)

            # send debug signal
            at_mu(time_stop_mu + 101)
            # self.ttl9.on()

            # wait for ramp-down to finish
            delay_mu(self.time_ramp_mu)

            # close DDS switch
            self.dds.sw.off()
            # self.ttl9.off()


            '''LOOP CLEANUP'''
            # disable ram
            self.dds.set_cfr1(ram_enable=0)
            self.dds.cpld.io_update.pulse_mu(8)

            # tmp remove
            self.ttl9.off()
            # tmp remove

            # check termination periodically
            if (i % 100) == 0:
                if self.scheduler.check_termination():
                    self.run_cleanup()
                    return

        # clean up
        self.run_cleanup()


    @kernel(flags={"fast-math"})
    def run_cleanup(self) -> TNone:
        # add slack
        self.core.break_realtime()

        # stop output
        self.dds.sw.off()
        self.dds.set_att_mu(0x00)
        self.dds.set_asf(0x00)
        self.dds.cpld.io_update.pulse_mu(8)

        # stop RAM mode
        self.dds.set_cfr1(ram_enable=0)
        self.dds.cpld.io_update.pulse_mu(8)
        self.core.break_realtime()

        # set clean output
        self.dds.cpld.set_profile(DEFAULT_PROFILE)
        self.dds.cpld.io_update.pulse_mu(8)
        self.dds.set(300. * MHz, amplitude=0., profile=DEFAULT_PROFILE)
        self.core.break_realtime()

        # synchronize timeline
        self.core.wait_until_mu(now_mu())

