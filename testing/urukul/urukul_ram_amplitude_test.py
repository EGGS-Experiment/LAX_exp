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
        self.setattr_device("ttl8")
        self.setattr_device("ttl9")

        self.setattr_argument("repetitions",            NumberValue(default=20, ndecimals=0, step=1, min=1, max=10000))

        # DDS parameters
        self.setattr_argument("dds_name",               StringValue(default='urukul0_ch2'), group='dds')
        self.setattr_argument("att_dds_db",             NumberValue(default=10, ndecimals=1, step=0.5, min=0., max=31.5), group='dds')
        self.setattr_argument("freq_dds_mhz",           NumberValue(default=85, ndecimals=6, step=0.5, min=0., max=400), group='dds')
        self.setattr_argument("ampl_dds_max_pct",       NumberValue(default=50., ndecimals=3, step=5., min=0., max=100.), group='dds')

        # modulation parameters
        self.setattr_argument("sample_rate_khz",        NumberValue(default=1000, ndecimals=1, step=1000, min=1., max=150000), group='modulation')
        self.setattr_argument("time_pulse_us",          NumberValue(default=250, ndecimals=1, step=1000, min=1., max=150000), group='modulation')


    def prepare(self):
        """
        Prepare values for speedy evaluation.
        """
        '''GET DEVICES'''
        try:
            self.dds = self.get_device(self.dds_name)
        except Exception as e:
            print("Error: invalid DDS channel.")
            raise e

        '''CONVERT VALUES TO MACHINE UNITS'''
        self.att_dds_mu =       self.dds.cpld.att_to_mu(self.att_dds_db * dB)
        self.freq_dds_ftw =     self.dds.frequency_to_ftw(self.freq_dds_mhz * MHz)
        self.time_pulse_mu =    self.core.seconds_to_mu(self.time_pulse_us * us)

        '''SPECFIY RAM PARAMETERS'''
        # specify RAM profile values
        # note: MUST BE PROFILE0 FOR BIDIRECTIONAL RAMP
        self.ram_profile =                  0
        self.ram_addr_start =               0x00

        # create modulation array
        self.data_modulation_length =       100

        # calculate step timings for ram profile (sample rate to num steps)
        self.freq_dds_sync_clk_hz =         1e9 / 4.
        self.data_modulation_step_rate =    int(self.freq_dds_sync_clk_hz / (self.sample_rate_khz * kHz))

        # calculate time for ram sequence
        self.time_ramp_mu =                 self.core.seconds_to_mu(self.data_modulation_step_rate *
                                                                    self.freq_dds_sync_clk_hz *
                                                                    self.data_modulation_length)

        '''CALCULATE WAVEFORM'''
        # create x-axis
        _data_x_vals =                      np.linspace(0., 1., 512)
        # calculate sine-squared window in fractional units
        self.data_modulation_arr =          (np.sin(_data_x_vals * (np.pi / 2.)) ** 2.) * (self.ampl_dds_max_pct / 100.)

        # convert fractional amplitude values to asf
        self.data_modulation_arr =          np.array([self.dds.amplitude_to_asf(ampl_frac) for ampl_frac in self.data_modulation_arr])


    """
    MAIN FUNCTIONS
    """
    @kernel(flags={"fast-math"})
    def run_initialize(self):
        # reset core
        self.core.reset()
        self.core.break_realtime()

        # set up essential DDS parameters
        self.dds.sw.off()
        self.dds.set_att_mu(self.att_dds_mu)
        self.dds.set_cfr2(matched_latency_enable=1)

        # set up debug TTLs
        self.ttl8.off()
        self.ttl9.off()
        self.core.break_realtime()


    @kernel
    def run(self):
        # initialize
        self.run_initialize()

        # configure parameters for RAM profile
        self.dds.set_profile_ram(
            start=self.ram_addr_start, end=self.ram_addr_start + (self.data_modulation_length - 1),
            step=self.data_modulation_step_rate,
            profile=self.ram_profile, mode=RAM_MODE_BIDIR_RAMP
        )

        # set RAM profile
        self.dds.cpld.set_profile(self.ram_profile)
        self.dds.cpld.io_update.pulse_mu(8)
        self.core.break_realtime()

        # write waveform to RAM
        self.dds.write_ram(self.data_modulation_arr)
        delay_mu(10000000)
        self.core.break_realtime()

        # set desired waveform parameters
        self.dds.set_ftw(self.freq_dds_ftw)
        self.dds.cpld.io_update.pulse_mu(8)
        self.core.break_realtime()


        """
        MAIN LOOP
        """
        for i in range(self.repetitions):
            self.core.break_realtime()

            # tmp remove - leave switch open to see dynamics
            self.dds.sw.on()

            # prepare DDS output & enable RAM
            self.dds.write32(_AD9910_REG_CFR1,
                               (1 << 31) |              # ram_enable
                               (RAM_DEST_ASF << 29) |   # ram_destination
                               (1 << 16)              # select_sine_output
                               # (1 << 13)                # phase_autoclear
                               )


            '''START RAMP-UP'''
            # coarse align to previous SYNC_CLK period
            time_start_mu = now_mu() & ~7

            # set start profile
            at_mu(time_start_mu)
            self.dds.cpld.set_profile(1)
            # start ramp-up (or really, restart it, since setting profile 0 should have set it)
            at_mu(time_start_mu)
            self.dds.cpld.io_update.pulse_mu(8)

            # open DDS switch at appropriate time
            at_mu(time_start_mu + 416 + 63 - 140)
            self.dds.sw.on()
            # send debug trigger
            self.ttl8.on()

            # wait for ramp-up to finish
            delay_mu(self.time_ramp_mu)


            '''PULSE DELAY'''
            # wait for main pulse
            delay_mu(self.time_pulse_mu)


            '''START RAMP-DOWN'''
            # coarse align to previous SYNC_CLK period
            time_stop_mu = now_mu() & ~7

            # set stop profile
            at_mu(time_stop_mu)
            self.dds.cpld.set_profile(0)
            # start ramp-down (or really, restart it, since setting profile 0 should have set it)
            at_mu(time_stop_mu)
            self.dds.cpld.io_update.pulse_mu(8)

            # wait for ramp-down to finish
            delay_mu(self.time_ramp_mu)

            # close DDS switch
            self.dds.sw.off()
            # send debug trigger
            self.ttl8.off()


            # # disable ram
            # self.dds.set_cfr1(ram_enable=0)
            # self.dds.cpld.io_update.pulse_mu(8)

        # clean up
        self.run_cleanup()

    @kernel(flags={"fast-math"})
    def run_cleanup(self):
        # ensure clean output
        self.core.break_realtime()
        self.dds.sw.off()
        self.dds.set_att_mu(0x00)

        # stop RAM mode
        self.dds.set_cfr1(ram_enable=0)
        self.dds.cpld.io_update.pulse_mu(8)
        self.core.break_realtime()

        # set clean output
        self.dds.cpld.set_profile(DEFAULT_PROFILE)
        self.dds.cpld.io_update.pulse_mu(8)
        self.dds.set(250 * MHz, amplitude=0., profile=DEFAULT_PROFILE)
        self.core.break_realtime()


