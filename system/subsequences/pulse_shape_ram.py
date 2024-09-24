from artiq.experiment import *
from artiq.coredevice.ad9910 import *
from artiq.coredevice.ad9910 import _AD9910_REG_CFR1
from artiq.coredevice.urukul import DEFAULT_PROFILE

from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence


class PulseShapeRAM(LAXSubsequence):
    """
    Subsequence: Pulse Shape RAM

    Transfer the population in the S-1/2 mj=-1/2 state to the D-5/2 mj=-5/2 state,
    and vice versa using the polarized 729nm beam.
    """
    name = 'pulse_shape_ram'
    kernel_invariants = {
        "time_rabiflop_mu"
    }

    def build_subsequence(self):
        # modulation parameters
        self.setattr_argument("sample_rate_khz",        NumberValue(default=250, ndecimals=1, step=1000, min=1., max=150000), group='modulation')
        self.setattr_argument("time_rolloff_us",        NumberValue(default=1000, ndecimals=1, step=1000, min=1., max=150000), group='modulation')
        self.setattr_argument("time_body_us",           NumberValue(default=100, ndecimals=1, step=1000, min=1., max=150000), group='modulation')

        # get device
        self.setattr_device('qubit')
        self.setattr_device('ttl8')
        self.setattr_device('ttl9')

    def prepare_subsequence(self):
        # timing
        self.time_pulse_mu =    self.core.seconds_to_mu(self.time_pulse_us * us)
        self.time_body_mu =     self.core.seconds_to_mu(self.time_body_us * us)


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


    @kernel(flags={"fast-math"})
    def run(self) -> TNone:
        # population transfer pulse
        self.qubit.on()
        delay_mu(self.time_rabiflop_mu)
        self.qubit.off()

    @kernel(flags={"fast-math"})
    def thkim(self, data_waveform_mu: TList(TInt32)) -> TNone:
        # configure parameters for RAM profile
        self.dds.set_profile_ram(
            start=self.ram_addr_start, end=self.ram_addr_start + (len(data_waveform_mu) - 1),
            step=self.data_modulation_step_num_clks,
            profile=self.ram_profile, mode=RAM_MODE_BIDIR_RAMP
        )
        self.core.break_realtime()

        # set RAM profile
        self.dds.cpld.set_profile(self.ram_profile)
        self.dds.cpld.io_update.pulse_mu(8)
        # write waveform to RAM profile
        delay_mu(1000000)   # 1ms
        self.dds.write_ram(data_waveform_mu)
        delay_mu(1000000)   # 1ms
        self.core.break_realtime()


    @kernel(flags={"fast-math"})
    def fire(self) -> TNone:
        # initialize as profile 0 (necessary for bidirectional ramp mode)
        self.dds.cpld.set_profile(0)
        self.dds.cpld.io_update.pulse_mu(8)

        # enable RAM mode and clear DDS phase accumulator
        self.dds.write32(_AD9910_REG_CFR1,
                         (1 << 31) |  # ram_enable
                         (RAM_DEST_ASF << 29) |  # ram_destination
                         (1 << 16) |  # select_sine_output
                         (1 << 13)  # phase_autoclear
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
                         (1 << 31) |  # ram_enable
                         (RAM_DEST_ASF << 29) |  # ram_destination
                         (1 << 16)  # select_sine_output
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
        self.ttl9.on()

        # wait for ramp-down to finish
        delay_mu(self.time_ramp_mu)

        # close DDS switch
        self.dds.sw.off()
        self.ttl9.off()

        '''LOOP CLEANUP'''
        # disable ram
        self.dds.set_cfr1(ram_enable=0)
        self.dds.cpld.io_update.pulse_mu(8)

