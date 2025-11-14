from artiq.experiment import *
from artiq.coredevice import ad9910
from artiq.coredevice import urukul
from artiq.coredevice import spi2 as spi

import numpy as np

CS_MASK_NU = 3


class urukul_masknu(EnvExperiment):
    """
    urukul_masknu
    Testing masknu
    """
    kernel_invariants = {
        "dds", "ttl",

        "dds_list", "dds_mask_list",
        "att_db", "att_mu",

        "profile0",
        "freq_mhz0", "ampl_pct0", "phas_turns0",
        "freq_mhz1", "ampl_pct1", "phas_turns1",

        "profile1",
        "freq_ftw0", "ampl_asf0", "phas_pow0",
        "freq_ftw1", "ampl_asf1", "phas_pow1",

        "_word_mask_nu",
    }

    def build(self):
        self.setattr_device("core")
        self.setattr_device("tt15")

        for i in range(2):
            self.setattr_device("urukul{:d}_cpld".format(i))
            for j in range(4):
                self.setattr_device("urukul{:d}_ch{:d}".format(i, j))

        self.att_db =       10.

        self.profile0 =     5
        self.freq_mhz0 =    100.
        self.ampl_pct0 =    50.
        self.phas_turns0 =  0.

        self.profile1 =     6
        self.freq_mhz1 =    100.
        self.ampl_pct1 =    10.
        self.phas_turns1 =  0.

        self.ttl =  self.ttl15
        self.dds =  self.urukul0_ch0

        self.dds_list =         [self.urukul0_ch0, self.urukul0_ch1, self.urukul0_ch2, self.urukul0_ch3]
        self.dds_mask_list =    [self.urukul0_ch0, self.urukul0_ch1, self.urukul0_ch2]

    def prepare(self):
        # convert values to machine units
        self.att_mu = self.dds.cpld.att_to_mu(self.att_db * dB)

        self.freq_ftw0 = self.dds.frequency_to_ftw(self.freq_mhz0 * MHz)
        self.ampl_asf0 = self.dds.amplitude_to_asf(self.ampl_pct0 / 100.)
        self.phas_pow0 = self.dds.turns_to_pow(self.phas_turns0)

        self.freq_ftw1 = self.dds.frequency_to_ftw(self.freq_mhz1 * MHz)
        self.ampl_asf1 = self.dds.amplitude_to_asf(self.ampl_pct1 / 100.)
        self.phas_pow1 = self.dds.turns_to_pow(self.phas_turns1)

        # convert dds_mask_list to mask_nu word
        self._word_mask_nu = 0x0
        for dds in self.dds_mask_list:
            self._word_mask_nu |= (1 << dds.chip_select - 4)
        # update DDS with chip_select (?)

    @kernel(flags={"fast-math"})
    def _run_prepare(self) -> TNone:
        self.core.wait_until_mu(now_mu())
        self.core.reset()

        # prepare everything else
        self.ttl.off()
        self.dds.get_att_mu()
        self.core.break_realtime()

        # prepare DDS waveforms
        time_ref_mu = now_mu() & ~0x7
        for dds in self.dds_list:
            self.core.break_realtime()
            dds.sw.off()
            dds.set_mu(self.freq_ftw0, asf=self.ampl_asf0, _pow=self.phas_pow0,
                       profile=self.profile0,
                       phase_mode=ad9910.PHASE_MODE_TRACKING,
                       ref_time_mu=time_ref_mu)
            delay_mu(5000)
            dds.set_mu(self.freq_ftw0, asf=self.ampl_asf0, _pow=self.phas_pow0,
                       profile=self.profile1,
                       phase_mode=ad9910.PHASE_MODE_TRACKING,
                       ref_time_mu=time_ref_mu)

            self.core.break_realtime()
            dds.set_att_mu(self.att_mu)
            dds.set_cfr2(matched_latency_enable=1)
            delay_mu(1000)
            dds.sw.on()

            self.core.wait_until_mu(now_mu())

        # set target dds profile
        self.core.break_realtime()
        self.dds.cpld.set_profile(self.profile0)
        self.dds.cpld.io_update.pulse_mu(8)

    @kernel(flags={"fast-math"})
    def run(self) -> TNone:
        self._run_prepare()

        '''
        MAIN SEQUENCE
        '''
        ### mask_nu protocol ###
        # set mask_nu for relevant devices
        self.dds.cpld.cfg_write(self.cfg_reg | (self._word_mask_nu << urukul.CFG_MASK_NU))

        # write to masked DDSs
        delay_mu(16)
        at_mu(now_mu() & ~0x7)
        self.write64(ad9910._AD9910_REG_PROFILE0 + self.profile1,
                     (self.ampl_asf1 << 16) | (self.phas_pow0 & 0xffff), self.freq_ftw0)
        self.write32(ad9910._AD9910_REG_CFR1, 2 | (1 << 13)) # phase_autoclear


        ### latch change and observe: mask_nu (via profile update) ###
        self.dds.cpld.set_profile(self.profile1)
        self.ttl.on()   # set TTL to trigger scope
        delay_mu(10000) # 10us
        self.ttl.off()  # clear TTL for next update

        ### latch change and observe: io_update effect w/mask_nu ###
        delay_mu(25000)
        self.dds.cpld.io_update.pulse_mu(8)
        self.ttl.on()   # set TTL to trigger scope
        delay_mu(10000) # 10us
        self.ttl.off()  # clear TTL for next update

        ### clean up ###
        self.dds.cpld.cfg_write(self.cfg_reg | (0 << urukul.CFG_MASK_NU))
        self.core.wait_until_mu(now_mu())
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def write32(self, addr: TInt32, data: TInt32) -> TNone:
        self.bus.set_config_mu(urukul.SPI_CONFIG, 8,
                               urukul.SPIT_DDS_WR, CS_MASK_NU)
        self.bus.write(addr << 24)
        self.bus.set_config_mu(urukul.SPI_CONFIG | spi.SPI_END, 32,
                               urukul.SPIT_DDS_WR, CS_MASK_NU)
        self.bus.write(data)

    @kernel(flags={"fast-math"})
    def write64(self, addr: TInt32, data_high: TInt32, data_low: TInt32) -> TNone:
        self.bus.set_config_mu(urukul.SPI_CONFIG, 8,
                               urukul.SPIT_DDS_WR, CS_MASK_NU)
        self.bus.write(addr << 24)
        self.bus.set_config_mu(urukul.SPI_CONFIG, 32,
                               urukul.SPIT_DDS_WR, CS_MASK_NU)
        self.bus.write(data_high)
        self.bus.set_config_mu(urukul.SPI_CONFIG | spi.SPI_END, 32,
                               urukul.SPIT_DDS_WR, CS_MASK_NU)
        self.bus.write(data_low)






