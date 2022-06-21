import numpy as np
from artiq.experiment import *
from artiq.coredevice.ad9910 import RAM_DEST_FTW, RAM_MODE_DIRECTSWITCH


class urukul_frequency_switch(EnvExperiment):
    """
    Switch Urukul frequency quickly by switching profile.
    """

    def build(self):
        self.setattr_device('core')
        self.setattr_argument("device", StringValue(default="urukul0_ch0"))

    def prepare(self):
        # devices
        self.setattr_device(self.device)
        self.dds = self.get_device(self.device)
        # waveform
        self.waveform_asf = self.dds.amplitude_to_asf(0.7)
        # frequencies
        self.waveform_freq0 = self.dds.frequency_to_ftw(100 * MHz)
        self.waveform_freq1 = self.dds.frequency_to_ftw(200 * MHz)
        # timing
        self.time_delay_mu = self.core.seconds_to_mu(500 * ms)
        self.time_profile_mu = self.core.seconds_to_mu(1 * ms)
        self.time_ram_mu = self.core.seconds_to_mu(10 * ms)

        self.asf = self.urukul0_ch0.amplitude_to_asf(1)
        self.ftw = self.urukul0_ch0.frequency_to_ftw(50 * MHz)

    @kernel
    def run(self):
        self.core.reset()

        # initialize devices
        #self.dds.cpld.init()
        self.core.break_realtime()
        #self.dds.init()
        self.core.break_realtime()

        self.urukul0_ch0.set_mu(self.ftw, asf=self.asf)
        self.core.break_realtime()
        self.core.break_realtime()

        # set up waveforms on each profile
        self._setup_ram()

        # set up rf switch and attenuator
        self.dds.set_att_mu(0xff)
        self.dds.cfg_sw(1)
        self.core.break_realtime()

        # switch profiles periodically
        while True:
            delay_mu(self.time_delay_mu)
            self.dds.cpld.set_profile(0)
            delay_mu(self.time_delay_mu)
            self.dds.cpld.set_profile(1)


    @kernel
    def _setup_ram(self):
        # set ram profile 0
        self.dds.set_profile_ram(start=0, end=1, profile=0, mode=RAM_MODE_DIRECTSWITCH)
        self.dds.cpld.set_profile(0)
        self.dds.cpld.io_update.pulse_mu(8)
        delay_mu(self.time_profile_mu)
        self.dds.write_ram([self.waveform_freq0])
        delay_mu(self.time_ram_mu)
        self.dds.set_asf(self.waveform_asf)
        delay_mu(self.time_ram_mu)

        # set ram profile 1
        self.dds.set_profile_ram(start=2, end=3, profile=1, mode=RAM_MODE_DIRECTSWITCH)
        self.dds.cpld.set_profile(1)
        self.dds.cpld.io_update.pulse_mu(8)
        delay_mu(self.time_profile_mu)
        self.dds.write_ram([self.waveform_freq1])
        delay_mu(self.time_ram_mu)

        # set cfr
        self.dds.set_cfr1(ram_enable=1, ram_destination=RAM_DEST_FTW)

        # set amplitude
        self.dds.set_asf(self.waveform_asf)
        # todo: maybe have to set asf for all?
        # todo: maybe have to set_mu first?
        self.dds.cpld.io_update.pulse_mu(8)

    def analyze(self):
        pass
