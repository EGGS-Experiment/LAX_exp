from artiq.experiment import *

# todo: duty cycle


class DDSClock(EnvExperiment):
    """
    DDS Clock
    Use the DIO-BNC TTLs to create a clock pulse.
    """

    def build(self):
        """
        Set devices and arguments for the experiment.
        """
        self.setattr_device("core")
        self.setattr_device("core_dma")
        self.time_trigger_delay_mu =                                self.core.seconds_to_mu(50.88 * us)

        # timing
        self.setattr_argument('frequency_clock_khz',                NumberValue(default=1, ndecimals=3, step=1, min=0.001, max=1000))
        self.setattr_argument('duty_cycle_frac',                    NumberValue(default=0.5, ndecimals=6, step=1, min=0.0000001, max=1000))
        self.setattr_argument('time_total_s',                       NumberValue(default=10, ndecimals=6, step=1, min=0.00001, max=100000000))
        #self.setattr_argument('time_off_us', NumberValue(default=100, ndecimals=3, step=1, min=0, max=100))

        # TTL channel
        self.setattr_argument("ttl_channel",                        NumberValue(default=14, ndecimals=0, step=1, min=8, max=23))
        self.setattr_argument("ttl_trigger",                        NumberValue(default=15, ndecimals=0, step=1, min=8, max=23))


    def prepare(self):
        """
        Set up the dataset and prepare things such that the kernel functions have minimal overhead.
        """
        # devices
        self.ttl_clock =                                            self.get_device("ttl{:d}".format(self.ttl_channel))
        self.ttl_trigger =                                          self.get_device("ttl{:d}".format(self.ttl_trigger))

        # timing
        self.time_on_mu = self.core.seconds_to_mu(self.duty_cycle_frac / (self.frequency_clock_khz * kHz))
        self.time_on_mu = self.core.seconds_to_mu((1 - self.duty_cycle_frac) / (self.frequency_clock_khz * kHz))
        self.num_repetitions = int(self.frequency_clock_khz * self.time_total_ms)

        # self.time_off_mu = self.core.seconds_to_mu(self.time_reset_us * us)

        # tmp remove
        self.dds = self.get_device('urukul0_ch1')
        self.dds_sw = self.get_device('ttl22')


    @kernel(flags={"fast-math"})
    def run(self):
        """
        Run the experimental sequence.
        """
        self.core.break_realtime()

        # prepare devices
        self.ttl_trigger.off()

        # record dma
        # todo

        # MAIN LOOP
        for i in range(self.num_repetitions):

            # OUTPUT ON
            with parallel:

                # toggle devices
                with sequential:
                    self.ttl_trigger.on()
                    self.dds.cpld.cfg_sw(0b0100)

                # on time
                delay_mu(self.time_delay_mu)

            # OUTPUT OFF
            with parallel:

                # toggle devices
                with sequential:
                    self.ttl_trigger.off()
                    self.dds.cpld.cfg_sw(0b0000)

                # off time
                delay_mu(self.time_delay_mu)

    @kernel
    def configure_ram_mode(self):
        self.core.break_realtime()

        self.dds.set_cfr1(ram_enable=0)
        self.cpld.io_update.pulse_mu(8)
        self.cpld.set_profile(0)  # Enable the corresponding RAM profile
        # Profile 0 is the default
        self.dds.set_profile_ram(start=0, end=len(self.asf_ram) - 1,
                            step=250, profile=0, mode=RAM_MODE_CONT_RAMPUP)
        self.cpld.io_update.pulse_mu(8)
        self.dds.amplitude_to_ram(self.amp, self.asf_ram)
        self.dds.write_ram(self.asf_ram)
        self.core.break_realtime()
        self.dds.set(frequency=5 * MHz, ram_destination=RAM_DEST_ASF)
        # Pass osk_enable=1 to set_cfr1() if it is not an amplitude RAM
        self.dds.set_cfr1(ram_enable=1, ram_destination=RAM_DEST_ASF)
        self.cpld.io_update.pulse_mu(8)
