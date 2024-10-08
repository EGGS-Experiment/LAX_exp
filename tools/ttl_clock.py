from artiq.experiment import *

# todo: duty cycle


class TTLClock(EnvExperiment):
    """
    TTL Clock
    Use the DIO-BNC TTLs to create a clock pulse.
    """

    def build(self):
        """
        Set devices and arguments for the experiment.
        """
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # timing
        self.setattr_argument('frequency_clock_khz',            NumberValue(default=1, ndecimals=3, step=1, min=0.001, max=1000))
        self.setattr_argument('time_total_ms',                  NumberValue(default=10000, ndecimals=3, step=1, min=0.001, max=100000))
        #self.setattr_argument('time_off_us', NumberValue(default=100, ndecimals=3, step=1, min=0, max=100))

        # TTL channel
        self.setattr_argument("ttl_channel",                    NumberValue(default=9, ndecimals=0, step=1, min=4, max=23))


    def prepare(self):
        """
        Set up the dataset and prepare values to minimize kernel overhead.
        """
        # devices
        self.ttl_clock =                self.get_device("ttl{:d}".format(self.ttl_channel))

        # timing
        self.time_trigger_delay_mu =    self.core.seconds_to_mu(50.88 * us)
        self.time_delay_mu =            self.core.seconds_to_mu(0.5 / (self.frequency_clock_khz * kHz))
        self.num_repetitions =          int(self.frequency_clock_khz * self.time_total_ms)
        # self.time_off_mu = self.core.seconds_to_mu(self.time_reset_us * us)


    @kernel
    def run(self):
        """
        Run the experimental sequence.
        """
        # reset hardware
        self.core.wait_until_mu(now_mu())
        self.core.reset()
        self.ttl_clock.off()

        # MAIN LOOP
        for i in range(self.num_repetitions):
            # TTL ON
            with parallel:
                self.ttl_clock.on()
                delay_mu(self.time_delay_mu)

            # TTL OFF
            with parallel:
                self.ttl_clock.off()
                delay_mu(self.time_delay_mu)
