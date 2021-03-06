from artiq.experiment import *


def input_led_state() -> TBool:
    return input("LED state: ") == "1"


class LED(EnvExperiment):
    """
    LED
    Switch on some LEDs according to input.
    """

    def build(self):
        self.setattr_device("core")
        self.setattr_device("zotino0")
        self.setattr_device("ttl4")

    @kernel
    def run(self):
        self.core.reset()
        self.ttl4.output()
        self.ttl4.on()
        led_state = input_led_state()
        self.core.break_realtime()
        if led_state:
            self.ttl4.on()
            self.zotino0.set_leds(0b11111111)
        else:
            self.ttl4.off()
            self.zotino0.set_leds(0b00000000)

