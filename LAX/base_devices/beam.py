from artiq.experiment import *


class Beam_Urukul(HasEnvironment):
    """
    A generic "beam" object based off an Urukul DDS channel.
    """

    def build(self):
        self.setattr_argument("sc1_scan",
            Scannable(default=[NoScan(3250), RangeScan(10, 20, 6, randomize=True)],
                      unit="kHz"),
            "Flux capacitor")
        self.setattr_argument("sc1_enum", EnumerationValue(["1", "2", "3"]),
                              "Flux capacitor")

    def do(self):
        print("SC1:")
        for i in self.sc1_scan:
            print(i)
        print(self.sc1_enum)
