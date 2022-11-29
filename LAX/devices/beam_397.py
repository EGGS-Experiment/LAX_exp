from artiq.experiment import *

DDS_CHANNEL = 'urukul1_'


class beam_397_pump(HasEnvironment):
    """
    A beam
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
