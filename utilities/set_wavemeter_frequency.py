from artiq.experiment import *
from LAX_exp.base import LAXExperiment

class SetWavemeterFrequency(LAXExperiment, Experiment):

    def build_experiment(self):

        self.setattr_argument("wavemeter_channel", NumberValue(default=5))

        self.setattr_argument("wavemeter_frequency_thz", NumberValue(default=755.2217))

        self.setattr_device("wavemeter")
        self.setattr_device("core")
        self.setattr_device("scheduler")

    def prepare_experiment(self):


        ### safety checks ###
        # check if wavemeter channel is available
        if self.wavemeter_channel not in [5,4, 9, 14, 13]:
            raise ValueError("Please insert a proper channel - below are the channels available: \n"
                             "Channel 4: 423nm \n"
                             "Channel 5: 397nm \n"
                             "Channel 9: 729nm \n"
                             "Channel 13: 866nm \n"
                             "Channel 14: 854nm")

        # check if wavemeter frequency is reasonable for a given channel:
        if self.wavemeter_channel == 4 and (709.0776 < self.wavemeter_frequency_thz
                                            or self.wavemeter_frequency_thz > 709.0777):
            raise ValueError("Wavemeter frequency for channel 4 (423 nm) must be between 755.2217 and 755.2219 THz")


        if self.wavemeter_channel == 5 and (755.2217 < self.wavemeter_frequency_thz
                                            or self.wavemeter_frequency_thz > 755.2219):
            raise ValueError("Wavemeter frequency for channel 5 (397nm) must be between 755.2217 and 755.2219 THz")


        if self.wavemeter_channel == 9:
            raise ValueError("Do not change wavemeter value of channel 9 (729nm)")

        if self.wavemeter_channel == 13 and (345.9998 < self.wavemeter_frequency_thz
                                            or self.wavemeter_frequency_thz > 346):
            raise ValueError("Wavemeter frequency for channel 13 (866nm) must be between 345.9998 and 346 THz")

        if self.wavemeter_channel == 14 and (350.86 < self.wavemeter_frequency_thz
                                             or self.wavemeter_frequency_thz > 250.87):
            raise ValueError("Wavemeter frequency for channel 14 (854nm) must be between 345.9998 and 346 THz")

    def initialize_experiment(self):
        pass

    @property
    def results_shape(self):
        return (2,2)

    def run_main(self):

        status = self.scheduler.get_status()
        exp_status = [status[rid]['status'] for rid in status]

        # prevents experiment from running if other experiments are still analyzing data
        while 'analyzing' not in exp_status:
            exp_status = [status[rid]['status'] for rid in status]
            if self.scheduler.check_termination():
                break

        if self.scheduler.check_termination():
            raise TerminationRequested

    def analyze_experiment(self):
        pass

