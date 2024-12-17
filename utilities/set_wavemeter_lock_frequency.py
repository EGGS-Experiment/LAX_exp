from artiq.experiment import *
from LAX_exp.base import LAXExperiment
from EGGS_labrad.config.multiplexerclient_config import multiplexer_config

class SetWavemeterFrequency(LAXExperiment, Experiment):

    """
    Experiment: Set Wavemeter Lock Frequency

    Set the lock frequency for a wavemeter channel
    """

    DATASET_KEY = 'calibration/wavemeter'

    def build_experiment(self):

        self.setattr_argument("wavemeter_channel", NumberValue(default=5))

        self.setattr_argument("wavemeter_lock_frequency_thz", NumberValue(default=755.2217))

        self.setattr_argument("check_dataset", BooleanValue(False))

        self.setattr_device("wavemeter")
        self.setattr_device("core")
        self.setattr_device("scheduler")

    def prepare_experiment(self):

        available_channels = [channels[channel_key][0] for channel_key in list(channels.keys())]
        ### safety checks ###
        # check if wavemeter channel is available
        if self.wavemeter_channel not in available_channels:
            raise ValueError("Please insert a proper channel - below are the channels available: \n"
                             "Channel 4: 423nm \n"
                             "Channel 5: 397nm \n"
                             "Channel 9: 729nm \n"
                             "Channel 13: 866nm \n"
                             "Channel 14: 854nm")

        # check if wavemeter frequency is reasonable for a given channel:
        if self.wavemeter_channel == 4:
            if self.check_dataset:
                self.DATASET_KEY += '/423'

            if (709.0776 < self.wavemeter_frequency_thz or self.wavemeter_frequency_thz > 709.0777):
                raise ValueError("Wavemeter frequency for channel 4 (423 nm) must be between 755.2217 and 755.2219 THz")


        if self.wavemeter_channel == 5:
            if self.check_dataset:
                self.DATASET_KEY += '/397'

            if (755.2217 < self.wavemeter_frequency_thz or self.wavemeter_frequency_thz > 755.2219):
                raise ValueError("Wavemeter frequency for channel 5 (397nm) must be between 755.2217 and 755.2219 THz")


        if self.wavemeter_channel == 9:
            raise ValueError("Do not change wavemeter value of channel 9 (729nm)")

        if self.wavemeter_channel == 13:

            if self.check_dataset:
                self.DATASET_KEY += '/866'

            if (345.9998 < self.wavemeter_frequency_thz or self.wavemeter_frequency_thz > 346):
                raise ValueError("Wavemeter frequency for channel 13 (866nm) must be between 345.9998 and 346 THz")

        if self.wavemeter_channel == 14:
            if self.check_dataset:
                self.DATASET_KEY += '/854'

            if (350.86 < self.wavemeter_frequency_thz or self.wavemeter_frequency_thz > 250.87):
                raise ValueError("Wavemeter frequency for channel 14 (854nm) must be between 345.9998 and 346 THz")


        if self.check_dataset:
            self.lock_frequency_thz = self.get_dataset(self.DATASET_KEY)
        else:
            self.lock_frequency_thz = self.wavemeter_lock_frequency_thz

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

        # self.wavemeter.set_channel_lock_frequency(self.wavemeter_channel, self.lock_frequency_thz)
    def analyze_experiment(self):
        pass

