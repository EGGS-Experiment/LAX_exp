import extensions
import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment


class IonLoad(LAXExperiment, Experiment):
    """
    Experiment: Ion Load

    Gets the number of counts as a function of frequency for a fixed time.
    """
    name = 'Ion Load'

    def build_experiment(self):

        # laser arguments
        self.setattr_argument('freq_397',
                              NumberValue(default=110., ndecimals=1, step=0.1, min=90., max=120., unit="MHz"), group='397')
        self.setattr_argument('ampl_397', NumberValue(default=50., ndecimals=1, step=0.1, min=0., max=100.), group='397')
        self.setattr_argument('att_397', NumberValue(default=14., ndecimals=1, step=0.1, min=0., max=31.5, unit="dB"), group='397')

        self.setattr_argument('freq_866',
                              NumberValue(default=110., ndecimals=1, step=0.1, min=90., max=120., unit="MHz"), group='866')
        self.setattr_argument('ampl_866', NumberValue(default=50., ndecimals=1, step=0.1, min=0., max=100.), group='866')
        self.setattr_argument('att_866', NumberValue(default=14., ndecimals=1, step=0.1, min=0., max=31.5, unit="dB"), group='866')

        # pmt arguments
        self.setattr_argument('ion_count_threshold', NumberValue(default=120, ndecimals=0, step=1, min=80, max=250), group='phonton_counting')
        self.setattr_argument('pmt_sample_time_us',
                              NumberValue(default=3e-3, ndecimals=0, step=1, min=1e-6, max=5e-3, unit='us'), group='phonton_counting')

        # relevant devices
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')
        self.setattr_device('shutters')
        self.setattr_device('oven')
        self.setattr_device('pmt')
        self.setattr_device('scheduler')

    def prepare_experiment(self):

        # convert 397 parameters
        self.ftw_397 = extensions.mhz_to_ftw(self.freq_397)
        self.asf_397 = extensions.pct_to_asf(self.ampl_397)
        self.att_397 = extensions.att_to_mu(self.att_397)

        # convert 866 parameters
        self.ftw_866 = extensions.mhz_to_ftw(self.freq_866)
        self.asf_866 = extensions.pct_to_asf(self.ampl_866)
        self.att_866 = extensions.att_to_mu(self.att_866)

        # open shutters
        self.shutters.open_377_shutter()
        self.shutters.open_423_shutter()

        # turn on oven
        self.gpp3060.turn_oven_on()

    @property
    def results_shape(self):
        return None

    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self):
        self.core.break_realtime()

        # set 397 parameters
        self.pump.beam.set_ftw(self.ftw_397)
        self.pump.beam.set_asf(self.asf_397)
        self.pump.beam.set_att_mu(self.att_397)

        # set 866 parameters
        self.repump_cooling.beam.set_ftw(self.ftw_866)
        self.repump_cooling.beam.set_asf(self.asf_866)
        self.repump_cooling.beam.set_att_mu(self.att_866)

        # set pmt sample time
        self.pmt.count(self.pmt_sample_time_us)

        self.start_time = now_mu()

    @kernel(flags={"fast-math"})
    def run_main(self):
        self.core.break_realtime()  # add slack

        counts = 0
        idx = 0
        breaker = False

        while counts < self.ion_count_threshold:
            self.core.break_realtime()  # add slack
            counts = self.pmt.fetch_count()  # grab counts from PMT
            self.check_termination()  # check if termination is over threshold
            idx += 1

            if idx >= 49:
                idx = 0
                breaker = self.check_time(now_mu())

            if breaker: break

    # ANALYSIS
    def analyze_experiment(self):

        self.cleanup_devices()

    @rpc
    def cleanup_devices(self):
        # turn off oven
        self.gpp3060.turn_oven_off()

        # close shutters
        self.shutters.open_377_shutter()
        self.shutters.open_423_shutter()

        # disconnect from labjack
        self.shutters.close_labjack()

    @rpc(flags={"async"}) # does this need to be async 
    def check_time(self, time):
        return 600 > self.core.mu_to_seconds(time - self.start_time)  # check if longer than 10 min

    @rpc
    def check_termination(self):
        """
        Check whether termination of the experiment has been requested.
        """
        if self.scheduler.check_termination():
            self.cleanup_devices()
            raise TerminationRequested
