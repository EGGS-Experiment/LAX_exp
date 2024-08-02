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
        self.setattr_argument('freq_397_mhz',
                              NumberValue(default=110., ndecimals=1,
                                          step=0.1, min=90., max=120., unit="MHz"), group='397')
        self.setattr_argument('ampl_397', NumberValue(default=50., ndecimals=1,
                                                      step=0.1, min=0., max=100.), group='397')
        self.setattr_argument('att_397_dB', NumberValue(default=14., ndecimals=1,
                                                        step=0.1, min=0., max=31.5, unit="dB"), group='397')

        self.setattr_argument('freq_866_mhz', NumberValue(default=110., ndecimals=1,
                                                          step=0.1, min=90., max=120., unit="MHz"), group='866')
        self.setattr_argument('ampl_866', NumberValue(default=50., ndecimals=1,
                                                      step=0.1, min=0., max=100.), group='866')
        self.setattr_argument('att_866_dB', NumberValue(default=14., ndecimals=1,
                                                        step=0.1, min=0., max=31.5, unit="dB"), group='866')

        self.setattr_argument('freq_854_mhz',
                              NumberValue(default=110., ndecimals=1,
                                          step=0.1, min=90., max=120., unit="MHz"), group='854')
        self.setattr_argument('ampl_854', NumberValue(default=50., ndecimals=1,
                                                      step=0.1, min=0., max=100.), group='854')
        self.setattr_argument('att_854_dB', NumberValue(default=14., ndecimals=1,
                                                        step=0.1, min=0., max=31.5, unit="dB"), group='854')

        # pmt arguments
        self.setattr_argument('ion_count_threshold', NumberValue(default=120, ndecimals=0,
                                                                 step=1, min=80, max=250), group='Photon Counting')
        self.setattr_argument('pmt_sample_time_us',
                              NumberValue(default=3e-3, ndecimals=0,
                                          step=1, min=1e-6, max=5e-3, unit='us'), group='Photon Counting')

        # starting trap arguments
        self.setattr_argument('starting_east_endcap_voltage',
                              NumberValue(default=19., ndecimals=1, step=0.1, min=0., max=300.),
                              group='Starting Trap Parameters')
        self.setattr_argument('starting_west_endcap_voltage',
                              NumberValue(default=25., ndecimals=1, step=0.1, min=0., max=300.),
                              group='Starting Trap Parameters')

        # starting trap arguments
        self.setattr_argument('ending_east_endcap_voltage',
                              NumberValue(default=19., ndecimals=1, step=0.1, min=0., max=300.),
                              group='Ending Trap Parameters')
        self.setattr_argument('ending_west_endcap_voltage',
                              NumberValue(default=25., ndecimals=1, step=0.1, min=0., max=300.),
                              group='Ending Trap Parameters')
        self.setattr_argument('ending_v_shim_voltage',
                              NumberValue(default=60., ndecimals=1, step=0.1, min=0., max=150.),
                              group='Ending Trap Parameters')
        self.setattr_argument('ending_h_shim_voltage',
                              NumberValue(default=51.5, ndecimals=1, step=0.1, min=0., max=150.),
                              group='Ending Trap Parameters')
        self.setattr_argument('ending_a_ramp2_voltage',
                              NumberValue(default=2., ndecimals=1, step=0.1, min=0., max=100.),
                              group='Ending Trap Parameters')

        # relevant devices
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')
        self.setattr_device('repump_qubit')
        self.setattr_device('shutters')
        self.setattr_device('gpp3060')
        self.setattr_device('pmt')
        self.setattr_device('aperture')
        self.setattr_device('trap_dc')
        self.setattr_device('scheduler')

    def prepare_experiment(self):

        # convert 397 parameters
        self.ftw_397 = extensions.mhz_to_ftw(self.freq_397_mhz)
        self.asf_397 = extensions.pct_to_asf(self.ampl_397)
        self.att_397 = extensions.att_to_mu(self.att_397_dB)

        # convert 854 parameters
        self.ftw_854 = extensions.mhz_to_ftw(self.freq_854_mhz)
        self.asf_854 = extensions.pct_to_asf(self.ampl_854)
        self.att_854 = extensions.att_to_mu(self.att_854_dB)

        # convert 866 parameters
        self.ftw_866 = extensions.mhz_to_ftw(self.freq_866_mhz)
        self.asf_866 = extensions.pct_to_asf(self.ampl_866)
        self.att_866 = extensions.att_to_mu(self.att_866_dB)

    @property
    def results_shape(self):
        return (2, 2)

    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self):
        self.core.break_realtime()
        self.pump.beam.cpld.get_att_mu()

        # set 397 parameters
        self.pump.set_mu(self.ftw_397, asf=self.asf_397, profile=6)
        self.pump.set_att_mu(self.att_397)

        # set 854 parameters
        self.repump_qubit.set_mu(self.ftw_854, asf=self.asf_854, profile=6)
        self.repump_qubit.set_att_mu(self.att_854)

        # set 866 parameters
        self.repump_cooling.set_mu(self.ftw_866, asf=self.asf_866, profile=6)
        self.repump_cooling.set_att_mu(self.att_866)

        # turn on lasers
        self.pump.on()
        self.repump_qubit.on()
        self.repump_cooling.on()
        self.core.break_realtime()

        # set endcaps to loading voltages
        self.trap_dc.set_east_endcap_voltage(self.starting_east_endcap_voltage)
        self.trap_dc.set_west_endcap_voltage(self.starting_west_endcap_voltage)

        # turn on endcap channels and ensure others are off
        self.trap_dc.east_endcap_on()
        self.trap_dc.west_endcap_on()
        self.trap_dc.h_shim_off()
        self.trap_dc.v_shim_off()
        self.trap_dc.a_ramp2_off()
        self.core.break_realtime()

        # open shutters
        self.shutters.open_377_shutter()
        self.shutters.open_423_shutter()
        self.core.break_realtime()

        self.gpp3060.turn_oven_on()  # turn on oven
        self.core.break_realtime()

        self.aperture.open_aperture()
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def run_main(self):
        self.core.break_realtime()  # add slack
        self.start_time = self.get_rtio_counter_mu()

        counts = 0
        idx = 0
        breaker = False
        count_successes = 0

        while count_successes < 10:
            self.core.break_realtime()  # add slack
            self.pmt.count(self.pmt_sample_time_us)  # set pmt sample time
            counts = self.pmt.fetch_count()  # grab counts from PMT
            self.core.break_realtime()

            if counts >= self.ion_count_threshold:
                count_successes += 1

            idx += 1

            if idx >= 49:
                idx = 0
                breaker = self.check_time()
                self.check_termination()  # check if termination is over threshold
                self.core.break_realtime()

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

        # close aperture
        self.aperture.close_aperture()

        # set trap parameters as if ion was loaded
        self.trap_dc.set_east_endcap_voltage(self.starting_east_endcap_voltage)
        self.trap_dc.set_west_endcap_voltage(self.starting_west_endcap_voltage)
        self.trap_dc.set_h_shim_voltage(self.ending_h_shim_voltage)
        self.trap_dc.set_v_shim_voltage(self.ending_v_shim_voltage)
        self.trap_dc.set_a_ramp2_voltage(self.ending_a_ramp2_voltage)

        # turn on the endcap channels
        self.trap_dc.east_endcap_on()
        self.trap_dc.west_endcap_on()
        self.trap_dc.h_shim_on()
        self.trap_dc.v_shim_on()
        self.trap_dc.a_ramp2_on()

        # ramp endcaps to end values
        self.trap_dc.ramp_both_endcaps([self.ending_east_endcap_voltage, self.ending_west_endcap_voltage],
                                       [100, 100])

    @kernel
    def check_time(self) -> TBool:
        self.core.break_realtime()
        return 600 > self.core.mu_to_seconds(
            self.get_rtio_counter_mu() - self.start_time)  # check if longer than 10 min

    @rpc
    def check_termination(self):
        """
        Check whether termination of the experiment has been requested.
        """
        if self.scheduler.check_termination():
            self.cleanup_devices()
            raise TerminationRequested
