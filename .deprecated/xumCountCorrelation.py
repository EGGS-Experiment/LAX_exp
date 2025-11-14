import labrad
from os import environ
from EGGS_labrad.config.dc_config import dc_config

import numpy as np
from artiq.experiment import *
from artiq.coredevice.ad9910 import PHASE_MODE_CONTINUOUS

from LAX_exp.analysis import *
from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import ParametricExcite, RescueIon


class InsufficientCounts(Exception):
    """
    Raised when the PMT is not getting sufficient counts for some reason.
    Prevents fixed PMT count-type experiments (e.g. ParametricSweep) from
    blocking due to insufficient PMT counts.
    """
    pass


class MicromotionCountCorrelation(LAXExperiment, Experiment):
    """
    Experiment: Micromotion Count Correlation

    Correlate readout counts with the trap RF to compensate axial micromotion.
    """
    name = 'Micromotion Count Correlation'
    kernel_invariants = {
        # hardware parameters
        "dc_channel_num", "dc_voltages_v_list", "time_dc_synchronize_delay_mu",
        "ampl_cooling_asf", "freq_cooling_ftw", "time_cooling_holdoff_mu",

        # labrad objects
        "cxn", "dc",

        # subsequences
        "rescue_subsequence"
    }

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions",            NumberValue(default=2, precision=0, step=1, min=1, max=10000))
        self.setattr_argument("freq_trap_rf_mhz",       NumberValue(default=19.05316, precision=6, step=0.01, min=0.01, max=100, scale=1., unit='MHz'))

        # shimming voltages
        self.dc_channeldict =                           dc_config.channeldict
        self.setattr_argument("dc_channel",             EnumerationValue(list(self.dc_channeldict.keys()), default='V Shim'), group='voltage')
        self.setattr_argument("dc_voltages_v_list",     Scannable(
                                                            default=[
                                                                CenterScan(74.4, 4., 0.6, randomize=True),
                                                                ExplicitScan([75.5]),
                                                            ],
                                                            global_min=0, global_max=400, global_step=1,
                                                            unit="V", scale=1, precision=1
                                                        ), group='voltage')

        # cooling
        self.setattr_argument("ampl_cooling_pct",       NumberValue(default=23, precision=2, step=5, min=0.01, max=50, scale=1., unit='%'), group='cooling')
        self.setattr_argument("freq_cooling_mhz",       NumberValue(default=105, precision=6, step=1, min=1, max=500, scale=1., unit='MHz'), group='cooling')

        # get relevant devices
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')
        self.setattr_device('repump_qubit')
        self.setattr_device('pmt')
        self.setattr_device('trigger_rf')

        # get relevant subsequences
        self.rescue_subsequence = RescueIon(self)

    def prepare_experiment(self):
        # get voltage parameters
        self.dc_channel_num =               self.dc_channeldict[self.dc_channel]['num']
        self.dc_voltages_v_list =           np.array(list(self.dc_voltages_v_list))
        self.time_dc_synchronize_delay_mu = self.core.seconds_to_mu(888 * ms)

        # convert cooling parameters to machine units
        self.ampl_cooling_asf =         self.pump.amplitude_to_asf(self.ampl_cooling_pct / 100)
        self.freq_cooling_ftw =         self.pump.frequency_to_ftw(self.freq_cooling_mhz * MHz)
        self.time_cooling_holdoff_mu =  self.core.seconds_to_mu(3 * ms)

        # connect to labrad
        self.cxn = labrad.connect(environ['LABRADHOST'],
                                  port=7682, tls_mode='off',
                                  username='', password='lab')
        self.dc = self.cxn.dc_server

    @property
    def results_shape(self):
        return (self.repetitions * len(self.dc_voltages_v_list),
                3)


    # LABRAD FUNCTIONS
    @rpc
    def voltage_set(self, channel: TInt32, voltage_v: TFloat) -> TNone:
        """
        Set the channel to the desired voltage.
        """
        # set desired voltage
        self.dc.voltage_fast(channel, voltage_v)

    @rpc(flags={"async"})
    def prepareDevicesLabrad(self) -> TNone:
        """
        Prepare LabRAD devices for the experiment via RPC.
        """
        # set up amo8
        self.dc.polling(False)
        self.dc.alarm(False)
        self.dc.serial_write('remote.w 1\r\n')
        self.dc.serial_read('\n')


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self) -> TNone:
        # set up labrad devices via RPC
        self.prepareDevicesLabrad()
        self.core.break_realtime()

        # set cooling beams
        self.pump.set_mu(self.freq_cooling_ftw, asf=self.ampl_cooling_asf, profile=0)
        self.pump.set_profile(0)
        self.pump.on()
        self.repump_cooling.on()
        self.repump_qubit.on()
        delay_mu(10000)

    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        # run given number of repetitions
        for trial_num in range(self.repetitions):

            # sweep voltage
            for voltage_v in self.dc_voltages_v_list:

                # set DC voltage and synchronize
                self.voltage_set(self.dc_channel_num, voltage_v)
                # synchronize hardware clock with timeline, then add delay for voltages to settle
                # note: delay has added advantage of recooling the ion
                self.core.wait_until_mu(now_mu())
                delay_mu(self.time_dc_synchronize_delay_mu)

                '''COUNT CORRELATION'''
                # synchronize to same phase of trap RF
                self.trigger_rf.trigger(self.time_rf_gating_mu, self.time_rf_holdoff_mu)

                # get timestamped photon counts
                pmt_timestamp_list = self.pmt.timestamp_counts(self.timestamp_mu_list, self.time_pmt_gating_mu)
                self.core.break_realtime()

                # process results (stores them in our results dataset for us)
                self._process_results(pmt_timestamp_list)
                self.core.reset()

            # rescue ion as needed
            self.rescue_subsequence.run(trial_num)

            # support graceful termination
            self.check_termination()
            self.core.break_realtime()


    @rpc(flags={"async"})
    def _process_results(self, voltage_v: TFloat, timestamp_mu_list: TArray(TInt64, 1)) -> TNone:
        """
        Convert modulation frequency and timestamps from machine units and demodulate.
        Arguments:
            voltage_v           (float)         : the current shim voltage (in volts).
            timestamp_mu_list   (list(int64))   : the list of timestamps (in machine units) to demodulate.
        """
        # todo: fix - need to process differently
        # convert timestamps and digitally demodulate counts
        timestamps_s = self.core.mu_to_seconds(np.array(timestamp_mu_list))
        correlated_signal = np.mean(np.exp((2.j * np.pi * self.freq_trap_rf_mhz * 1e6) * timestamps_s))
        # convert demodulated signal to polar coordinates (i.e. abs and angle)
        correlated_ampl = np.abs(correlated_signal)
        correlated_phase = np.angle(correlated_signal)

        # update dataset
        self.update_results(voltage_v, correlated_ampl, correlated_phase)


    # ANALYSIS
    def analyze_experiment(self):
        # todo
        pass
