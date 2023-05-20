import labrad
from os import environ
from EGGS_labrad.config.dc_config import dc_config

import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import ParametricExcite
# from LAX_exp.analysis.signals import demodulateCounts


class ParametricSweep(LAXExperiment, Experiment):
    """
    Experiment: Parametric Sweep

    Modulate the trap RF close to a secular frequency while sweeping shim voltgaes
    to measure micromotion.
    """
    name = 'Parametric Sweep'


    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions",                        NumberValue(default=1, ndecimals=0, step=1, min=1, max=10000))
        self.setattr_argument("num_counts",                         NumberValue(default=1000, ndecimals=0, step=1, min=1, max=10000000))

        # modulation
        self.setattr_argument("mod_att_db",                         NumberValue(default=31, ndecimals=1, step=0.5, min=0, max=31.5), group='modulation')
        self.setattr_argument("mod_freq_khz_list",                  Scannable(
                                                                        default=CenterScan(1207, 8, 0.2, randomize=True),
                                                                        global_min=1, global_max=200000, global_step=1,
                                                                        unit="kHz", scale=1, ndecimals=4
                                                                    ), group='modulation')

        # shimming voltages
        self.dc_micromotion_channeldict =                           dc_config.channeldict
        self.setattr_argument("dc_micromotion_channel",             EnumerationValue(list(self.dc_micromotion_channeldict.keys()), default='V Shim'), group='voltage')
        self.setattr_argument("dc_micromotion_voltages_v_list",     Scannable(
                                                                        default=ExplicitScan([30.5]),
                                                                        global_min=0, global_max=400, global_step=1,
                                                                        unit="V", scale=1, ndecimals=4
                                                                    ), group='voltage')

        # cooling
        self.setattr_argument("ampl_cooling_pct",                   NumberValue(default=50, ndecimals=2, step=5, min=0.01, max=50), group='cooling')
        self.setattr_argument("freq_cooling_mhz",                   NumberValue(default=105, ndecimals=6, step=1, min=1, max=500), group='cooling')

        # get devices
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')
        self.setattr_device('repump_qubit')
        self.setattr_device('dds_modulation')

        # subsequences
        self.parametric_subsequence =                               ParametricExcite(self)


    def prepare_experiment(self):
        # get voltage parameters
        self.dc_micromotion_channel_num =                           self.dc_micromotion_channeldict[self.dc_micromotion_channel]['num']
        self.dc_micromotion_voltages_v_list =                       np.array(list(self.dc_micromotion_voltages_v_list))

        # convert cooling parameters to machine units
        self.ampl_cooling_asf =                                     self.pump.amplitude_to_asf(self.ampl_cooling_pct / 100)
        self.freq_cooling_ftw =                                     self.pump.frequency_to_ftw(self.freq_cooling_mhz * MHz)
        self.att_cooling_mu =                                       att_to_mu(14 * dB)
        self.time_cooling_holdoff_mu =                              self.core.seconds_to_mu(3 * ms)

        # modulation control and synchronization
        self.att_modulation_mu =                                    att_to_mu(self.mod_att_db * dB)
        self.freq_modulation_list_mu =                              np.array([
                                                                        self.dds_modulation.frequency_to_ftw(freq_mhz * kHz)
                                                                        for freq_mhz in self.mod_freq_khz_list
                                                                    ])

        # connect to labrad
        self.cxn =                                                  labrad.connect(environ['LABRADHOST'], port=7682, tls_mode='off', username='', password='lab')
        self.dc =                                                   self.cxn.dc_server

    @property
    def results_shape(self):
        return (self.repetitions * len(self.dc_micromotion_voltages_v_list) * len(self.mod_freq_khz_list),
                5)


    # LABRAD FUNCTIONS
    @rpc
    def voltage_set(self, channel: TInt32, voltage_v: TFloat):
        """
        Set the channel to the desired voltage.
        """
        # set desired voltage
        voltage_set_v = self.dc.voltage_fast(channel, voltage_v)
        print('\tvoltage set: {}'.format(voltage_set_v))

    @rpc
    def prepareDevicesLabrad(self):
        """
        Prepare LabRAD devices for the experiment via RPC.
        """
        # set up amo8
        self.dc.polling(False)
        self.dc.alarm(False)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self):
        self.core.break_realtime()

        with parallel:
            # set cooling beams
            with sequential:
                self.pump.set_mu(self.freq_cooling_ftw, asf=self.ampl_cooling_asf)
                self.pump.on()
                self.repump_cooling.on()
                self.repump_qubit.on()

            # set modulation attenuation
            self.dds_modulation.set_att_mu(self.att_modulation_mu)

            # set up labrad devices via RPC
            self.prepareDevicesLabrad()


    @kernel(flags={"fast-math"})
    def run_main(self):
        self.core.reset()

        # run given number of repetitions
        for trial_num in range(self.repetitions):

            # sweep voltage
            for voltage_v in self.dc_micromotion_voltages_v_list:

                # set DC voltage
                self.voltage_set(self.dc_micromotion_channel_num, voltage_v)
                self.core.break_realtime()

                # sweep modulation frequencies
                for freq_mu in self.freq_modulation_list_mu:

                    # add holdoff period for recooling the ion
                    at_mu(now_mu() + self.time_cooling_holdoff_mu)

                    # set modulation frequency and run parametric excitation
                    self.dds_modulation.set_mu(freq_mu, asf=self.dds_modulation.ampl_modulation_asf)
                    pmt_timestamp_list = self.parametric_subsequence.run(self.num_counts)

                    # process results (stores them in our results dataset for us)
                    self._process_results(freq_mu, voltage_v, pmt_timestamp_list)
                    self.core.break_realtime()

    @rpc(flags={"async"})
    def _process_results(self, freq_mu: TInt32, voltage_v: TFloat, timestamp_mu_list: TArray(TInt64, 1)):
        """
        Convert modulation frequency and timestamps from machine units and demodulate.

        Arguments:
            freq_mu             (int)           : the modulation frequency (in machine units).
            voltage_v           (float)         : the current shim voltage (in volts).
            timestamp_mu_list   (list(int64))   : the list of timestamps (in machine units) to demodulate.
        """
        # convert frequency to mhz
        freq_mhz = self.dds_modulation.ftw_to_frequency(freq_mu) / MHz

        # remove starting time and digitally demodulate counts
        timestamps_s = self.core.mu_to_seconds(np.array(timestamp_mu_list))
        correlated_signal = np.mean(np.exp((2.j * np.pi * freq_mhz * 1e6) * timestamps_s))

        # get count rate in seconds
        count_rate_hz = len(timestamps_s) / (timestamps_s[-1] - timestamps_s[0])

        # update dataset
        self.update_dataset(freq_mhz, voltage_v, np.abs(correlated_signal), np.angle(correlated_signal), count_rate_hz)
