import numpy as np
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import DopplerCool

import labrad
from os import environ
from datetime import datetime




class ReadoutCalibration(LAXExperiment, Experiment):
    """
    Calibration: Readout

    Get amplitude scaling factors to compensate for frequency dependence.
    """
    # kernel_invariants = {}
    name = 'Readout Calibration'


    def build_experiment(self):
        """
        Set devices and arguments for the experiment.
        """
        # search parameters
        self.setattr_argument("repetitions",                                NumberValue(default=100, ndecimals=0, step=10, min=1, max=10000))

        # readout scan parameters
        self.setattr_argument("time_readout_us",                            NumberValue(default=3000, ndecimals=3, step=100, min=1, max=100000), group='Readout')
        self.setattr_argument("freq_readout_mhz_list",                      Scannable(
                                                                                default=[
                                                                                    RangeScan(100, 110, 21, randomize=True),
                                                                                    ExplicitScan([105])
                                                                                ],
                                                                                global_min=70, global_max=140, global_step=1,
                                                                                unit="MHz", scale=1, ndecimals=5
                                                                            ), group='Readout')
        self.setattr_argument("ampl_readout_pct_list",                      Scannable(
                                                                                default=[
                                                                                    RangeScan(20, 50, 61, randomize=True),
                                                                                    ExplicitScan([45])
                                                                                ],
                                                                                global_min=1, global_max=50, global_step=1,
                                                                                unit="pct", scale=1, ndecimals=1
                                                                            ), group='Readout')

        # relevant devices
        self.setattr_device('pmt')
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')

        # get relevant subsequences
        self.dopplercool_subsequence =                                      DopplerCool(self)

    def prepare_experiment(self):
        """
        Prepare things such that kernel functions have minimal overhead.
        """
        # convert DDS values
        self.time_readout_mu =                                              self.core.seconds_to_mu(self.time_readout_us * us)
        self.freq_readout_ftw_list =                                        np.array([self.pump.frequency_to_ftw(freq_mhz * MHz) for freq_mhz in list(self.freq_readout_mhz_list)])
        self.ampl_readout_asf_list =                                        np.array([self.pump.amplitude_to_asf(ampl_pct / 100.) for ampl_pct in list(self.ampl_readout_pct_list)])

        # create experimental config
        self.config_readout_calibration_list =                              np.stack(np.meshgrid(self.freq_readout_ftw_list, self.ampl_readout_asf_list), -1).reshape(-1, 2)
        np.random.shuffle(self.config_readout_calibration_list)

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_readout_calibration_list),
                4)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self):
        self.core.break_realtime()

        # record calibration sequence onto DMA
        with self.core_dma.record("PMT_READOUT_CALIB"):

            # run doppler cooling
            self.dopplercool_subsequence.run()

            # run readout
            self.pump.readout()
            self.pump.on()
            self.pmt.count(self.time_readout_mu)

            # get background counts
            self.repump_cooling.off()
            self.pmt.count(self.time_readout_mu)

    @kernel(flags={"fast-math"})
    def run_main(self):
        """
        Run the experimental sequence.
        """
        self.core.reset()

        # retrieve calibration sequence from DMA
        handle_dma = self.core_dma.get_handle("PMT_READOUT_CALIB")
        self.core.break_realtime()

        for trial_num in range(self.repetitions):

            # sweep experiment config: readout frequency and amplitude
            for config_vals in self.config_readout_calibration_list:

                # extract values from config list
                freq_readout_ftw = config_vals[0]
                ampl_readout_asf = config_vals[1]
                self.core.break_realtime()

                # set up readout profile for DDS
                self.pump.set_mu(freq_readout_ftw, asf=ampl_readout_asf, profile=1)

                # run calibration sequence
                self.core_dma.playback_handle(handle_dma)

                # retrieve PMT counts and update datasets
                counts_signal = self.pmt.fetch_count()
                counts_background = self.pmt.fetch_count()
                with parallel:
                    self.update_results(freq_readout_ftw, ampl_readout_asf, counts_signal, counts_background)
                    self.core.break_realtime()


    # ANALYSIS
    def analyze(self):
        """
        Analyze the results from the experiment.
        """
        print('idk')
        # # set calibration timestamp
        # calib_timestamp = datetime.timestamp(datetime.now())
        # self.set_dataset('calibration.readout.calibration_timestamp', calib_timestamp, broadcast=True, persist=True)
        #
        # # sort results
        # _indices_sorted = np.argsort(self.results, axis=0)[:, 0]
        # results_tmp = self.results[_indices_sorted].transpose()
        #
        # # convert results to desired output form
        # calib_freq_ftw, calib_ampl_asf = results_tmp
        # calib_freq_mhz = np.array(self.dds.ftw_to_frequency(calib_freq_ftw)) / MHz
        # calib_ampl_frac = np.array(self.dds.asf_to_amplitude(calib_ampl_asf)) * 100
        # calib_final = np.array([calib_freq_mhz, calib_ampl_frac]).transpose()
        #
        # # print run time
        # time_stop = datetime.timestamp(datetime.now())
        # print()
        # print('\t\t\tDONE')
        # print('\t\t\t\tTOTAL RUN TIME: {:.2f}'.format(time_stop-self.time_start))
        # print()
        #
        # # add calibration values to dataset manager
        # self.set_dataset('calibration.temperature.asf_calibration_curve_mhz_pct', calib_final, broadcast=True, persist=True)
        #
        #
        # ### LABRAD UPLOAD ###
        # # upload data to labrad for visualization in RealSimpleGrapher
        # try:
        #     # create connections to labrad servers
        #     cxn =                   labrad.connect(environ['LABRADHOST'], port=7682, tls_mode='off', username='', password='lab')
        #     dv =                    cxn.data_vault
        #     cr =                    cxn.context()
        #
        #     # create labrad dataset title
        #     date =                  datetime.now()
        #     dataset_title_tmp =     'DDS Amplitude Calibration'
        #     trunk =                 '{0:s}_{1:02d}:{2:02d}'.format(dataset_title_tmp, date.hour, date.minute)
        #     trunk_tmp =             ['', 'labrad', str(date.year), '{:02d}'.format(date.month), '{0:02d}'.format(date.day), trunk]
        #
        #     # create labrad dataset
        #     dv.cd(trunk_tmp, True, context=cr)
        #     dv.new(
        #         dataset_title_tmp,
        #         [('DDS Frequency', 'MHz')],
        #         [('DDS Amplitude', 'Amplitude', 'pct')],
        #         context=cr
        #     )
        #
        #     # upload data to labrad's Data Vault
        #     dv.add(calib_final, context=cr)
        #     print("\tLabRAD upload successful.")
        #
        # except Exception as e:
        #     print(e)
        #     print("Warning: unable to upload data to labrad.")
