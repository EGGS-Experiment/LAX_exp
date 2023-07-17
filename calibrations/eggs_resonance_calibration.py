import labrad
import numpy as np
from os import environ
from time import sleep
from datetime import datetime
from artiq.experiment import *


class EGGSResonanceCalibration(EnvExperiment):
    """
    Calibration: EGGS Resonance Calibration

    Examines the resonance for the EGGS feedthrough and calculates the appropriate amplitude scalings.
    """
    # kernel_invariants = {}

    def build(self):
        """
        Set devices and arguments for the experiment.
        """
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # experiment runs
        self.setattr_argument("repetitions",                                    NumberValue(default=5, ndecimals=0, step=1, min=1, max=10000))

        # eggs heating
        self.setattr_argument("freq_eggs_heating_mhz_list",                     Scannable(
                                                                                    default=RangeScan(75, 90, 151, randomize=False),
                                                                                    global_min=30, global_max=400, global_step=1,
                                                                                    unit="MHz", scale=1, ndecimals=5
                                                                                ))
        # spectrum analyzer
        self.setattr_argument("spectrum_analyzer_bandwidth_khz",                NumberValue(default=10, ndecimals=5, step=1, min=0.00001, max=10000))
        self.setattr_argument("spectrum_analyzer_attenuation_internal_db",      NumberValue(default=10, ndecimals=5, step=1, min=0.00001, max=10000))
        self.setattr_argument("spectrum_analyzer_attenuation_external_db",      NumberValue(default=0, ndecimals=5, step=1, min=0.00001, max=10000))


    def prepare(self):
        """
        Prepare things such that kernel functions have minimal overhead.
        """
        # eggs heating devices
        self.awg_board =                                                        self.get_device("phaser0")
        self.awg_eggs =                                                         self.awg_board.channel[0]

        # eggs heating timings
        self.time_phaser_sample_mu =                                            np.int64(40)

        # calculate eggs frequency
        self.freq_eggs_heating_mhz_list =                                       list(self.freq_eggs_heating_mhz_list)
        #tmp remove
        # self.freq_eggs_heating_mhz_list =                                       np.array([90, 87.5])
        #tmp remove clear
        freq_eggs_heating_center_mhz =                                          np.median(self.freq_eggs_heating_mhz_list)
        freq_eggs_heating_range_mhz =                                           np.max(self.freq_eggs_heating_mhz_list) - np.min(self.freq_eggs_heating_mhz_list)

        # connect to labrad
        self.cxn =                                                              labrad.connect(environ['LABRADHOST'], port=7682, tls_mode='off', username='', password='lab')
        self.sa =                                                               self.cxn.spectrum_analyzer_server
        self.dv =                                                               self.cxn.data_vault
        self.cr =                                                               self.cxn.context()

        # set up spectrum analyzer
        # get list of spectrum analyzers
        sa_dev_list = self.sa.list_devices()
        sa_dev_dict = dict(tuple(sa_dev_list))

        # select correct spectrum analyzer
        dev_exists = False
        for dev_num, dev_desc in sa_dev_dict.items():
            if 'MY500' in dev_desc:
                dev_exists = True
                self.sa.select_device(dev_num)

        # raise error if function generator doesn't exist
        if not dev_exists:
            raise Exception("Error: sexy spectrum analyzer not detected.")

        # set up spectrum analyzer sweep
        self.sa.attenuation(self.spectrum_analyzer_attenuation_internal_db)
        self.sa.frequency_span(1.25 * freq_eggs_heating_range_mhz * MHz)
        self.sa.frequency_center(freq_eggs_heating_center_mhz * MHz)
        self.sa.bandwidth_resolution(self.spectrum_analyzer_bandwidth_khz * 1000)
        # set up spectrum analyzer marker
        self.sa.peak_threshold(-90)
        self.sa.peak_excursion(15)
        self.sa.marker_toggle(1, True)
        # todo: move to actual labrad functions
        self.sa.gpib_write('DET:TRAC1 AVER')
        self.sa.gpib_write('DISP:ENAB 1')

        # set up data vault
        # create labrad dataset title
        date =                  datetime.now()
        dataset_title_tmp =     'EGGS Resonance Calibration'
        trunk =                 '{0:s}_{1:02d}:{2:02d}'.format(dataset_title_tmp, date.hour, date.minute)
        trunk_tmp =             ['', 'labrad', str(date.year), '{:02d}'.format(date.month), '{0:02d}'.format(date.day), trunk]

        # create labrad dataset
        self.dv.cd(trunk_tmp, True, context=self.cr)
        self.dv.new(
            dataset_title_tmp,
            [('Frequency', 'MHz')],
            [
                ('Transmission',            'Power',        'dBm'),
                ('Normalized Amplitude',    'Fraction',     'frac')
            ],
            context=self.cr
        )
        print("Data vault setup successful.")

        # set up artiq datasets
        self._iter_dataset = 0
        self.set_dataset("eggs_resonance_calibration", np.zeros([len(self.freq_eggs_heating_mhz_list), 2]))
        self.setattr_dataset("eggs_resonance_calibration")


    @kernel(flags={"fast-math"})
    def run(self):
        """
        Run the experimental sequence.
        """
        self.core.reset()

        # prepare phaser
        self.phaser_prepare()
        self.core.break_realtime()


        # MAIN SEQUENCE
        # for trial_num in range(self.repetitions):

        # sweep eggs rf frequencies
        for freq_mhz in self.freq_eggs_heating_mhz_list:

            # add extra delay
            self.core.break_realtime()

            # set eggs carrier via the DUC
            at_mu(self.awg_board.get_next_frame_mu())
            self.awg_eggs.set_duc_frequency((freq_mhz - 85) * MHz)
            self.awg_board.duc_stb()
            self.core.break_realtime()

            # wait given time
            self.record_power(freq_mhz)
            self.core.break_realtime()

        # disable rf output
        self.awg_eggs.en_trf_out(rf=0, lo=0)
        self.core.break_realtime()


    @kernel(flags={"fast-math"})
    def phaser_prepare(self):
        """
        Prepare phaser for the calibration.
        """
        # initialize phaser
        self.awg_board.init(debug=True)
        self.core.break_realtime()

        # set nco to center eggs rf around the center frequency (should be 85 MHz) exactly
        at_mu(self.awg_board.get_next_frame_mu())
        self.awg_eggs.set_nco_frequency(-217.083495 * MHz)
        self.awg_eggs.set_nco_phase(0.)
        self.awg_board.dac_sync()
        self.core.break_realtime()

        # trf setup, and disable rf output while we set things up
        self.awg_eggs.set_att(3 * dB)
        self.awg_eggs.en_trf_out(rf=0, lo=0)
        self.core.break_realtime()

        # duc
        self.awg_eggs.set_duc_frequency(0 * MHz)
        self.awg_eggs.set_duc_cfg()
        self.awg_board.duc_stb()
        self.core.break_realtime()

        # oscillators (i.e. sidebands)
        at_mu(self.awg_board.get_next_frame_mu())
        self.awg_eggs.oscillator[0].set_frequency(0)
        delay_mu(self.time_phaser_sample_mu)
        self.awg_eggs.oscillator[0].set_amplitude_phase(0.49, clr=0)
        delay_mu(self.time_phaser_sample_mu)
        self.awg_eggs.oscillator[1].set_amplitude_phase(0., clr=1)
        delay_mu(self.time_phaser_sample_mu)
        self.awg_eggs.oscillator[2].set_amplitude_phase(0., clr=1)
        delay_mu(self.time_phaser_sample_mu)
        self.awg_eggs.oscillator[3].set_amplitude_phase(0., clr=1)
        delay_mu(self.time_phaser_sample_mu)
        self.awg_eggs.oscillator[4].set_amplitude_phase(0., clr=1)
        self.core.break_realtime()

        # re-enable rf output
        self.awg_eggs.en_trf_out(rf=1, lo=0)
        self.core.break_realtime()


    @rpc
    def record_power(self, peak_freq_mhz):
        """
        Get and record the given number of peaks, sorted by amplitude.
        """
        # zoom in on band of interest
        self.sa.frequency_span(4 * self.spectrum_analyzer_bandwidth_khz * kHz)
        self.sa.frequency_center(peak_freq_mhz * MHz)

        # attempt to minimize measurement time
        sweep_time_s = float(self.sa.gpib_query('SWE:TIME?'))
        self.sa.gpib_write('CALC:MARK1:X ' + str(peak_freq_mhz * MHz))

        # zoom in on amplitude
        # todo

        # average multiple measurements
        peak_power_dbm_list = np.zeros(self.repetitions)
        num_vals = 0
        while num_vals < self.repetitions:

            # wait for spectrum analyzer sweep to finish before recording data point
            sleep(1.25 * sweep_time_s)
            power_dbm_tmp = self.sa.marker_amplitude(1)

            # ensure reading is valid
            if power_dbm_tmp < 1e10:
                peak_power_dbm_list[num_vals] = power_dbm_tmp
                num_vals += 1

        # average readings
        peak_power_dbm_avg = np.mean(peak_power_dbm_list)

        # save to dataset (artiq)
        self.mutate_dataset('eggs_resonance_calibration', self._iter_dataset, np.array([peak_freq_mhz, peak_power_dbm_avg]))
        self._iter_dataset += 1


    def analyze(self):
        """
        Analyze the results from the experiment.
        """
        # reenable spec anal display
        self.sa.gpib_write('DISP:ENAB 1')

        # store calibration timestamp
        calib_timestamp = datetime.timestamp(datetime.now())
        self.set_dataset('calibration.eggs.calibration_timestamp', calib_timestamp, broadcast=True, persist=True)

        # copy calibration results to temporary array for convenience
        power_dataset_tmp = np.array(self.eggs_resonance_calibration)

        # convert power to normalized amplitude values
        norm_ampl_dataset = np.zeros((2, len(power_dataset_tmp)))
        norm_ampl_dataset[:, 0] = power_dataset_tmp[:, 0]
        norm_ampl_dataset[:, 1] = 10 ** (power_dataset_tmp[:, 1] / 20)
        norm_ampl_dataset[:, 1] /= np.max(power_dataset_tmp[:, 1])
        self.set_dataset('calibration.eggs.resonance_ratio_curve_mhz', norm_ampl_dataset, persist=True, broadcast=True)

        # add data to data vault for visualization
        self.dv.add(np.array([power_dataset_tmp[:, 0], power_dataset_tmp[:, 1], norm_ampl_dataset[:, 1]]).transpose(), context=self.cr)
