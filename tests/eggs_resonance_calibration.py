import labrad
import numpy as np
from os import environ
from time import sleep
from datetime import datetime
from artiq.experiment import *


class EGGSResonanceCalibration(EnvExperiment):
    """
    EGGS Resonance Calibration

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
        self.setattr_argument("repetitions",                                    NumberValue(default=8, ndecimals=0, step=1, min=1, max=10000))

        # eggs heating
        self.setattr_argument("freq_eggs_heating_mhz_list",                     Scannable(
                                                                                    default=RangeScan(70, 95, 2501, randomize=False),
                                                                                    global_min=30, global_max=400, global_step=1,
                                                                                    unit="MHz", scale=1, ndecimals=5
                                                                                ))
        # spectrum analyzer
        self.setattr_argument("spectrum_analyzer_bandwidth_khz",                NumberValue(default=1, ndecimals=5, step=1, min=0.00001, max=10000))
        self.setattr_argument("spectrum_analyzer_attenuation_internal_db",      NumberValue(default=10, ndecimals=5, step=1, min=0.00001, max=10000))
        self.setattr_argument("spectrum_analyzer_attenuation_external_db",      NumberValue(default=0, ndecimals=5, step=1, min=0.00001, max=10000))


    def prepare(self):
        """
        Prepare things such that kernel functions have minimal overhead.
        """
        # eggs heating devices
        self.awg_board =                                                self.get_device("phaser0")
        self.awg_eggs =                                                 self.awg_board.channel[0]

        # eggs heating timings
        self.time_phaser_sample_mu =                                    np.int64(40)

        # calculate eggs frequency
        self.freq_eggs_heating_mhz_list =                               list(self.freq_eggs_heating_mhz_list)
        freq_eggs_heating_center_mhz =                                  np.median(self.freq_eggs_heating_mhz_list)
        freq_eggs_heating_range_mhz =                                   np.max(self.freq_eggs_heating_mhz_list) - np.min(self.freq_eggs_heating_mhz_list)

        # connect to labrad
        self.cxn = labrad.connect(environ['LABRADHOST'], port=7682, tls_mode='off', username='', password='lab')
        self.sa = self.cxn.spectrum_analyzer_server
        self.dv = self.cxn.data_vault
        self.cr = self.cxn.context()

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
        print('\tmedian: {}'.format(freq_eggs_heating_center_mhz))
        # set up spectrum analyzer marker
        self.sa.peak_threshold(-90)
        self.sa.peak_excursion(15)
        self.sa.marker_toggle(1, True)
        # todo: move to actual labrad functions
        self.sa.gpib_write('DET:TRAC1 AVER')
        self.sa.gpib_write('DISP:ENAB 1')

        # set up data vault
        date = datetime.now()
        dataset_title_tmp = 'EGGS Resonance Calibration'
        trunk = '{0:s}_{1:02d}:{2:02d}'.format(dataset_title_tmp, date.hour, date.minute)
        trunk_tmp = ['', 'labrad', str(date.year), '{:02d}'.format(date.month), '{0:02d}'.format(date.day), trunk]
        self.dv.cd(trunk_tmp, True, context=self.cr)
        self.dv.new(
            dataset_title_tmp,
            [('Frequency', 'MHz')],
            [('Pickoff Power', 'Power', 'dBm')],
            context=self.cr
        )
        print("Data vault setup successful.")

        # set up artiq datasets
        self._iter_dataset = 0
        self.set_dataset("eggs_resonance_calibration", np.zeros([len(self.freq_eggs_heating_mhz_list), 2]))
        self.setattr_dataset("eggs_resonance_calibration")
        #self.set_dataset('management.completion_pct', (trial_num + 1) / self.repetitions * 100., broadcast=True, persist=True, archive=False)


    @kernel(flags={"fast-math"})
    def run(self):
        """
        Run the experimental sequence.
        """
        self.core.reset()

        # prepare phaser
        self.phaser_prepare()
        delay_mu(500000000)
        self.core.wait_until_mu(now_mu())
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
        self.awg_eggs.set_att(20 * dB)
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

        # average multiple measurements
        peak_power_dbm_list = np.zeros(self.repetitions)
        for i in range(self.repetitions):
            sleep(1.25 * sweep_time_s)
            peak_power_dbm_list[i] = self.sa.marker_amplitude(1)
        peak_power_dbm_avg = np.mean(peak_power_dbm_list)

        # zoom in on amplitude
        # todo

        # save to dataset (artiq)
        self.mutate_dataset('eggs_resonance_calibration', self._iter_dataset, np.array([peak_freq_mhz, peak_power_dbm_avg]))
        self._iter_dataset += 1
        # save to dataset (labrad)
        self.dv.add(peak_freq_mhz, peak_power_dbm_avg, context=self.cr)


    def analyze(self):
        """
        Analyze the results from the experiment.
        """
        # print results
        print(self.eggs_resonance_calibration)
        # reenable spec anal display
        self.sa.gpib_write('DISP:ENAB 1')
        # turn dataset into numpy array for ease of use
        #self.eggs_heating = np.array(self.eggs_heating)
