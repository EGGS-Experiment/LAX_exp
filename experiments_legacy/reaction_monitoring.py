import labrad
import numpy as np
from os import environ
from time import sleep

from artiq.experiment import *

_DMA_HANDLE_SEQUENCE = "reaction_monitoring_sequence"


class ReactionMonitoring(EnvExperiment):
    """
    Reaction Monitoring

    Monitor the HCl-H2Cl reaction dependence on intensity.
    """

    # kernel_invariants = {}
    global_parameters = [
        "pmt_input_channel",
        "pmt_gating_edge",

        "dds_board_num",
        "dds_pump_channel",
        "dds_repump_cooling_channel",

        "freq_repump_cooling_mhz",
        "ampl_repump_cooling_pct",
    ]

    def build(self):
        """
        Set devices and arguments for the experiment.
        """
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # experiment runs
        self.setattr_argument("repetitions",                            NumberValue(default=1, ndecimals=0, step=1, min=1, max=10000))

        # cooling parameters
        self.setattr_argument("freq_pump_cooling_mhz",                  NumberValue(default=110, ndecimals=6, step=1, min=1, max=10000))
        self.setattr_argument("ampl_pump_cooling_pct_list",             Scannable(
                                                                            default=RangeScan(30.0, 40.0, 2, randomize=True),
                                                                            global_min=0, global_max=1000, global_step=1,
                                                                            unit="%", scale=1, ndecimals=2
                                                                        ))


        # tmp remove
        self.setattr_argument("amplitude_pct_time_s_dict",              PYONValue({30.0: 60, 40.0: 100, 2:100}), tooltip='gthkim')

        # timing
        self.setattr_argument("time_cooling_s",                         NumberValue(default=120, ndecimals=5, step=1, min=1, max=10000))
        self.setattr_argument("time_control_s",                         NumberValue(default=30, ndecimals=5, step=1, min=1, max=10000))

        # fluorescence
        self.setattr_argument("sample_time_pmt_us",                     NumberValue(default=3000, ndecimals=3, step=1, min=1, max=1000000))
        self.setattr_argument("sample_rate_pmt_hz",                     NumberValue(default=100, ndecimals=3, step=1, min=1, max=10000))


        # tickling parameters
        self.setattr_argument("tickle_freq_khz",                        NumberValue(default=602, ndecimals=4, step=1, min=1, max=10000))
        self.setattr_argument("tickle_ampl_v",                          NumberValue(default=0.2, ndecimals=4, step=0.1, min=0.001, max=5))

        # get global parameters
        for param_name in self.global_parameters:
            self.setattr_dataset(param_name, archive=True)


    def prepare(self):
        """
        Set up the dataset and prepare things such that
        the kernel functions have minimal overhead.
        """
        # PMT devices
        self.pmt_counter =                                              self.get_device("ttl{:d}_counter".format(self.pmt_input_channel))
        self.pmt_gating_edge =                                          getattr(self.pmt_counter, 'gate_{:s}_mu'.format(self.pmt_gating_edge))

        # convert PMT values
        self.sample_time_pmt_mu =                                       self.core.seconds_to_mu(self.sample_time_pmt_us * us)
        delay_time_pmt_s_tmp =                                          (1 / self.sample_rate_pmt_hz) - (self.sample_time_pmt_us * us)
        assert delay_time_pmt_s_tmp > 0, "Error: Invalid PMT count timings."
        self.delay_time_pmt_mu =                                        self.core.seconds_to_mu(delay_time_pmt_s_tmp)

        # set up pmt sampling iterators
        self.iter_pmt_test =                                            np.arange(self.time_cooling_s * self.sample_rate_pmt_hz)
        self.iter_pmt_control =                                         np.arange(self.time_control_s * self.sample_rate_pmt_hz)

        # tmp remove
        # self.iter_pmt_test_list_tmp =                                       [np.arange(time_s * self.sample_rate_pmt_hz) for time_s in list(self.amplitude_pct_time_s_dict.values())]
        # self.ampl_pump_cooling_asf_list_tmp =                               np.array([self.dds_pump.amplitude_to_asf(ampl_pct / 100)
        #                                                                           for ampl_pct in list(self.amplitude_pct_time_s_dict.keys())])
        # num_points =                                                    np.sum(np.array([len(iter_arr) for iter_arr in self.iter_pmt_test_list_tmp]))

        # DDS devices
        self.dds_board =                                                self.get_device("urukul{:d}_cpld".format(self.dds_board_num))
        self.dds_pump =                                                 self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_pump_channel))
        self.dds_repump_cooling =                                       self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_repump_cooling_channel))
        self.asf_to_pct =                                               1e2 / (2 ** 14 - 1)

        # get repump values
        self.freq_repump_cooling_ftw =                                  self.dds_pump.frequency_to_ftw(self.freq_repump_cooling_mhz * MHz)
        self.ampl_repump_cooling_asf =                                  self.dds_pump.amplitude_to_asf(self.ampl_repump_cooling_pct / 100)

        # convert cooling values
        self.freq_pump_cooling_ftw =                                    self.dds_pump.frequency_to_ftw(self.freq_pump_cooling_mhz * MHz)
        self.ampl_pump_cooling_asf_list =                               np.array([self.dds_pump.amplitude_to_asf(ampl_pct / 100) for ampl_pct in list(self.ampl_pump_cooling_pct_list)])

        # set up datasets
        self._iter_dataset_test =                                       0
        self.set_dataset("reaction_monitoring",                         np.zeros((
                                                                            len(self.iter_pmt_test) * len(self.ampl_pump_cooling_asf_list) * self.repetitions,
                                                                            3
                                                                        )))
        # self.set_dataset("reaction_monitoring",                         np.zeros((
        #                                                                     len(self.iter_pmt_test) * len(self.ampl_pump_cooling_asf_list) * self.repetitions,
        #                                                                     3
        #                                                                 )))
        self.setattr_dataset("reaction_monitoring")

        self._iter_dataset_control =                                    0
        self.set_dataset("reaction_monitoring_control",                 np.zeros((
                                                                            len(self.iter_pmt_control) * len(self.ampl_pump_cooling_asf_list) * self.repetitions * 2,
                                                                            4
                                                                        )))
        self.setattr_dataset("reaction_monitoring_control")


        # connect to labrad
        self.cxn = labrad.connect(environ['LABRADHOST'], port=7682, tls_mode='off', username='', password='lab')
        self.fg = self.cxn.function_generator_server

        # set up function generator
        # get list of function generators
        fg_dev_list = self.fg.list_devices()
        fg_dev_dict = dict(tuple(fg_dev_list))

        # select correct function generator
        dev_exists = False
        for dev_num, dev_desc in fg_dev_dict.items():
            if 'MY4800' in dev_desc:
                dev_exists = True
                self.fg.select_device(dev_num)

        # raise error if function generator doesn't exist
        if not dev_exists:
            raise Exception("Error: modulation function generator not detected.")

        # set up function generator
        # set corect mode
        self.fg.burst_mode('TRIG')
        self.fg.burst(False)
        # set waveform
        self.fg.toggle(0)
        self.fg.frequency(self.tickle_freq_khz * kHz)
        self.fg.amplitude(self.tickle_ampl_v)
        self.fg.toggle(1)


    @kernel(flags={"fast-math"})
    def run(self):
        """
        Run the experimental sequence.
        """
        self.core.reset()

        # prepare devices
        self.prepareDevices()

        # record dma and get handle
        self.DMArecord()
        handle_sequence = self.core_dma.get_handle(_DMA_HANDLE_SEQUENCE)
        self.core.break_realtime()

        # MAIN SEQUENCE
        for trial_num in range(self.repetitions):
            self.core.break_realtime()

            # sweep cooling amplitudes
            for ampl_asf in self.ampl_pump_cooling_asf_list:

                # set cooling beam amplitude
                self.dds_pump.set_mu(self.freq_pump_cooling_ftw, asf=ampl_asf)
                self.core.break_realtime()

                # get control data - before
                self.tickle_toggle(False)
                self.core.break_realtime()
                for i in self.iter_pmt_control:
                    # run pulse sequence from core DMA
                    self.core_dma.playback_handle(handle_sequence)

                    # record pmt counts to dataset
                    with parallel:
                        self.update_dataset_control(ampl_asf, i, self.pmt_counter.fetch_count(), 0)
                        delay_mu(self.delay_time_pmt_mu)

                # get treatment data w/tickle on
                self.tickle_toggle(True)
                self.core.break_realtime()
                for i in self.iter_pmt_test:
                    # run pulse sequence from core DMA
                    self.core_dma.playback_handle(handle_sequence)

                    # record pmt counts to dataset
                    with parallel:
                        self.update_dataset_test(ampl_asf, i, self.pmt_counter.fetch_count())
                        delay_mu(self.delay_time_pmt_mu)

                # get control data - after
                self.tickle_toggle(False)
                self.core.break_realtime()
                for i in self.iter_pmt_control:
                    # run pulse sequence from core DMA
                    self.core_dma.playback_handle(handle_sequence)

                    # record pmt counts to dataset
                    with parallel:
                        self.update_dataset_control(ampl_asf, i, self.pmt_counter.fetch_count(), 1)
                        delay_mu(self.delay_time_pmt_mu)

        # reset board profiles
        self.core.break_realtime()
        self.dds_board.set_profile(0)

        # reset AOMs after experiment
        self.dds_board.cfg_switches(0b1110)


    @kernel(flags={"fast-math"})
    def prepareDevices(self):
        """
        Prepare devices for the experiment.
        """
        self.core.break_realtime()

        # set cooling repump waveform
        self.dds_repump_cooling.set_mu(self.freq_repump_cooling_ftw, asf=self.ampl_repump_cooling_asf)
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def DMArecord(self):
        """
        Record onto core DMA the AOM sequence for a single data point.
        """
        # reset sequence
        with self.core_dma.record(_DMA_HANDLE_SEQUENCE):
            self.pmt_gating_edge(self.sample_time_pmt_mu)

    def tickle_toggle(self, status):
        """
        Toggle the output state of the tickle.
        """
        self.fg.toggle(status)
        #sleep(0.1)
        status_set = self.fg.toggle()
        # tmp remove
        print('\tstatus set: {}'.format(status_set))

    @rpc(flags={"async"})
    def update_dataset_test(self, ampl_asf, pmt_iter, pmt_counts):
        """
        Records values via rpc to minimize kernel overhead.
        """
        data_tmp = np.array([ampl_asf * self.asf_to_pct, pmt_iter / self.sample_rate_pmt_hz, pmt_counts])
        self.mutate_dataset('reaction_monitoring', self._iter_dataset_test, data_tmp)
        self._iter_dataset_test += 1

    @rpc(flags={"async"})
    def update_dataset_control(self, ampl_asf, pmt_iter, pmt_counts, before_or_after):
        """
        Records values via rpc to minimize kernel overhead.
        """
        data_tmp = np.array([ampl_asf * self.asf_to_pct, pmt_iter / self.sample_rate_pmt_hz, pmt_counts, before_or_after])
        self.mutate_dataset('reaction_monitoring_control', self._iter_dataset_control, data_tmp)
        self._iter_dataset_control += 1


    def analyze(self):
        """
        Analyze the results from the experiment.
        """
        pass
