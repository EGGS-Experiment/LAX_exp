import labrad
import numpy as np
from os import environ
from time import sleep
from artiq.experiment import *

_DMA_HANDLE_PARAMETRIC_SEQUENCE = "parametric_heating_sequence"


class ParametricHeating(EnvExperiment):
    """
    Parametric Heating Scan

    Modulates the trap RF to try and detect the secular frequency.
    """

    # kernel_invariants = {}
    global_parameters = [
        "pmt_input_channel",
        "pmt_gating_edge",

        "dds_board_num",
        "dds_pump_channel",
        "dds_repump_cooling_channel",

        "freq_pump_cooling_mhz",
        "freq_pump_readout_mhz",
        "freq_pump_rescue_mhz",
        "freq_repump_cooling_mhz",

        "ampl_pump_cooling_pct",
        "ampl_pump_readout_pct",
        "ampl_pump_rescue_pct",
        "ampl_repump_cooling_pct",

        "time_profileswitch_delay_us"
    ]

    def build(self):
        """
        Set devices and arguments for the experiment.
        """
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # experiment runs
        self.setattr_argument("repetitions",                                NumberValue(default=100, ndecimals=0, step=1, min=1, max=10000))

        # parametric heating values
        self.setattr_argument("time_parametric_us",                         NumberValue(default=100, ndecimals=3, step=10, min=1, max=1000000))
        self.setattr_argument("time_recovery_us",                           NumberValue(default=100, ndecimals=3, step=10, min=1, max=1000000))
        self.setattr_argument("ampl_parametric_vpp",                        NumberValue(default=0.05, ndecimals=5, step=0.1, min=0, max=10000))
        self.setattr_argument("freq_parametric_mhz_list",                   Scannable(
                                                                                default=RangeScan(0.9, 1.2, 31, randomize=True),
                                                                                global_min=0, global_max=1000, global_step=1,
                                                                                unit="MHz", scale=1, ndecimals=6
                                                                            ))

        # get global parameters
        for param_name in self.global_parameters:
            self.setattr_dataset(param_name, archive=True)


    def prepare(self):
        """
        Set up the dataset and prepare things such that
        the kernel functions have minimal overhead.
        """
        # PMT devices
        self.pmt_counter =                                                  self.get_device("ttl{:d}_counter".format(self.pmt_input_channel))
        self.pmt_gating_edge =                                              getattr(self.pmt_counter, 'gate_{:s}_mu'.format(self.pmt_gating_edge))

        # function generator
        self.ttl_function_generator =                                       self.get_device("ttl8")

        # DDS devices
        self.dds_board =                                                    self.get_device("urukul{:d}_cpld".format(self.dds_board_num))
        self.dds_pump =                                                     self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_pump_channel))
        self.dds_repump_cooling =                                           self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_repump_cooling_channel))

        # get DDS waveforms
        self.freq_pump_cooling_ftw =                                        self.dds_qubit.frequency_to_ftw(self.freq_pump_cooling_mhz * MHz)
        self.freq_pump_readout_ftw =                                        self.dds_qubit.frequency_to_ftw(self.freq_pump_readout_mhz * MHz)
        self.freq_pump_rescue_ftw =                                         self.dds_qubit.frequency_to_ftw(self.freq_pump_rescue_mhz * MHz)
        self.freq_repump_cooling_ftw =                                      self.dds_pump.frequency_to_ftw(self.freq_repump_cooling_mhz * MHz)

        # convert amplitude to asf
        self.ampl_pump_cooling_asf =                                        self.dds_qubit.amplitude_to_asf(self.ampl_pump_cooling_pct / 100)
        self.ampl_pump_readout_asf =                                        self.dds_qubit.amplitude_to_asf(self.ampl_pump_readout_pct / 100)
        self.ampl_pump_rescue_asf =                                         self.dds_qubit.amplitude_to_asf(self.ampl_pump_rescue_pct / 100)
        self.ampl_repump_cooling_asf =                                      self.dds_pump.amplitude_to_asf(self.ampl_repump_cooling_pct / 100)

        # timings
        self.time_profileswitch_delay_mu =                                  self.core.seconds_to_mu(self.time_profileswitch_delay_us * us)
        self.time_fuction_generator_delay_mu =                              self.core.seconds_to_mu(300 * ns)
        self.time_parametric_mu =                                           self.core.seconds_to_mu(self.time_parametric_us * us)
        self.time_recovery_mu =                                             self.core.seconds_to_mu(self.time_recovery_us * us)

        # parametric heating frequencies
        self.freq_parametric_hz_list =                                      2 * MHz * np.array(list(self.freq_tickle_mhz))


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
            if 'DG2P' in dev_desc:
                dev_exists = True
                self.fg.select_device(dev_num)

        # raise error if function generator doesn't exist
        if not dev_exists:
            raise Exception("Error: modulation function generator not detected.")

        # set up function generator
        self.fg.gpib_write(':OUTP:IMP 50')
        self.fg.toggle(0)
        self.fg.amplitude(self.ampl_parametric_vpp)
        self.fg.burst(True)
        self.fg.burst_mode('GAT')
        self.fg.toggle(1)
        # todo: do we have to set trigger stuff?

        # set up datasets
        self.set_dataset("parametric_heating", np.zeros((self.repetitions * len(self.freq_parametric_hz_list), 2)))
        self.setattr_dataset("parametric_heating")


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
        handle_sequence = self.core_dma.get_handle(_DMA_HANDLE_PARAMETRIC_SEQUENCE)
        self.core.break_realtime()

        # MAIN SEQUENCE
        for freq_mhz in self.freq_parametric_hz_list:

            # set frequency
            self.frequency_set(freq_mhz)
            self.core.break_realtime()

            # repeat experiment
            for trial_num in range(self.repetitions):

                # run sequence
                self.core_dma.playback_handle(handle_sequence)

                # add data to dataset
                with parallel:
                    delay_mu(self.time_recovery_mu)
                    self.update_dataset(freq_mhz, self.pmt_counter.fetch_count())

        # reset devices
        self.dds_board.set_profile(0)
        self.dds_board.cpld.io_update.pulse_mu(8)
        self.core.break_realtime()

        self.dds_board.cfg_switches(0b1110)
        self.core.break_realtime()
        self.ttl_function_generator.off()


    @kernel(flags={"fast-math"})
    def DMArecord(self):
        """
        Record onto core DMA the AOM sequence for a single data point.
        """
        # parametric heating sequence
        with self.core_dma.record(_DMA_HANDLE_PARAMETRIC_SEQUENCE):

            # set readout profile
            self.dds_board.set_profile(1)
            delay_mu(self.time_profileswitch_delay_mu)

            # start parametric heating
            self.ttl_function_generator.on()
            delay_mu(self.time_fuction_generator_delay_mu)

            # read pmt
            self.pmt_gating_edge(self.time_parametric_mu)

            # stop parametric heating, start rescuing ion
            with parallel:
                self.ttl_function_generator.off()
                self.dds_board.set_profile(2)
            delay_mu(self.time_profileswitch_delay_mu)


    @kernel(flags={"fast-math"})
    def prepareDevices(self):
        """
        Prepare devices for the experiment.
        """
        # set up AOMs
        self.dds_pump.set_mu(self.freq_pump_cooling_ftw, asf=self.ampl_pump_cooling_asf, profile=0)
        self.dds_pump.set_mu(self.freq_pump_readout_ftw, asf=self.ampl_pump_readout_asf, profile=1)
        self.dds_pump.set_mu(self.freq_pump_rescue_ftw, asf=self.ampl_pump_rescue_asf, profile=2)
        self.core.break_realtime()

        self.dds_repump_cooling.set_mu(self.freq_repump_cooling_ftw, asf=self.ampl_repump_cooling_asf, profile=0)
        self.dds_repump_cooling.set_mu(self.freq_repump_cooling_ftw, asf=self.ampl_repump_cooling_asf, profile=1)
        self.dds_repump_cooling.set_mu(self.freq_repump_cooling_ftw, asf=self.ampl_repump_cooling_asf, profile=2)
        self.core.break_realtime()

        # set up ttl for function generator trigger
        self.ttl_function_generator.off()

    @rpc
    def frequency_set(self, freq_hz):
        """
        Set the RF to the desired frequency.
        """
        freq_set_hz = self.fg.frequency(freq_hz)
        print('\tfrequency set: {}'.format(freq_set_hz))

    def fg_write(self, msg):
        """
        Write a GPIB message to the function generator.
        """
        self.fg.gpib_write(msg)

    @rpc(flags={"async"})
    def update_dataset(self, freq_mhz, pmt_counts):
        """
        Records values via rpc to minimize kernel overhead.
        """
        self.parametric_heating[self._result_iter] = np.array([freq_mhz, pmt_counts])

    def analyze(self):
        """
        Analyze the results from the experiment.
        """
        self.fg.toggle(0)
        #self.ttl_function_generator.on()
        #self.micromotion_compensation[:, 0] = float(self.micromotion_compensation[:, 0] / 2**16)
