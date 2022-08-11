import numpy as np
from artiq.experiment import *
# todo: upload data to labrad
# todo: check scannable works correctly
_DMA_HANDLE_READOUT = "time_sweep_readout"


class TimeSweep(EnvExperiment):
    """
    729nm Laser Time Sweep
    Measures ion fluorescence vs AOM detuning
    """

    #kernel_invariants = {}

    def build(self):
        print('kk123')
        """
        Set devices and arguments for the experiment.
        """
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # experiment runs
        self.setattr_argument("repetitions",            NumberValue(default=150, ndecimals=0, step=1, min=1, max=10000))

        # timing
        self.setattr_argument("time_854_us",            NumberValue(default=400, ndecimals=5, step=1, min=1, max=10000))
        self.setattr_argument("time_readout_us",        NumberValue(default=500, ndecimals=5, step=1, min=1, max=10000))

        # AOM DDS channels
        self.setattr_argument("dds_board_num",          NumberValue(default=0, ndecimals=0, step=1, min=0, max=1))
        self.setattr_argument("dds_qubit_channel",      NumberValue(default=0, ndecimals=0, step=1, min=0, max=3))

        # qubit frequency scan  #RangeScan(100, 120, 5)
       # self.setattr_argument("freq_qubit_scan_mhz",    Scannable(default=RangeScan(109, 111, 401),
         #                                                         global_min=60, global_max=200, global_step=1,
           #                                                        unit="MHz", scale=1, ndecimals=3))
        self.freq_qubit_scan_mhz = [110]
        # qubit time scan
        self.setattr_argument("time_sweep_us",          Scannable(default=RangeScan(1, 2000, 200, randomize=True),
                                                                  global_min=1, global_max=10000, global_step=1,
                                                                  unit="us", scale=1, ndecimals=0))

        # PMT
        self.setattr_argument("pmt_input_channel",      NumberValue(default=0, ndecimals=0, step=1, min=0, max=3))
        self.setattr_argument("pmt_gating_edge",        EnumerationValue(["rising", "falling", "both"], default="rising"))

    def prepare(self):
        """
        Set up the dataset and prepare things such that
        the kernel functions have minimal overhead.
        """
        # PMT devices
        self.pmt_counter = self.get_device("ttl_counter{:d}".format(self.pmt_input_channel))
        self.pmt_gating_edge = getattr(self.pmt_counter, 'gate_{:s}_mu'.format(self.pmt_gating_edge))

        # convert time values to machine units
        self.time_854_mu = self.core.seconds_to_mu(self.time_854_us * us)
        self.time_readout_mu = self.core.seconds_to_mu(self.time_readout_us * us)

        # DDS devices
        self.dds_board = self.get_device("urukul{:d}_cpld".format(self.dds_board_num))
        #self.dds_qubit = self.get_device("urukul{:d}_ch{:d}".format(self.dds_board_num, self.dds_qubit_channel))
        self.dds_qubit = self.get_device("urukul{:d}_ch{:d}".format(0, 0))

        # convert dds values to machine units - qubit
        self.ftw_to_frequency = 1e9 / (2**32 - 1)
        self.freq_qubit_scan_mhz2 = list(self.freq_qubit_scan_mhz)

        # convert time values for sweep to machine values
        self.time_sweep_mu2 = [self.core.seconds_to_mu(val * us) for val in self.time_sweep_us]

        # convert dds values to machine units - everything else
        self.ampl_qubit_asf = self.dds_qubit.amplitude_to_asf(0.5)

        # set up datasets
        self.set_dataset("laser_scan", [], broadcast=True)
        #self.set_dataset("tmp", [], broadcast=True)

        # tmp remove:
        self.setattr_device('urukul1_cpld')

    @kernel(flags={"fast-math"})
    def run(self):
        """
        Run the experimental sequence.
        """
        self.core.reset()

        # prepare devices
        self.prepareDevices()
        self.core.break_realtime()

        # record dma and get handle
        self.DMArecord()
        handle = self.core_dma.get_handle(_DMA_HANDLE_READOUT)
        self.core.break_realtime()

        # MAIN SEQUENCE
        for trial_num in range(self.repetitions):

            # set frequencies
            for freq_mhz in self.freq_qubit_scan_mhz2:
                self.core.break_realtime()

                # set freq and ampl for qubit
                freq_mu = self.dds_qubit.frequency_to_ftw(freq_mhz * MHz)
                self.dds_qubit.set_mu(freq_mu, asf=self.ampl_qubit_asf)
                self.core.break_realtime()

                # sweep time
                for time_mu in self.time_sweep_mu2:

                    # turn on 854 to pump back down
                    with parallel:
                        self.urukul1_cpld.cfg_switches(0b1100)
                        delay_mu(self.time_854_mu)
                    self.urukul1_cpld.cfg_switches(0b0100)

                    # get pmt counts w/397 onto calibrate
                    self.core_dma.playback_handle(handle)
                    pmt_calib = self.pmt_counter.fetch_count()
                    self.core.break_realtime()

                    # rabi flopping w/qubit laser
                    with parallel:
                        self.dds_qubit.cfg_sw(1)
                        delay_mu(time_mu)
                    self.dds_qubit.cfg_sw(0)

                    # get pmt counts (actual)
                    self.core_dma.playback_handle(handle)
                    pmt_actual = self.pmt_counter.fetch_count()

                    # update dataset
                    with parallel:
                        #self.update_dataset(time_mu, freq_mhz, self.pmt_counter.fetch_count(), pmt_calib)
                        self.update_dataset(time_mu, freq_mhz, pmt_actual, pmt_calib)
                        self.core.break_realtime()

        # reset after experiment
        self.urukul1_cpld.cfg_switches(0b1110)
        self.dds_qubit.cfg_sw(0)

    @kernel(flags={"fast-math"})
    def DMArecord(self):
        """
        Record onto core DMA the AOM sequence for a single data point.
        """
        with self.core_dma.record(_DMA_HANDLE_READOUT):
            self.urukul1_cpld.cfg_switches(0b0110)
            self.pmt_gating_edge(self.time_readout_mu)
            self.urukul1_cpld.cfg_switches(0b0100)

    @kernel(flags={"fast-math"})
    def prepareDevices(self):
        """
        Prepare devices for the experiment.
        """
        self.core.break_realtime()

        # initialize dds board
        self.core.break_realtime()

        # initialize qubit AOM and set waveform
        #self.dds_qubit.init()
        # tmp remove: dds set_att cfg sw stuff
        #self.dds_qubit.set_att(2 * dB)
        self.core.break_realtime()
        self.dds_qubit.cfg_sw(1)

        # turn all core lasers on
        self.urukul1_cpld.cfg_switches(0b1110)

    @rpc(flags={"async"})
    def update_dataset(self, time_mu, freq_mhz, pmt_counts, pmt_calib):
        """
        Records values via rpc to minimize kernel overhead.
        """
        self.append_to_dataset('laser_scan', [self.core.mu_to_seconds(time_mu), freq_mhz, pmt_counts, pmt_calib])

    def analyze(self):
        """
        Analyze the results from the experiment.
        """
        pass
        # data = self.get_dataset("laser_scan")
        # pmt_data = data[2500:, 1]
        #
        # print("results:")
        # print("\tpmt counts: {:.2f} +/- {:.2f}".format(np.mean(pmt_data), np.std(pmt_data)))
        # print("\tphotodiode value: {:.2f} +/- {:.2f}".format(np.mean(pd_data), np.std(pd_data)))
