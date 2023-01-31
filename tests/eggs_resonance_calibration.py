import numpy as np
from random import shuffle
from artiq.experiment import *

_DMA_HANDLE_EGGS_OFF =          "eggs_heating_eggs_off"


class EGGSResonanceCalibration(EnvExperiment):
    """
    EGGS Resonance Calibration

    Examines the resonance for the EGGS feedthrough and calculates the appropriate amplitude scalings.
    """

    # kernel_invariants = {}
    global_parameters = [
        "pmt_input_channel",
        "pmt_gating_edge",

        "time_profileswitch_delay_us",
        "time_repump_qubit_us",
        "time_doppler_cooling_us",
        "time_readout_us",
        "time_redist_us",

        "dds_board_num",
        "dds_board_qubit_num",
        "dds_probe_channel",
        "dds_pump_channel",
        "dds_repump_cooling_channel",
        "dds_repump_qubit_channel",
        "dds_qubit_channel",

        "freq_redist_mhz",
        "freq_pump_cooling_mhz",
        "freq_pump_readout_mhz",
        "freq_pump_rescue_mhz",
        "freq_repump_cooling_mhz",
        "freq_repump_qubit_mhz",

        "ampl_redist_pct",
        "ampl_pump_cooling_pct",
        "ampl_pump_readout_pct",
        "ampl_pump_rescue_pct",
        "ampl_repump_cooling_pct",
        "ampl_repump_qubit_pct",
        "ampl_qubit_pct"
    ]

    def build(self):
        """
        Set devices and arguments for the experiment.
        """
        self.setattr_device("core")
        self.setattr_device("core_dma")

        # experiment runs
        self.setattr_argument("repetitions",                            NumberValue(default=1, ndecimals=0, step=1, min=1, max=10000))

        # eggs heating
        self.setattr_argument("time_eggs_heating_ms",                   NumberValue(default=20, ndecimals=5, step=1, min=0.00001, max=10000))
        self.setattr_argument("freq_eggs_heating_secular_mhz",          NumberValue(default=1.6, ndecimals=5, step=0.1, min=0.001, max=1000000))
        self.setattr_argument("freq_eggs_heating_mhz_list",             Scannable(
                                                                            default=CenterScan(85, 5, 0.01, randomize=True),
                                                                            global_min=30, global_max=400, global_step=1,
                                                                            unit="MHz", scale=1, ndecimals=5
                                                                        ))

        # get global parameters
        for param_name in self.global_parameters:
            self.setattr_dataset(param_name, archive=True)


    def prepare(self):
        """
        Prepare things such that kernel functions have minimal overhead.
        """
        # ensure input has correct dimensions
        min_time_length = len(list(self.time_min_sideband_cooling_us_list))
        max_time_length = len(list(self.time_max_sideband_cooling_us_list))
        modes_length = len(list(self.freq_sideband_cooling_mhz_list))
        assert min_time_length == max_time_length == modes_length

        # PMT devices
        self.pmt_counter =                                              self.get_device("ttl{:d}_counter".format(self.pmt_input_channel))
        self.pmt_gating_edge =                                          getattr(self.pmt_counter, 'gate_{:s}_mu'.format(self.pmt_gating_edge))

        # eggs heating
        self.awg_board =                                                self.get_device("phaser0")
        self.awg_eggs =                                                 self.awg_board.channel[0]
        self.time_phaser_sample_mu =                                    np.int64(40)

        self.freq_eggs_heating_center_mhz =                             85
        self.freq_eggs_heating_mhz_list =                               np.array([(freq_mhz - self.freq_eggs_heating_center_mhz) for freq_mhz in self.freq_eggs_heating_mhz_list])
        self.time_eggs_heating_mu =                                     self.core.seconds_to_mu(self.time_eggs_heating_ms * ms)

        # ensure eggs heating time is a multiple of the phaser frame period
        # 4 ns/clock * 8 clock cycles * 10 words = 320ns
        if self.time_eggs_heating_mu % self.awg_board.t_frame:
            t_frame_multiples = round(self.time_eggs_heating_mu / self.awg_board.t_frame + 0.5)
            self.time_eggs_heating_mu = np.int64(self.awg_board.t_frame * t_frame_multiples)

        # calculate eggs sidebands amplitude, adjusted for eggs coupling resonance curve
        # todo: set correctly
        freq_eggs_sidebands_mhz_list =                                  np.zeros((len(self.freq_eggs_heating_mhz_list), 2))
        freq_eggs_sidebands_mhz_list[:, 0] =                            self.freq_eggs_heating_center_mhz + self.freq_eggs_heating_mhz_list - self.freq_eggs_heating_secular_mhz
        freq_eggs_sidebands_mhz_list[:, 1] =                            self.freq_eggs_heating_center_mhz + self.freq_eggs_heating_mhz_list + self.freq_eggs_heating_secular_mhz

        # create interpolated coupling resonance curve
        from scipy.interpolate import Akima1DInterpolator
        ampl_calib_points =                                             self.get_dataset('calibration.eggs.resonance_ratio_curve_mhz')
        ampl_calib_curve =                                              Akima1DInterpolator(ampl_calib_points[:, 0], ampl_calib_points[:, 1])

        # get scaled amplitudes
        self.ampl_eggs_heating_frac_list =                              np.zeros((len(freq_eggs_sidebands_mhz_list), 2), dtype=float)
        for i, freqs_mhz in enumerate(freq_eggs_sidebands_mhz_list):
            a1, a2 = ampl_calib_curve(freqs_mhz)
            # print('\t: {:f}/{:f} = {:f}'.format(a1, a2, a1+a2))
            self.ampl_eggs_heating_frac_list[i] = np.array([a2, a1]) / (a1+a2)

        print(self.ampl_eggs_heating_frac_list)
        print(freq_eggs_sidebands_mhz_list)


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
        handle_eggs_off =       self.core_dma.get_handle(_DMA_HANDLE_EGGS_OFF)
        self.core.break_realtime()

        # MAIN SEQUENCE
        for trial_num in range(self.repetitions):

            # sweep eggs rf frequencies
            for i in range(len(self.freq_eggs_heating_mhz_list)):

                # get eggs carrier frequency
                freq_eggs_mhz = self.freq_eggs_heating_mhz_list[i] * MHz

                # get eggs amplitude values
                ampl_eggs_rsb_frac = self.ampl_eggs_heating_frac_list[i][0]
                ampl_eggs_bsb_frac = self.ampl_eggs_heating_frac_list[i][1]

                # set eggs carrier via the DUC
                self.awg_eggs.set_duc_frequency(freq_eggs_mhz)
                self.awg_board.duc_stb()
                self.core.break_realtime()

                # sweep final pi-pulse frequency
                for freq_ftw in self.freq_qubit_scan_ftw:

                    # set readout frequency in advance
                    #self.dds_qubit.set_mu(freq_ftw, asf=self.ampl_qubit_asf, profile=0)
                    #self.core.break_realtime()

                    # initialize by running doppler cooling and spin polarization
                    #self.core_dma.playback_handle(handle_initialize)

                    # run sideband cooling cycles and repump afterwards
                    #self.core_dma.playback_handle(handle_sideband)

                    # enable eggs heating output
                    # self.awg_eggs.oscillator[0].set_amplitude_phase(amplitude=ampl_eggs_rsb_frac, clr=0)
                    self.awg_eggs.oscillator[0].set_amplitude_phase(amplitude=0.49, clr=0)
                    delay_mu(self.time_phaser_sample_mu)
                    self.awg_eggs.oscillator[1].set_amplitude_phase(amplitude=0.49, clr=0)
                    # self.awg_eggs.oscillator[1].set_amplitude_phase(amplitude=ampl_eggs_bsb_frac, clr=0)

                    # let eggs heating run, then turn off
                    at_mu(self.awg_board.get_next_frame_mu())
                    self.core_dma.playback_handle(handle_eggs_off)

                    # read out
                    self.core_dma.playback_handle(handle_readout)

                    # record data
                    #self.update_dataset(freq_ftw, self.pmt_counter.fetch_count(), freq_eggs_mhz + self.freq_eggs_heating_center_mhz)
                    #self.core.break_realtime()

        # disable rf output
        self.awg_eggs.en_trf_out(rf=0, lo=0)
        self.core.break_realtime()


    @kernel(flags={"fast-math"})
    def DMArecord(self):
        """
        Record onto core DMA the AOM sequence for a single data point.
        """
        # eggs sequence
        with self.core_dma.record(_DMA_HANDLE_EGGS_OFF):
            # apply eggs heating for given time;
            # no need to align here since we ensure time_eggs_heating_mu is an integer
            # multiple of phaser0.t_frame in prepare()
            delay_mu(self.time_eggs_heating_mu)

            # disable output
            self.awg_eggs.oscillator[0].set_amplitude_phase(amplitude=0., clr=1)
            delay_mu(self.time_phaser_sample_mu)
            self.awg_eggs.oscillator[1].set_amplitude_phase(amplitude=0., clr=1)
            delay_mu(self.time_phaser_sample_mu)
            self.awg_eggs.oscillator[0].set_amplitude_phase(amplitude=0., clr=0)
            delay_mu(self.time_phaser_sample_mu)
            self.awg_eggs.oscillator[1].set_amplitude_phase(amplitude=0., clr=0)


    @kernel(flags={"fast-math"})
    def prepareDevices(self):
        """
        Prepare devices for the experiment.
        """
        # initialize phaser
        self.awg_board.init(debug=True)
        self.core.break_realtime()

        # set nco to center eggs rf around the center frequency (should be 85 MHz) exactly
        self.awg_eggs.set_nco_frequency((-217.083495) * MHz)
        self.awg_eggs.set_nco_phase(0.)
        self.awg_board.dac_sync()
        self.core.break_realtime()

        # trf setup, and disable rf output while we set things up
        self.awg_eggs.set_att(0 * dB)
        self.awg_eggs.en_trf_out(rf=0, lo=0)
        self.core.break_realtime()

        # duc
        self.awg_eggs.set_duc_frequency(0 * MHz)
        self.awg_eggs.set_duc_cfg()
        self.awg_board.duc_stb()
        self.core.break_realtime()

        # oscillators (i.e. sidebands)
        self.awg_eggs.oscillator[0].set_frequency(-self.freq_eggs_heating_secular_mhz * MHz)
        self.awg_eggs.oscillator[0].set_amplitude_phase(0., clr=0)
        self.core.break_realtime()
        self.awg_eggs.oscillator[1].set_frequency(self.freq_eggs_heating_secular_mhz * MHz)
        self.awg_eggs.oscillator[1].set_amplitude_phase(0., clr=0)
        self.core.break_realtime()

        # re-enable rf output
        self.awg_eggs.en_trf_out(rf=1, lo=0)
        self.core.break_realtime()


    @rpc(flags={"async"})
    def update_dataset(self, freq_ftw, pmt_counts, time_mu):
        """
        Records values via rpc to minimize kernel overhead.
        """
        self.append_to_dataset('eggs_heating', [freq_ftw * self.ftw_to_mhz, pmt_counts, self.core.mu_to_seconds(time_mu)])

    @rpc
    def get_peaks(self, num_peaks):
        """
        Get the number of peaks
        """
        self.append_to_dataset('eggs_heating', [freq_ftw * self.ftw_to_mhz, pmt_counts, self.core.mu_to_seconds(time_mu)])


    def analyze(self):
        """
        Analyze the results from the experiment.
        """
        # turn dataset into numpy array for ease of use
        self.eggs_heating = np.array(self.eggs_heating)
