from artiq.experiment import *
from numpy import array, zeros, int32, int64, nan
from artiq.coredevice.ad9910 import PHASE_MODE_CONTINUOUS

from LAX_exp.language import *
from LAX_exp.system.subsequences import (
    InitializeQubit, QubitPulseShape, SidebandCoolContinuousRAM, NoOperation, ReadoutAdaptive, Readout
)
from LAX_exp.system.objects.MeanFilter import MeanFilter
# todo: finish implementing changes


class QubitAlignmentRDX(LAXExperiment, Experiment):
    """
    Utility: Qubit Alignment RDX

    Excite a dark-state qubit transition for a short, fixed time
        to examine alignment of the qubit beam.
    Note: uses Adaptive Readout to speed up readout procedure, which requires
        knowledge of bright/dark count rates. These are "count_rate_bright_3ms" and "count_rate_dark_3ms,"
        and are set in the dataset manager under "sequences.adaptive_readout."
    """
    name = 'Qubit Alignment RDX'
    kernel_invariants = {
        # configs, conversions, etc
        "repetitions", "_filter_arr", "time_readout_s",

        # hardware parameters etc.
        "freq_qubit_ftw", "att_qubit_mu", "ampl_qubit_asf", "time_qubit_mu_list",
        "profile_729_readout", "profile_729_SBC",

        # subsequences & objects
        "initialize_subsequence", "readout_subsequence", "pulseshape_subsequence",
        "sbc_subsequence", "doppler_subsequence",
    }

    def build_experiment(self):
        """
        Set devices and arguments for the experiment.
        """
        # general
        self.setattr_argument("cooling_type",       EnumerationValue(["Doppler", "SBC"], default="Doppler"))
        self.setattr_argument('time_total_s',       NumberValue(default=200, precision=0, step=100, min=5, max=100000, scale=1., unit="s"), group='timing')
        self.setattr_argument('samples_per_point',  NumberValue(default=50, precision=0, step=10, min=15, max=500), group='timing')

        # qubit
        self.setattr_argument("time_qubit_us_list",   Scannable(
                                                        default=[
                                                            ExplicitScan([1.5, 15]),
                                                        ],
                                                        global_min=0.001, global_max=100000, global_step=1,
                                                        unit="us", scale=1, precision=5
                                                    ), group="qubit")
        self.setattr_argument("freq_qubit_mhz", NumberValue(default=101.0598, precision=5, step=1, min=1, max=10000, scale=1., unit="MHz"), group='qubit')
        self.setattr_argument("ampl_qubit_pct", NumberValue(default=50., precision=3, step=5., min=0.01, max=50., scale=1., unit="%"), group='qubit')
        self.setattr_argument("att_qubit_db",   NumberValue(default=8, precision=1, step=0.5, min=8, max=31.5, scale=1., unit="dB"), group='qubit')
        self.setattr_argument("enable_pulseshaping",    BooleanValue(default=False), group='qubit')

        # allocate profiles on 729nm for different subsequences
        self.profile_729_readout = 0
        self.profile_729_SBC =     1

        # instantiate subsequences
        self.sbc_subsequence =  SidebandCoolContinuousRAM(
            self, profile_729=self.profile_729_SBC, profile_854=3,
            ram_addr_start_729=0, ram_addr_start_854=0, num_samples=200
        )
        self.pulseshape_subsequence =   QubitPulseShape(
            self, ram_profile=self.profile_729_readout, ram_addr_start=500, num_samples=100,
            ampl_max_pct=self.ampl_qubit_pct,
        )
        self.initialize_subsequence =   InitializeQubit(self)
        # self.readout_subsequence =      ReadoutAdaptive(self, time_bin_us=10, error_threshold=1e-2)
        self.readout_subsequence = Readout(self)
        self.doppler_subsequence =      NoOperation(self)

        # relevant devices
        self.setattr_device('qubit')
        self.setattr_device('pmt')

    def prepare_experiment(self):
        """
        Prepare & precompute experimental values.
        """
        # choose correct cooling subsequence
        if self.cooling_type == "Doppler":  self.cooling_subsequence = self.doppler_subsequence
        elif self.cooling_type == "SBC":    self.cooling_subsequence = self.sbc_subsequence

        # convert qubit parameters
        self.freq_qubit_ftw =   self.qubit.frequency_to_ftw(self.freq_qubit_mhz * MHz)
        self.ampl_qubit_asf =   self.qubit.amplitude_to_asf(self.ampl_qubit_pct / 100.)
        self.att_qubit_mu =     att_to_mu(self.att_qubit_db * dB)
        self.time_qubit_mu_list =   [self.core.seconds_to_mu(time_us * us)
                                     for time_us in self.time_qubit_us_list]

        # calculate configuration and timings
        time_per_point_us = 3000    # guess time per shot lol
        self.repetitions =  round(self.time_total_s / (self.samples_per_point * time_per_point_us * us))
        self.num_times =    len(self.time_qubit_mu_list)
        self._state_array =         zeros(self.samples_per_point, dtype=int32) # todo: document
        self._time_exp_shot_s =     0.   # store time taken per shot
        # get readout time for later conversion
        self.time_readout_s =   self.get_parameter('time_readout_us', group='timing', override=False) * us

        # create filter objects & prepare them (b/c only LAXExperiment classes call their own children)
        self._filter_arr = [MeanFilter(self, filter_length=self.samples_per_point) for _ in self.time_qubit_us_list]
        for filter in self._filter_arr: filter.prepare()

    @rpc
    def initialize_plotting(self) -> TNone:
        """
        Configure datasets and applets for plotting.
        """
        # create datasets for storing counts
        # workaround: set first element to 0 to avoid "RuntimeWarning: All-NaN slice encountered"
        state_x_arr = zeros(self.repetitions) * nan
        state_x_arr[0] = 0
        state_y_arr = zeros(self.repetitions) * nan
        state_y_arr[0] = 0

        # prepare datasets for storing counts
        self.set_dataset('temp.qubit_align.counts_x', state_x_arr, broadcast=True, persist=True, archive=False)
        for i in range(self.num_times):
            self.set_dataset('temp.qubit_align.counts_y_{:d}'.format(i), state_y_arr,
                             broadcast=True, persist=True, archive=False)
            # initialize plotting applet
            stra = '${artiq_applet}' + (
                'plot_xy temp.qubit_align.counts_y_{:d}'
                ' --x temp.qubit_align.counts_x'
                ' --title "Qubit Alignment - {:f}"'
            ).format(i, self.time_qubit_mu_list[i])
            self.ccb.issue(
                # name of broadcast & applet name
                "create_applet", "qubit_alignment_{:d}".format(i),
                stra, # command
                group=["alignment"]  # folder directory for applet
            )

    @property
    def results_shape(self):
        return (self.repetitions, 3)


    '''
    MAIN SEQUENCE
    '''
    @kernel(flags={"fast-math"})
    def initialize_experiment(self) -> TNone:
        # prepare plotting
        # note: do it here instead of prepare to prevent overriding other experiments
        self.initialize_plotting()
        self.core.break_realtime()

        # prepare qubit beam for readout
        if self.enable_pulseshaping:
            self.pulseshape_subsequence.configure(self.time_qubit_mu_list[0])
            self.qubit.set_ftw(self.freq_qubit_ftw)
        else:
            self.qubit.set_mu(self.freq_qubit_ftw, asf=self.ampl_qubit_asf,
                              profile=self.profile_729_readout,
                              phase_mode=PHASE_MODE_CONTINUOUS)
        delay_mu(25000)

        # record alignment sequence
        with self.core_dma.record('_QUBIT_ALIGNMENT'):
            t0 = now_mu()   # record sequence begin time
            # initialize ion in S-1/2 state
            self.initialize_subsequence.run()
            self.cooling_subsequence.run()

            # prepare 729nm waveform for readout
            self.qubit.set_att_mu(self.att_qubit_mu)
            self.qubit.set_profile(self.profile_729_readout)
            self.qubit.cpld.io_update.pulse_mu(8)
            t1 = now_mu()   # record sequence end time

        # record exp shot period
        self._time_exp_shot_s = self.core.mu_to_seconds(t1 - t0) + self.time_readout_s

    @kernel(flags={"fast-math"})
    def run_main(self):
        # retrieve DMA handles for qubit alignment
        _handle_alignment = self.core_dma.get_handle('_QUBIT_ALIGNMENT')
        ion_state = (-1, 0, int64(0))   # stores ion state as (state, counts, elapsed_time_mu)

        # MAIN LOOP
        for num_rep in range(self.repetitions):
            d_state_prob = 0.

            # burst readout
            d_state_probs = [0.] * self.samples_per_point
            for idx in range(self.samples_per_point):
                self.core.break_realtime()

                ion_states = [int32(0)] * self.num_times

                # loop over target times
                for idx_time in range(self.num_times):
                    delay_mu(50000) # add slack - 50us

                    # prepare pulse shaping for target time
                    if self.enable_pulseshaping:
                        self.pulseshape_subsequence.configure(self.time_qubit_mu_list[idx_time])
                        delay_mu(25000)

                    # run qubit alignment sequence
                    self.core_dma.playback_handle(_handle_alignment)
                    if self.enable_pulseshaping:
                        self.pulseshape_subsequence.run()
                    else:
                        self.qubit.on()
                        delay_mu(self.time_qubit_mu_list[idx_time])
                        self.qubit.off()

                    # determine ion state & store results
                    self.readout_subsequence.run()
                    counts_res = self.readout_subsequence.fetch_count()

                    if ion_state[0] != -1:  # only store determinate results (indeterminate is -1)
                        self._filter_arr[idx_time].update_single(1 - ion_state[0]) # convert to D-state probability


                    if counts_res < 35:
                        ion_states[idx_time] = 1
                    else:
                        ion_states[idx_time] = 0

                d_state_prob = 0.
                for counts in ion_states:
                    d_state_prob += counts
                d_state_prob /= self.num_times


            # update dataset & clean up
            for idx_time in range(self.num_times):
                self.update_results(num_rep * (self.samples_per_point * self._time_exp_shot_s),
                                    idx_time,
                                    # self._filter_arr[idx_time].get_current() / self.samples_per_point)
                                    dstate_probability=d_state_prob)
            self.check_termination()    # check termination

    @rpc(flags={"async"})
    def update_results(self, time_dj: TFloat, idx_time: TInt32, dstate_probability: TFloat) -> TNone:
        """
        Overload the update_results function to allow real-time plot updating.
        :param time_dj: ***todo document***
        :param idx_time: ***todo document***
        :param dstate_probability: ***todo document***
        """
        # update datasets for broadcast
        self.mutate_dataset('temp.qubit_align.counts_x', self._result_iter, time_dj)
        self.mutate_dataset('temp.qubit_align.counts_y_{:d}'.format(idx_time), self._result_iter, dstate_probability)

        # update dataset for HDF5 storage
        # todo: slicing?
        # self.mutate_dataset('results', self._result_iter, array([time_dj, dstate_probability]))

        # update completion monitor
        if idx_time % self.num_times == 0:
            self.set_dataset('management.dynamic.completion_pct',
                             round(100. * self._result_iter / len(self.results), 3),
                             broadcast=True, persist=True, archive=False)
            self._result_iter += 1

