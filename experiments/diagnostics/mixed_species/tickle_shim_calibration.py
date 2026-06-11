from artiq.experiment import *
from artiq.coredevice import ad9910
from sipyco import pyon

from numpy import array, int32, int64, zeros
from scipy.interpolate import RegularGridInterpolator
from numpy import zeros, mean, linspace, array, sort, unique, transpose, reshape, take_along_axis, argsort
from numpy import interp, argmin, argmax, arange, shape, meshgrid, unravel_index, min, max, shape

from LAX_exp.language import *
from LAX_exp.system.subsequences import (
    InitializeQubit, SidebandCoolContinuousRAM, Readout, ReadoutAdaptive,
    RescueIon, QubitRAP
)

from LAX_exp.system.objects.PulseShaper import available_pulse_shapes
from LAX_exp.system.objects.dds_pulse_shaper import DDSPulseShaper
from LAX_exp.analysis.artiq_conversions import *
import matplotlib.pyplot as plt


class MixedSpeciesStrayFieldCalibration(LAXExperiment, Experiment):
    """
    Experiment: Mixed Species Stray Field Calibration

    Change shimming voltages and determine secular frequency parameters
    """
    name = 'Mixed Species Stray Field Calibration'
    kernel_invariants = {
        # subsequences & objects
        'initialize_subsequence', 'sidebandcool_subsequence', 'readout_subsequence',
        'rescue_subsequence', 'dds_pulse_shaper',

        # hardware values - tickle

        # hardware values - intensity servo
        'enable_servo_relock', 'time_servo_relock_mu',

        # configs
        'profile_729_SBC',
        'profile_tickle_RAM',
        'config_experiment_list',
    }

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions", NumberValue(default=50, precision=0, step=1, min=1, max=100000))

        # allocate relevant beam profiles
        self.profile_729_SBC = 0

        # allocate profiles for dds tickle
        self.profile_tickle_RAM = 0

        # get subsequences
        self.sidebandcool_subsequence = SidebandCoolContinuousRAM(
            self, profile_729=self.profile_729_SBC, profile_854=3,
            ram_addr_start_729=0, ram_addr_start_854=0, num_samples=200
        )
        self.initialize_subsequence = InitializeQubit(self)
        self.readout_subsequence = Readout(self)
        self.rescue_subsequence = RescueIon(self)

        # relevant devices
        self.setattr_device('qubit')
        self.setattr_device('pump')
        self.setattr_device('repump_qubit')
        self.setattr_device('dds_dipole')
        self.setattr_device('trap_dc')

        # set build arguments
        self._build_arguments_ion_parameters()
        self._build_arguments_tickle()
        self._build_arguments_shim()
        self._build_arguments_rap()
        self._build_arguments_intensity_servo()

        self.rap_subsequences = []
        self.att_reg_readout_rap_list = []

    def _build_arguments_ion_parameters(self):
        """
        Build arguments for ion frequencies
        """
        _argstr = "ion_parameters"

        self.setattr_argument("freq_secular_khz_list",
                             PYONValue([710]),
                             group=_argstr,
                             tooltip="Secular frequency used for tickle pulse (in kHz) applied via the urukul dds.")

        self.setattr_argument("freq_carrier_mhz", NumberValue(
            default=100.481,
            min=60., max=400, step=1,
            unit="MHz", scale=1, precision=6
        ), group=_argstr,
                              tooltip="Carrier frequency of the ion.\n"
                                      "Note: this is applied via the main doublepass DDS.\n")

    def _build_arguments_shim(self):
        """
        Build arguments for shim voltages
        """
        self.setattr_argument("h_shim_voltage_list", Scannable(
            default=[
                RangeScan(40, 60, 20, randomize=False),
                ExplicitScan([55, 60, 65]),
            ],
            global_min=0.0, global_max=120.0, global_step=1,
            unit="V", scale=1, precision=2
        ),
                              group='shim_voltages',
                              tooltip="The list of voltages to scan for the h shim.")

        self.setattr_argument("v_shim_voltage_list", Scannable(
            default=[
                RangeScan(40, 60, 20, randomize=False),
                ExplicitScan([55, 60, 65]),
            ],
            global_min=0.0, global_max=120.0, global_step=1,
            unit="V", scale=1, precision=2
        ),
                              group='shim_voltages',
                              tooltip="The list of voltages to scan for the v shim.")

    def _build_arguments_rap(self):
        _argstr = 'rap'

        self.setattr_argument("att_rap_db_list",
                             PYONValue([8]),
                             group=_argstr,
                             tooltip="rap att.")

        self.setattr_argument("freq_rap_dev_khz_list",
                             PYONValue([72]),
                             group=_argstr,
                             tooltip="rap att.")

        self.setattr_argument("time_rap_us_list",
                             PYONValue([1000]),
                             group=_argstr,
                             tooltip="rap att.")

    def _build_arguments_tickle(self):
        """
        Build core sweep arguments for the tickle pulse.
        """
        _argstr = "tickle"  # string to use for arguments

        # waveform - parameter sweeps
        self.setattr_argument("att_tickle_db_list",
                              PYONValue([25]),
                              group=_argstr,
                              tooltip="Attenuation to be used for the urukul channel used for generating the tickle.")
        self.setattr_argument("ampl_tickle_pct_list",
                              PYONValue([50,]),
                              group=_argstr,
                              tooltip='Amplitude of tickle pulse.')

        self.setattr_argument("time_tickle_us_list",
                              PYONValue([1000]),
                              group=_argstr,
                              tooltip="Time for the total pulse (including pulse shape).")


        self.setattr_argument("type_pulse_shape",
                              EnumerationValue(list(available_pulse_shapes.keys()), default='square'),
                              group=_argstr,
                              tooltip="Pulse shape type to be used.")

        self.setattr_argument("freq_tickle_detuning_khz_list", Scannable(
            default=[
                ExplicitScan([0]),
                CenterScan(0., 10., 0.001, randomize=True),
                RangeScan(-10, 10, 26, randomize=True),
            ],
            global_min = -1000, global_max=1000, global_step=0.001,
            unit="kHz", scale=1, precision=6),
                              group=_argstr,
                              tooltip="Detuning from secular frequency of tickle pulse (in kHz) applied via the urukul dds.")

    def _build_arguments_intensity_servo(self):
        """
        Build arguements for intensity servo hold between shots
        """
        _argstr = 'intensity_servo_relock'
        self.setattr_argument('enable_servo_relock', BooleanValue(default=False), group=_argstr,
                              tooltip='Enables the servo to relock the intensity of the 729 beam after every shot')
        self.setattr_argument('time_servo_relock_us', NumberValue(default=2000, precision=3, step=1, min=1,
                                                                  max=10000, scale=1., unit='us'),
                              group=_argstr, tooltip='Length of time to let the servo relock before each shot')

    def prepare_experiment(self):
        """
        Prepare & precompute experimental values.
        """

        ### Build Pulse Shaper
        self.dds_pulse_shaper = DDSPulseShaper(self, dds_target= self.dds_dipole.dds,
                                              ram_profile=self.profile_tickle_RAM,
                                              ram_addr_start=202, num_samples=250,
                                              ampl_max_pct=self.ampl_tickle_pct_list[0],
                                               pulse_shape=self.type_pulse_shape,
                                               phase_autoclear = 1)

        # run component preparation
        self._prepare_experiment_readout()

        freq_secular_ftw_list = self._prepare_experiment_ion_parameters()
        freq_tickle_detuning_ftw_list = self._prepare_experiment_tickle()

        self.h_shim_initial_voltage = self.trap_dc.get_h_shim_voltage()
        self.v_shim_initial_voltage = self.trap_dc.get_v_shim_voltage()

        self.h_shim_voltage_list = array(list(self.h_shim_voltage_list))
        self.v_shim_voltage_list = array(list(self.v_shim_voltage_list))

        # create experiment config
        self.config_experiment_list = create_experiment_config(

            # tickle sweeps
            freq_secular_ftw_list, freq_tickle_detuning_ftw_list,

            config_type=float, shuffle_config=True
        )

        # set timing for intensity servo between shots
        self.time_servo_relock_mu = self.core.seconds_to_mu(self.time_servo_relock_us * us)

    def _prepare_experiment_ion_parameters(self):
        """
        Prepare general ion parameters
        :return: list of carrier frequencies for cat and MS
        """
        self.freq_secular_ftw_list = [self.qubit.frequency_to_ftw(freq_secular_khz*kHz)
                                      for freq_secular_khz in self.freq_secular_khz_list]
        freq_secular_ftw_list = list(self.freq_secular_ftw_list)

        self.freq_carrier_ftw = self.qubit.frequency_to_ftw(self.freq_carrier_mhz * MHz)

        return freq_secular_ftw_list

    def _prepare_experiment_readout(self):
        """
        Prepare experiment values for state readout.
        """

        self.att_rap_mu = [att_to_mu(att_rap_db*dB) for att_rap_db in self.att_rap_db_list]
        self.rap_freq_dev_ftw_list = [self.qubit.frequency_to_ftw(freq_rap_dev_khz*kHz)
                                      for freq_rap_dev_khz in self.freq_rap_dev_khz_list]
        self.time_rap_mu_list = [self.core.seconds_to_mu(time_rap_us*us)
                                      for time_rap_us in self.time_rap_us_list]


        for idx, rap_profile in enumerate(self.freq_secular_khz_list):

            # attenuation register - readout (RAP): singlepasses set to default
            att_reg_readout_rap = 0x00000000 | (
                    (self.att_rap_mu[idx] << ((self.qubit.beam.chip_select - 4) * 8)) |
                    (self.qubit.att_singlepass0_default_mu << ((self.qubit.singlepass0.chip_select - 4) * 8)) |
                    (self.qubit.att_singlepass1_default_mu << ((self.qubit.singlepass1.chip_select - 4) * 8)) |
                    (self.qubit.att_singlepass2_default_mu << ((self.qubit.singlepass2.chip_select - 4) * 8))
            )

            self.rap_subsequences.append(QubitRAP(
                self, ram_profile=idx+1, ram_addr_start=202, num_samples=250,
                ampl_max_pct=50, pulse_shape="blackman"
            ))

            self.att_reg_readout_rap_list.append(att_reg_readout_rap)

    def _prepare_experiment_tickle(self):
        """
        Prepare general experiment values for the tickle pulse.
        :return: tuple of (freq_tickle_detuning_hz_list)
        """
        # convert values to convenience units
        self.att_tickle_mu_list = [att_to_mu(att_tickle_db * dB) for att_tickle_db in self.att_tickle_db_list]
        freq_tickle_detuning_ftw_list = [self.dds_pulse_shaper.dds_target.frequency_to_ftw(freq_tickle_detuning_khz*kHz)
                                         for freq_tickle_detuning_khz in self.freq_tickle_detuning_khz_list]

        self.time_tickle_mu_list = [self.core.seconds_to_mu(time_heating_us * us) for
                                time_heating_us in self.time_tickle_us_list]

        return freq_tickle_detuning_ftw_list


    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_experiment_list) * len(self.h_shim_voltage_list) * len(self.v_shim_voltage_list),
                5)

    '''
    MAIN SEQUENCE
    '''

    @kernel(flags={"fast-math"})
    def initialize_experiment(self) -> TNone:
        # record general subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.sidebandcool_subsequence.record_dma()
        self.readout_subsequence.record_dma()
        self.core.break_realtime()
        self.dds_pulse_shaper.dds_target.sw.off()
        self.qubit.off()


    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        # MAIN LOOP
        _loop_iter = 0
        self.trap_dc.h_shim_toggle(True)
        self.trap_dc.v_shim_toggle(True)
        self.core.break_realtime()

        for h_shim_voltage in self.h_shim_voltage_list:
            self.trap_dc.set_h_shim_voltage(h_shim_voltage)
            for v_shim_voltage in self.v_shim_voltage_list:
                self.trap_dc.set_v_shim_voltage(v_shim_voltage)
                self.core.break_realtime()
                for trial_num in range(self.repetitions):
                    for config_vals in self.config_experiment_list:
                        '''
                        PREPARE & CONFIGURE
                        '''
                        # extract values from config list
                        freq_secular_ftw = int32(config_vals[0])
                        freq_tickle_detuning_ftw = int32(config_vals[1])

                        '''
                        BEGIN MAIN SEQUENCE
                        '''
                        self.core.break_realtime()  # add slack for execution
                        delay_mu(125000)  # add even more slack lol
                        time_tickle_mu = 0

                        for idx in range(len(self.freq_secular_ftw_list)):
                            if freq_secular_ftw == self.freq_secular_ftw_list[idx]:

                                time_tickle_mu = self.time_tickle_mu_list[idx]

                                # configure tickle
                                self.dds_pulse_shaper.set_ampl_max_pct(self.ampl_tickle_pct_list[idx])
                                self.dds_pulse_shaper.sequence_initialize()
                                self.dds_pulse_shaper.dds_target.set_att_mu(self.att_tickle_mu_list[idx])
                                self.dds_pulse_shaper.dds_target.sw.off()
                                delay_mu(5000)
                                break

                        '''
                        Relock Intensity Servo
                        '''
                        if self.enable_servo_relock:
                            self.qubit.relock_intensity_servo(self.time_servo_relock_mu)


                        '''
                        INITIALIZE ION STATE
                        '''
                        # initialize ion in S-1/2 state & SBC to ground state
                        self.initialize_subsequence.run_dma()
                        self.sidebandcool_subsequence.run_dma()
                        delay_mu(8)

                        # set tickle frequency/phases
                        self.dds_pulse_shaper.dds_target.set_ftw(freq_secular_ftw + freq_tickle_detuning_ftw)
                        self.dds_pulse_shaper.dds_target.set_pow(0)

                        # set up config of shaped pulses to be fired for tickling
                        # also sets up phase autoclear
                        self.dds_pulse_shaper.configure_train(time_tickle_mu)
                        self.dds_pulse_shaper.dds_target.cpld.io_update.pulse_mu(8)
                        # for ururuk channel used for tickling keep RAM enabled but ensure we don't clear phase on io_update
                        self.dds_pulse_shaper.dds_target.set_cfr1(ram_enable=1, phase_autoclear=0,
                                                                  ram_destination=ad9910.RAM_DEST_ASF)
                        self.dds_pulse_shaper.dds_target.cpld.io_update.pulse_mu(8)
                        delay_mu(100)

                        '''
                        TICKLE PULSE
                        '''
                        self.dds_pulse_shaper.run_train_single()
                        delay_mu(time_tickle_mu)

                        '''
                        READ OUT & STORE RESULTS
                        '''
                        self.pulse_readout_rap(freq_secular_ftw)

                        # read out fluorescence & clean up loop
                        self.readout_subsequence.run_dma()
                        counts_res = self.readout_subsequence.fetch_count()

                        # cleanup dds_pulse_shaper
                        self.dds_pulse_shaper.sequence_cleanup()

                        # store results
                        self.update_results(freq_secular_ftw,
                                            counts_res,
                                            freq_tickle_detuning_ftw,
                                            h_shim_voltage,
                                            v_shim_voltage
                                    )

                # check termination more frequently in case reps are low
                if _loop_iter % 100 == 0:
                    self.check_termination()
                _loop_iter += 1

            # rescue ion as needed & support graceful termination
            self.check_termination()

        '''
        Reset shims
        '''
        self.trap_dc.set_h_shim_voltage(self.h_shim_initial_voltage)
        self.trap_dc.set_v_shim_voltage(self.v_shim_initial_voltage)

    @kernel(flags={"fast-math"})
    def pulse_readout_rap(self, sec_freq_ftw: TInt32) -> TNone:
        """
        Run a RAP readout pulse.
        """
        # set up relevant beam waveforms
        # find right sec freq
        rap_idx = 0
        for idx in range(len(self.freq_secular_ftw_list)):
            if sec_freq_ftw == self.freq_secular_ftw_list[idx]:
                rap_idx = idx
                break

        sec_freq_ftw = self.freq_secular_ftw_list[rap_idx]
        self.rap_subsequences[rap_idx].configure(self.time_rap_mu_list[rap_idx], self.freq_carrier_ftw - (sec_freq_ftw >> 1),
                                                 self.rap_freq_dev_ftw_list[rap_idx])
        delay_mu(50000)

        self.qubit.off()
        self.qubit.singlepass0_on()
        self.qubit.singlepass1_off()
        self.qubit.singlepass2_off()
        self.qubit.cpld.set_all_att_mu(self.att_reg_readout_rap_list[rap_idx])
        # run RAP readout pulse
        # run RAP turns on qubit
        self.rap_subsequences[rap_idx].run_rap(self.time_rap_mu_list[rap_idx])

    def analyze_experiment(self):
        """
        Analyze data and determine optimal voltages to operate at
        """

        # determine the nbar for a range of coherent displacements and determine conversion for rsb rap
        alphas = linspace(0, 5, 1000)
        nbars = [alpha_to_nbar(alpha) for alpha in alphas]
        phonon_conversions = [phonon_conversion(alpha) for alpha in alphas]

        # get data
        data = array(self.results)

        # grab relevant data
        sec_freqs = ftw_to_frequency_khz(data[:, 0])
        counts_arr = array(data[:, 1])
        tickle_detunings = ftw_to_frequency_khz(data[:, 2])
        h_shim_voltages = data[:, 3]
        v_shim_voltages = data[:, 4]

        # get number of detuning points and their uniques values
        num_freq_points = len(unique(tickle_detunings))
        unique_tickle_freqs = sort(unique(tickle_detunings))
        unique_sec_freqs = sort(unique(sec_freqs))
        unique_h_shim_voltages = sort(unique(h_shim_voltages))
        unique_v_shim_voltages = sort(unique(v_shim_voltages))
        num_sec_freqs = len(unique_sec_freqs)

        reps = self.repetitions

        # get the lowest threshold to determine when all ions when dark
        threshold_list = array(findThresholdScikit(counts_arr, 60))
        threshold = min(threshold_list)

        # setup storage arrays for fitted secular frequency at a given h shim voltage and v shim voltage
        secular_freq_detuning_s_kHz_arr = zeros((len(unique_h_shim_voltages), len(unique_v_shim_voltages)))
        secular_freq_detuning_s_err_kHz_arr = zeros((len(unique_h_shim_voltages), len(unique_v_shim_voltages)))

        # generate detuning values for fitting to tickle linescan
        fit_x = linspace(min(unique_tickle_freqs), max(unique_tickle_freqs), 1000)

        # loop through all h shim and v shim voltages values
        for i, unique_h_shim_voltage in enumerate(unique_h_shim_voltages):
            for j, unique_v_shim_voltage in enumerate(unique_v_shim_voltages):

                # get photon counts and tickle detunings at the selected h shim AND v shim voltage
                counts_tmp = counts_arr[
                    (v_shim_voltages == unique_v_shim_voltage) & (h_shim_voltages == unique_h_shim_voltage)]
                tickle_detunings_tmp = tickle_detunings[
                    (v_shim_voltages == unique_v_shim_voltage) & (h_shim_voltages == unique_h_shim_voltage)]

                # reshape counts and detunings
                counts_reshaped = reshape(counts_tmp, (-1, num_freq_points, reps))
                tickle_detunings_reshaped = reshape(tickle_detunings_tmp, (-1, num_freq_points, reps))

                for idx, counts_tmp_arr in enumerate(counts_reshaped):

                    if not (counts_arr < threshold).all():
                        tickle_detunings_arr = tickle_detunings_reshaped[idx]

                        tickle_detunings_arr = transpose(tickle_detunings_arr)
                        counts_tmp_arr = transpose(counts_tmp_arr)

                        # sort counts based on tickle detuning
                        counts_tmp_arr_sorted = take_along_axis(counts_tmp_arr, argsort(tickle_detunings_arr), axis=1)

                        # use threshold to determine dark state of ion
                        state_vals = counts_tmp_arr_sorted < threshold
                        state_vals_ave = sum(state_vals, 0) / reps
                        # convert state vals to phonons from rsb rap
                        phonons = interp(state_vals_ave, phonon_conversions, nbars)

                        # set guess params for lineshape
                        max_phonon = max(phonons)
                        peak_freq = unique_tickle_freqs[argmax(phonons)]
                        linewidth_scalar = 1
                        contrast_loss = 0
                        try:
                            # exctract fitted parameters from fit
                            popt_s, pcov_s = curve_fit(sinc_fit, unique_tickle_freqs, phonons,
                                                       p0=[max_phonon, peak_freq, linewidth_scalar, contrast_loss],
                                                       maxfev=15000)

                            perr_s = sqrt(diag(pcov_s))

                            # add fitted parameters to storage arrays
                            secular_freq_detuning_s_kHz_arr[i, j] = popt_s[1]
                            secular_freq_detuning_s_err_kHz_arr[i, j] = perr_s[1]

                        except:
                            # print(f'Could not process data')
                            pass
                        # from numpy.random import randint
                        #
                        # secular_freq_detuning_s_kHz_arr[i, j] = (unique_h_shim_voltage - mean(
                        #     unique_h_shim_voltages)) ** 2 + (unique_v_shim_voltage - mean(unique_v_shim_voltages)) ** 2
                        # secular_freq_detuning_s_err_kHz_arr[i, j] = 1


        # subtract mean value of secular frequency to focus on offset
        sec_freq_drift = array(secular_freq_detuning_s_kHz_arr) - mean(secular_freq_detuning_s_kHz_arr)


        #
        if len(unique_v_shim_voltages) > 1 and len(unique_h_shim_voltages) == 1:
            plotting_results = {'x': unique_v_shim_voltages.flatten(),
                                'y': sec_freq_drift.flatten(),
                                'subplot_titles': f'Secular Freq Shimming Calibration',
                                'subplot_x_labels': 'H Shim Voltage (V)',
                                'subplot_y_labels': 'Secular Freq Offset (kHz)',
                                'rid': self.scheduler.rid, }
            projection_3d = False

        elif len(unique_v_shim_voltages) == 1 and len(unique_h_shim_voltages) > 1:
            plotting_results = {'x': unique_h_shim_voltages.flatten(),
                                'y': sec_freq_drift.flatten(),
                                'subplot_titles': f'Secular Freq Shimming Calibration',
                                'subplot_x_labels': 'H Shim Voltage (V)',
                                'subplot_y_labels': 'Secular Freq Offset (kHz)',
                                'rid': self.scheduler.rid, }
            projection_3d = False

        else:
            h_shim_voltages_meshed, v_shim_voltages_meshed = meshgrid(unique_h_shim_voltages, unique_v_shim_voltages, indexing='ij')

            ### interpolator ###
            grid_interp = RegularGridInterpolator((unique_h_shim_voltages, unique_v_shim_voltages), sec_freq_drift)
            fit_h_shim_voltages = linspace(min(unique_h_shim_voltages), max(unique_h_shim_voltages), 1000)
            fit_v_shim_voltages = linspace(min(unique_v_shim_voltages), max(unique_v_shim_voltages), 1000)
            fit_h_shim_voltages_meshed, fit_v_shim_voltages_meshed = meshgrid(fit_h_shim_voltages, fit_v_shim_voltages, indexing='ij')
            interpolation_data = grid_interp((fit_h_shim_voltages_meshed, fit_v_shim_voltages_meshed))
            optimal_voltages = unravel_index(argmin(interpolation_data), interpolation_data.shape)

            # format dictionary for applet plotting
            plotting_results = {'x': h_shim_voltages_meshed,
                                'y': v_shim_voltages_meshed,
                                'z': sec_freq_drift,
                                'fit_x': fit_h_shim_voltages_meshed,
                                'fit_y': fit_v_shim_voltages_meshed,
                                'fit_z': interpolation_data,
                                'subplot_titles': f'Secular Freq Shimming Calibration',
                                'subplot_x_labels': 'H Shim Voltage (V)',
                                'subplot_y_labels': 'V Shim Voltage (V)',
                                'rid': self.scheduler.rid,}
            projection_3d = True

        self.set_dataset('temp.plotting.results_mixed_species_tickle_shim_calibration', pyon.encode(plotting_results), broadcast=True)
        self.ccb.issue("create_applet", f"Mixed Species Tickle Shim Calibration",
                       '$python -m LAX_exp.applets.plot_matplotlib temp.plotting.results_mixed_species_tickle_shim_calibration'
                       f' --num-subplots 1 --projection_3d {projection_3d}',
                       group=['plotting', 'diagnostics', 'mixed_species'])

        return data


