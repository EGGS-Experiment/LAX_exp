import numpy as np
from artiq.experiment import *
from artiq.coredevice import ad9910
from sipyco import pyon

from scipy.optimize import curve_fit
import numpy as np
from numpy import array, int32, linspace, sort, unique, interp, meshgrid, unravel_index
from scipy.interpolate import griddata
from collections.abc import Iterable

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
        'rescue_subsequence', 'dds_pulse_shaper', 'rap_subsequence',

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

    def _build_arguments_ion_parameters(self):
        """
        Build arguments for ion frequencies
        """
        _argstr = "ion_parameters"

        self.setattr_argument("freq_secular_khz",
                             NumberValue(default=710,
                                         min=60., max=2500,
                                         step=1, unit="kHz",
                                         scale=1, precision=3),
                             group=_argstr,
                             tooltip="Secular frequency used for tickle pulse (in kHz) applied via the urukul dds.")

        self.setattr_argument("freq_carrier_mhz", NumberValue(
            default=100.481,
            min=60., max=400, step=1,
            unit="MHz", scale=1, precision=6
        ),
                              group=_argstr,
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
        """
        Build arguements for rap readout
        """
        _argstr = 'rap'

        self.setattr_argument("att_rap_db",
                             NumberValue(default=8,
                                         min=8, max=31.5,
                                         step=0.5, scale=1,
                                         precision=1, unit="dB"),
                             group=_argstr,
                             tooltip="rap att.")

        self.setattr_argument("freq_rap_dev_khz",
                              NumberValue(default=72,
                                          min=10, max=300,
                                          step=0.1, scale=1,
                                          precision=2, unit="kHz"),
                             group=_argstr,
                             tooltip="rap att.")

        self.setattr_argument("time_rap_us",
                              NumberValue(default=1000,
                                          min=50, max=5000,
                                          step=0.1, scale=1,
                                          precision=3, unit="us"),
                             group=_argstr,
                             tooltip="rap att.")

        self.setattr_argument("freq_rap_center_mhz",
                              NumberValue(default=100.755, precision=6, step=1e-2, min=1, max=200, unit="MHz",
                                          scale=1.),
                              group=_argstr)

    def _build_arguments_tickle(self):
        """
        Build arguments for the tickle pulse.
        """
        _argstr = "tickle"  # string to use for arguments

        # waveform - parameter sweeps
        self.setattr_argument("att_tickle_db",
                              NumberValue(
                                  default=25.,
                                  min=5., max=31.5,
                                  step=0.5, scale=1,
                                  unit='dB', precision=1,
                              ),
                              group=_argstr,
                              tooltip="Attenuation to be used for the urukul channel used for generating the tickle.")
        self.setattr_argument("ampl_tickle_pct",
                              NumberValue(
                                  default=50.,
                                  min=0.01, max=50,
                                  step=0.5, scale=1,
                                  unit='%', precision=3,
                              ),
                              group=_argstr,
                              tooltip='Amplitude of tickle pulse.')

        self.setattr_argument("time_tickle_us",
                              NumberValue(
                                  default=1000.,
                                  min=2., max=100000.,
                                  step=0.5, scale=1,
                                  unit='us', precision=3,
                              ),
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

        # Build Pulse Shaper
        self.dds_pulse_shaper = DDSPulseShaper(self, dds_target= self.dds_dipole.dds,
                                              ram_profile=self.profile_tickle_RAM,
                                              ram_addr_start=202, num_samples=250,
                                              ampl_max_pct=self.ampl_tickle_pct,
                                               pulse_shape=self.type_pulse_shape,
                                               phase_autoclear = 1)

        # run component preparation
        self._prepare_experiment_readout()

        self._prepare_experiment_ion_parameters()
        freq_tickle_detuning_ftw_list = self._prepare_experiment_tickle()
        # create experiment config
        self.config_experiment_list = create_experiment_config(

            # tickle sweeps
            freq_tickle_detuning_ftw_list,

            config_type=float, shuffle_config=True
        )

        # set timing for intensity servo between shots
        self.time_servo_relock_mu = self.core.seconds_to_mu(self.time_servo_relock_us * us)

        # prepare shim voltages list
        self.h_shim_voltage_list = array(list(self.h_shim_voltage_list))
        self.v_shim_voltage_list = array(list(self.v_shim_voltage_list))

        # set storage variables
        self.h_shim_initial_voltage = self.trap_dc.get_h_shim_voltage()
        self.v_shim_initial_voltage = self.trap_dc.get_v_shim_voltage()


    def _prepare_experiment_ion_parameters(self):
        """
        Prepare general ion parameters
        """
        self.freq_secular_ftw = self.qubit.frequency_to_ftw(self.freq_secular_khz*kHz)

        self.freq_carrier_ftw = self.qubit.frequency_to_ftw(self.freq_carrier_mhz * MHz)

    def _prepare_experiment_readout(self):
        """
        Prepare experiment values for state readout.
        """

        self.att_rap_mu = att_to_mu(self.att_rap_db*dB)
        self.rap_freq_dev_ftw = self.qubit.frequency_to_ftw(self.freq_rap_dev_khz*kHz)
        self.time_rap_mu = self.core.seconds_to_mu(self.time_rap_us*us)
        self.freq_rap_center_ftw = self.qubit.frequency_to_ftw(self.freq_rap_center_mhz * MHz)

        # attenuation register - readout (RAP): singlepasses set to default
        self.att_reg_readout_rap = 0x00000000 | (
                (self.att_rap_mu << ((self.qubit.beam.chip_select - 4) * 8)) |
                (self.qubit.att_singlepass0_default_mu << ((self.qubit.singlepass0.chip_select - 4) * 8)) |
                (self.qubit.att_singlepass1_default_mu << ((self.qubit.singlepass1.chip_select - 4) * 8)) |
                (self.qubit.att_singlepass2_default_mu << ((self.qubit.singlepass2.chip_select - 4) * 8))
        )

        self.rap_subsequence = QubitRAP(
            self, ram_profile=1, ram_addr_start=202, num_samples=250,
            ampl_max_pct=50, pulse_shape="blackman"
        )


    def _prepare_experiment_tickle(self):
        """
        Prepare general experiment values for the tickle pulse.
        Returns:
            A tuple containing
                - freq_tickle_detuniing_ftw_list (iterable): a list of tickle detuning in ftw
        """
        # convert values to convenience units
        self.att_tickle_mu = att_to_mu(self.att_tickle_db * dB)
        freq_tickle_detuning_ftw_list = [self.dds_pulse_shaper.dds_target.frequency_to_ftw(freq_tickle_detuning_khz*kHz)
                                         for freq_tickle_detuning_khz in self.freq_tickle_detuning_khz_list]

        self.time_tickle_mu = self.core.seconds_to_mu(self.time_tickle_us * us)

        return freq_tickle_detuning_ftw_list


    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_experiment_list) * len(self.h_shim_voltage_list) * len(self.v_shim_voltage_list),
                4)

    '''
    MAIN SEQUENCE
    '''

    @kernel(flags={"fast-math"})
    def initialize_experiment(self) -> TNone:
        """
        Initialize the experimental hardware
        """
        # record general subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.sidebandcool_subsequence.record_dma()
        self.readout_subsequence.record_dma()
        self.core.break_realtime()
        self.dds_pulse_shaper.dds_target.sw.off()

        # configure tickle
        self.dds_pulse_shaper.set_ampl_max_pct(self.ampl_tickle_pct)
        self.dds_pulse_shaper.sequence_initialize()
        self.dds_pulse_shaper.dds_target.set_att_mu(self.att_tickle_mu)
        self.dds_pulse_shaper.dds_target.sw.off()
        self.qubit.off()

        # get shim voltages
        self.h_shim_initial_voltage = self.trap_dc.get_h_shim_voltage()
        self.v_shim_initial_voltage = self.trap_dc.get_v_shim_voltage()


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
                delay_mu(100000000)
                self.core.break_realtime()
                for trial_num in range(self.repetitions):
                    for config_vals in self.config_experiment_list:
                        '''
                        PREPARE & CONFIGURE
                        '''
                        # extract values from config list
                        freq_tickle_detuning_ftw = int32(config_vals[0])

                        '''
                        BEGIN MAIN SEQUENCE
                        '''
                        self.core.break_realtime()  # add slack for execution
                        delay_mu(125000)  # add even more slack lol

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
                        self.dds_pulse_shaper.dds_target.set_ftw(self.freq_secular_ftw + freq_tickle_detuning_ftw)
                        self.dds_pulse_shaper.dds_target.set_pow(0)

                        # set up config of shaped pulses to be fired for tickling
                        # also sets up phase autoclear
                        self.dds_pulse_shaper.configure_train(self.time_tickle_mu)
                        self.dds_pulse_shaper.dds_target.cpld.io_update.pulse_mu(8)
                        # for ururuk channel used for tickling keep RAM enabled but ensure we don't clear phase on io_update
                        self.dds_pulse_shaper.dds_target.set_cfr1(ram_enable=1, phase_autoclear=0,
                                                                  ram_destination=ad9910.RAM_DEST_ASF)
                        self.dds_pulse_shaper.dds_target.cpld.io_update.pulse_mu(8)
                        delay_mu(100)

                        '''
                        TICKLE PULSE
                        '''
                        self.dds_pulse_shaper.run_train_all_dds()
                        delay_mu(self.time_tickle_mu)

                        '''
                        READ OUT & STORE RESULTS
                        '''
                        self.pulse_readout_rap()

                        # read out fluorescence & clean up loop
                        self.readout_subsequence.run_dma()
                        counts_res = self.readout_subsequence.fetch_count()

                        # cleanup dds_pulse_shaper
                        self.dds_pulse_shaper.sequence_cleanup()

                        # store results
                        self.update_results(freq_tickle_detuning_ftw,
                                            counts_res,
                                            h_shim_voltage,
                                            v_shim_voltage
                                    )

                        # check termination more frequently in case reps are low
                        if _loop_iter % 200 == 0:
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
    def pulse_readout_rap(self) -> TNone:
        """
        Run a RAP readout pulse.
        """
        # set up relevant beam waveforms
        self.rap_subsequence.configure(self.time_rap_mu, self.freq_rap_center_ftw,
                                                 self.rap_freq_dev_ftw)

        self.qubit.off()
        self.qubit.singlepass0_on()
        self.qubit.singlepass1_off()
        self.qubit.singlepass2_off()
        self.qubit.cpld.set_all_att_mu(self.att_reg_readout_rap)
        # run RAP readout pulse
        # run RAP turns on qubit
        self.rap_subsequence.run_rap(self.time_rap_mu)

    @kernel(flags={"fast-math"})
    def cleanup_experiment(self) -> TNone:
        self.trap_dc.set_h_shim_voltage(self.h_shim_initial_voltage)
        self.trap_dc.set_v_shim_voltage(self.v_shim_initial_voltage)


    def analyze_experiment(self):
        """
        Analyze data and determine optimal voltages to operate at
        """

        # get data
        data = array(self.results)

        # grab relevant data
        counts_arr = array(data[:, 1])
        tickle_detunings = ftw_to_frequency_khz(data[:, 0])
        h_shim_voltages = array(data[:, 2])
        v_shim_voltages = array(data[:, 3])

        # get number of detuning points and their uniques values
        num_freq_points = len(unique(tickle_detunings))
        unique_tickle_freqs = sort(unique(tickle_detunings))
        unique_h_shim_voltages = sort(unique(h_shim_voltages))
        unique_v_shim_voltages = sort(unique(v_shim_voltages))

        reps = self.repetitions

        # get the lowest threshold to determine when all ions when dark
        threshold_list = array(findThresholdScikit(counts_arr, 60))
        dark_state_threshold = min(threshold_list)

        # setup storage arrays for fitted secular frequency at a given h shim voltage and v shim voltage
        secular_detuning_opt = np.zeros((len(unique_h_shim_voltages), len(unique_v_shim_voltages)))
        secular_detuning_err = np.zeros((len(unique_h_shim_voltages), len(unique_v_shim_voltages)))

        # loop through all h shim and v shim voltages values
        for i, unique_h_shim_voltage in enumerate(unique_h_shim_voltages):
            for j, unique_v_shim_voltage in enumerate(unique_v_shim_voltages):

                loop_iter = i*len(unique_v_shim_voltages) + j

                mask = (
                        (v_shim_voltages == unique_v_shim_voltage)
                        & (h_shim_voltages == unique_h_shim_voltage)
                )

                unique_detuning_khz, phonons = self._get_processed_data(counts_arr, tickle_detunings, mask, dark_state_threshold)


                # add fitted parameters to storage arrays
                opt_detuning, opt_detuning_err = self._fit_plot_tickle_scan(unique_detuning_khz, phonons,
                                                                            unique_h_shim_voltage, unique_v_shim_voltage, loop_iter)
                secular_detuning_opt[i, j] = opt_detuning
                secular_detuning_err[i, j] = opt_detuning_err

        # subtract mean value of secular frequency to focus on offset (nanmean ignores nans)
        sec_freq_detuning = array(secular_detuning_opt) - np.nanmean(secular_detuning_opt)
        issue_final_plot = True

        if len(unique_v_shim_voltages) > 1 and len(unique_h_shim_voltages) == 1:
            plotting_results = self._format_one_dim_plotting_results(unique_v_shim_voltages, sec_freq_detuning, plot_xlabel = 'V Shim Voltage (V)')
            projection_3d = False

        elif len(unique_v_shim_voltages) == 1 and len(unique_h_shim_voltages) > 1:
            plotting_results = self._format_one_dim_plotting_results(unique_h_shim_voltages, sec_freq_detuning, plot_xlabel = 'H Shim Voltage (V)')
            projection_3d = False

        elif len(unique_v_shim_voltages) > 1 and len(unique_h_shim_voltages) > 1:
            plotting_results = self._plot_shim_calibration(unique_h_shim_voltages, unique_v_shim_voltages, sec_freq_detuning)
            projection_3d = True
        else:
            issue_final_plot = False

        if issue_final_plot:
            self.create_matplotlib_applet(plotting_results,
                                      name=f'Shim Calibration',
                                      group=['plotting', 'diagnostics', 'mixed_species'],
                                      projection_3d=projection_3d)

        return data

    def _get_processed_data(self, counts_arr: Iterable[int], tickle_detunings: Iterable[float],
                            mask, threshold: float):
        """
        Extract and process the data for a given horizontal shim voltage and vertical shim voltage configuration

        Args:
            counts_arr: number of photons seen by PMT for each shot of the experiment
            tickle_detunings (iterable): detunings in kHz of tickle tone from secular frequency given to ARTIQ experiment
            mask: list containing indices which contain the data for the current configuration
            threshold (float): boundary to determine when ion (or all ions) is (are) in the dark state (D_{5/2})

        Returns:
            A tuple containing two elements:
                - unique_detuning_khz (iterable): list of unique detunings, in khz, used for tickle tone
                - phonons (iterable): list of average motional quanta in mode of ion's motional harmonic oscillator
        """

        # determine the nbar for a range of coherent displacements and determine conversion for rsb rap
        alphas = linspace(0, 5, 1000)
        nbars = [alpha_to_nbar(alpha) for alpha in alphas]
        phonon_conversions = [phonon_conversion(alpha) for alpha in alphas]

        # get photon counts and tickle detunings at the selected h shim AND v shim voltage
        counts_tmp = counts_arr[mask]
        tickle_detunings_tmp = tickle_detunings[mask]
        state_vals = np.array(counts_tmp < threshold, dtype=np.int32)

        # group all shots by identical detuning
        unique_detuning_khz, unique_detuning_idx = np.unique(tickle_detunings_tmp, return_inverse=True)

        # create storage array for average state vals
        dark_probability = np.zeros(len(unique_detuning_khz))

        for idx in range(len(unique_detuning_khz)):
            dark_probability[idx] = np.mean(state_vals[unique_detuning_idx == idx])

        # convert state vals to phonons from rsb rap
        phonons = interp(dark_probability , phonon_conversions, nbars)
        return unique_detuning_khz, phonons

    def _fit_plot_tickle_scan(self, unique_detuning_khz: Iterable[float], phonons: Iterable[float],
                              h_voltage: float, v_voltage: float, loop_iter: int):
        """
        Fit and plot a tickle scan taken at a single horizontal shim voltage and vertical shim voltage configuration

        Args:
            unique_detuning_khz (iterable): detuning (in khz) from secular frequency provided to ARTIQ experiment
            phonons (iterable): average motional quantum in the harmonic oscillator mode
            h_voltage (float): horizontal shim voltage
            v_voltage (float): vertical shim voltage
            loop_iter (int): loop iteration

        Returns:
            A tuple containing two elements:
                - opt_detuning (float): secular frequency detuning (from secular frequency provided to ARTIQ experiment) in kHz from fit
                - opt_detuning (float): error from fit in calculating opt_detuning
        """
        try:
            # exctract fitted parameters from fit
            fitter = fitSincGeneric()
            popt_s, perr_s, _ = fitter.fit(unique_detuning_khz, phonons)

            fit_x = np.linspace(min(unique_detuning_khz), max(unique_detuning_khz), 1000)
            fit_y = fitter.fit_func(fit_x, *popt_s)
            opt_detuning = popt_s[1]
            err_detuning = popt_s[2]

        except Exception as exc:
            print(f"Could not process h={h_voltage}, v={v_voltage}: {exc}")
            opt_detuning = np.nan
            err_detuning = np.nan
            fit_x = None
            fit_y = None

        if loop_iter % 5 == 0:
            # format dictionary for applet plotting
            plotting_results = {'x': unique_detuning_khz.flatten(),
                                'y': phonons.flatten(),
                                'fit_x': None if fit_x is None else fit_x.flatten(),
                                'fit_y': None if fit_y is None else fit_y.flatten(),
                                'subplot_titles': f'Secular Freq Shimming Calibration H\\V: {np.round(h_voltage,2)}\\{np.round(v_voltage,2)})',
                                'subplot_x_labels': 'Secular Frequency Detuning (kHz)',
                                'subplot_y_labels': 'Phonon Number',
                                'rid': self.scheduler.rid}
            projection_3d = False
            self.create_matplotlib_applet(plotting_results,
                                          name=f'Shim Calibration Tickle Scan {np.round(h_voltage,2)}_{np.round(v_voltage, 2)}',
                                          group=['plotting', 'diagnostics', 'mixed_species'],
                                          projection_3d=projection_3d)

        return opt_detuning, err_detuning

    def _plot_shim_calibration(self, h_voltages: Iterable[float], v_voltages: Iterable[float],
                               sec_freq_detuning: Iterable[float]):
        """
        Helper function to plot a 3D plot displaying secular frequency as a function of the shim voltage

        Args:
            h_voltages (iterable): voltages applied to the horizontal shim rod
            v_voltages (iterable): voltages applied to the vertical shim rod
            sec_freq_detuning (iterable): detuning of secular frequency from mean value

        Returns:
            plotting_results (dict): dictionary containing information used to create plotting applet
        """
        h_shim_voltages_meshed, v_shim_voltages_meshed = meshgrid(h_voltages,
                                                                  v_voltages,
                                                                  indexing='ij')

        # find valid points
        valid = np.isfinite(sec_freq_detuning)
        points = np.column_stack([h_shim_voltages_meshed[valid],
                                  v_shim_voltages_meshed[valid]])

        values = sec_freq_detuning[valid]

        # dense grid for plotting
        fit_h_shim_voltages = linspace(
            min(h_voltages),
            max(h_voltages),
            1000)
        fit_v_shim_voltages = linspace(
            min(v_voltages),
            max(v_voltages),
            1000)

        fit_h_shim_voltages_meshed, fit_v_shim_voltages_meshed = meshgrid(
            fit_h_shim_voltages,
            fit_v_shim_voltages,
            indexing='ij')

        # use griddata function for scan where some tickle scans fail to fit
        try:
            interpolation_data = griddata(points,
                                          values,
                                          (fit_h_shim_voltages_meshed, fit_v_shim_voltages_meshed),
                                          method='cubic')
        except Exception as exc:
            interpolation_data = griddata(points,
                                          values,
                                          (fit_h_shim_voltages_meshed, fit_v_shim_voltages_meshed),
                                          method='linear')


        if not np.isfinite(interpolation_data).any():
            print("No valid interpolation data for shim calibration")

        else:
            optimal_min_h_voltages_ind, optimal_min_v_voltages_ind = unravel_index(
                np.nanargmin(interpolation_data),
                interpolation_data.shape)
            optimal_min_v_voltage = fit_v_shim_voltages[optimal_min_v_voltages_ind]
            optimal_min_h_voltage = fit_h_shim_voltages[optimal_min_h_voltages_ind]
            print(f"The minimium secular frequency occurs at H\\V: {np.round(optimal_min_h_voltage,2)}\\{np.round(optimal_min_v_voltage,2)}")

            optimal_max_h_voltages_ind, optimal_max_v_voltages_ind = unravel_index(
                np.nanargmax(interpolation_data),
                interpolation_data.shape)
            optimal_max_v_voltage = fit_v_shim_voltages[optimal_max_v_voltages_ind]
            optimal_max_h_voltage = fit_h_shim_voltages[optimal_max_h_voltages_ind]
            print(f"The maximium secular frequency occurs at H\\V: {np.round(optimal_max_h_voltage,2)}\\{np.round(optimal_max_v_voltage,2)}")


        # format dictionary for applet plotting
        plotting_results = {'x': h_shim_voltages_meshed,
                            'y': v_shim_voltages_meshed,
                            'z': sec_freq_detuning,
                            'fit_x': fit_h_shim_voltages_meshed,
                            'fit_y': fit_v_shim_voltages_meshed,
                            'fit_z': interpolation_data,
                            'subplot_titles': f'Secular Freq Shimming Calibration',
                            'subplot_x_labels': 'H Shim Voltage (V)',
                            'subplot_y_labels': 'V Shim Voltage (V)',
                            'subplot_z_labels': 'Secular Frequency Offset (kHz)',
                            'rid': self.scheduler.rid}

        return plotting_results

    def _format_one_dim_plotting_results(self, voltages, sec_freq_detuning, plot_xlabel = ''):
        plotting_results = {'x': voltages.flatten(),
                                'y': sec_freq_detuning.flatten(),
                                'subplot_titles': f'Secular Freq Shimming Calibration',
                                'subplot_x_labels': plot_xlabel,
                                'subplot_y_labels': 'Secular Freq Offset (kHz)',
                                'rid': self.scheduler.rid}

        return plotting_results