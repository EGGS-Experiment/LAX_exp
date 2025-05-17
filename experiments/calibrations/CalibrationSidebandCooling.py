import numpy as np
from artiq.experiment import *

from LAX_exp.analysis import *
from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import InitializeQubit, SidebandReadout, Readout, RescueIon

from itertools import product
from artiq.coredevice import ad9910

_VALID_BEAM_SCANS = [
    'doppler_397', 'doppler_866',
    'readout_397'
]


class CalibrationSidebandCooling(LAXExperiment, Experiment):
    """
    Calibration: Sideband Cooling

    Continuous sideband cooling but with scannable parameters for optimization.
    Non-RAM based to simplify parameter scanning.
    """
    name = 'Calibration Sideband Cooling'
    kernel_invariants = {
        # hardware values
        'att_sbc_mu', 'dds_beam', 'beam_update_profile',
        'num_modes', 'sbc_config_base_list', 'sbc_mode_time_frac_list',

        # subsequences
        'initialize_subsequence', 'sidebandreadout_subsequence', 'readout_subsequence', 'rescue_subsequence',

        # configs
        'profile_729_readout', 'profile_729_SBC', 'profile_854_SBC', 'config_experiment_list'
    }

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions",    NumberValue(default=40, precision=0, step=1, min=1, max=100000))

        # allocate profiles on 729nm for different subsequences
        self.profile_729_readout =  0
        self.profile_729_SBC =      1
        self.profile_854_SBC =      5

        # SBC - base config
        self.setattr_argument("sideband_cooling_config_list", PYONValue({100.7555: [26., 5.], 100.455: [37., 5.], 100.315: [37., 5.]}),
                              tooltip="{freq_mode_mhz: [sbc_mode_pct_per_cycle, ampl_quench_mode_pct]}", group='SBC.base')
        self.setattr_argument("sideband_cycles_continuous", NumberValue(default=10, precision=0, step=1, min=1, max=10000, unit="cycles", scale=1),
                              tooltip="number of times to loop over the SBC configuration sequence", group='SBC.base')
        self.setattr_argument("att_sidebandcooling_continuous_db",  NumberValue(default=8, precision=1, step=0.5, min=8, max=31.5, unit="dB", scale=1),
                              group='SBC_RAM.continuous', tooltip="Attenuation of 729nm beam for SBC.")
        self.setattr_argument("time_sbc_us_list",   Scannable(
                                                        default=[
                                                            ExplicitScan([4000, 5000]),
                                                            RangeScan(100, 2000, 100, randomize=True),
                                                            CenterScan(3.05, 5., 0.1, randomize=True),
                                                        ],
                                                        global_min=50, global_max=100000, global_step=1,
                                                        unit="us", scale=1, precision=5
                                                    ), group="SBC.base")
        self.setattr_argument("time_per_spinpol_us_list",   Scannable(
                                                        default=[
                                                            RangeScan(400, 2000, 5, randomize=True),
                                                            ExplicitScan([500]),
                                                            CenterScan(500, 1000, 100, randomize=True),
                                                        ],
                                                        global_min=50, global_max=100000, global_step=1,
                                                        unit="us", scale=1, precision=5
                                                    ), tooltip="time between spin polarization pulses (in us)", group='SBC.base')

        # SBC - mode target sweep
        self.setattr_argument("sbc_mode_target",    NumberValue(default=0, precision=0, step=1, min=0, max=100), group='SBC.mode',
                              tooltip="The mode on which to sweep mode-specific SBC parameters.")
        self.setattr_argument("freq_sbc_scan_khz_list",  Scannable(
                                                        default=[
                                                            ExplicitScan([0.]),
                                                            CenterScan(0, 20, 5, randomize=True),
                                                            RangeScan(-10, 10, 20, randomize=True),
                                                        ],
                                                        global_min=-10000, global_max=10000, global_step=1.,
                                                        unit="kHz", scale=1, precision=3
                                                    ), group="SBC.mode", tooltip="Relative frequency for scanning the target 729nm mode.")
        self.setattr_argument("ampl_quench_scan_pct_list",   Scannable(
                                                            default=[
                                                                ExplicitScan([3.5]),
                                                                RangeScan(1, 5, 15, randomize=True),
                                                                CenterScan(3.5, 4., 0.2, randomize=True),
                                                            ],
                                                            global_min=0.01, global_max=50., global_step=1,
                                                            unit="%", scale=1, precision=3
                                                        ), group="SBC.mode", tooltip="Absolute amplitude to set for the target 729nm mode.")

        # beam parameters
        self.setattr_argument("enable_beam_sweep",  BooleanValue(default=False), group='SBC.beam', tooltip="Enable scanning of some beam's parameter.")
        self.setattr_argument("beam_sweep_target",  EnumerationValue(_VALID_BEAM_SCANS, default=_VALID_BEAM_SCANS[0]),
                              group='SBC.mode', tooltip="Beam parameter to scan.")
        self.setattr_argument("freq_beam_mhz_list",  Scannable(
                                                        default=[
                                                            ExplicitScan([115.]),
                                                            CenterScan(115, 6, 20, randomize=True),
                                                            RangeScan(112, 118, 20, randomize=True),
                                                        ],
                                                        global_min=60, global_max=400, global_step=1.,
                                                        unit="MHz", scale=1, precision=6
                                                    ), group='SBC.beam', tooltip="Absolute frequency for the beam parameter.")
        self.setattr_argument("ampl_beam_pct_list",   Scannable(
                                                            default=[
                                                                ExplicitScan([15]),
                                                                RangeScan(9, 20, 11, randomize=True),
                                                                CenterScan(15, 6., 0.2, randomize=True),
                                                            ],
                                                            global_min=0.01, global_max=50., global_step=1,
                                                            unit="%", scale=1, precision=3
                                                        ), group='SBC.beam', tooltip="Absolute amplitude for the beam parameter.")

        # get subsequences
        self.initialize_subsequence =       InitializeQubit(self)
        self.sidebandreadout_subsequence =  SidebandReadout(self, profile_dds=self.profile_729_readout)
        self.readout_subsequence =          Readout(self)
        self.rescue_subsequence =           RescueIon(self)

        # get relevant devices
        self.setattr_device('probe')
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')
        self.setattr_device('repump_qubit')
        self.setattr_device('qubit')

    def prepare_experiment(self):
        """
        Prepare & precompute experimental values.
        """
        # validate inputs
        self._prepare_argument_checks()

        '''
        CONVERT VALUES TO MACHINE UNITS
        '''
        # hardware values
        self.att_sbc_mu =   att_to_mu(self.att_sidebandcooling_continuous_db * dB)

        # scan - sbc base & sbc mode
        time_sbc_mu_list =  np.array([self.core.seconds_to_mu(time_us * us)
                                      for time_us in self.time_sbc_us_list])
        time_per_spinpol_mu_list =  np.array([self.core.seconds_to_mu(time_us * us)
                                              for time_us in self.time_per_spinpol_us_list])
        freq_sbc_scan_ftw_list =    np.array([self.qubit.frequency_to_ftw(freq_khz * kHz)
                                              for freq_khz in self.freq_sbc_scan_khz_list])
        ampl_quench_scan_asf_list = np.array([self.repump_qubit.amplitude_to_asf(ampl_pct / 100.)
                                              for ampl_pct in self.ampl_quench_scan_pct_list])
        # scan - beam parameters
        if self.enable_beam_sweep:
            # process target beam profile and device
            if self.beam_sweep_target == 'doppler_397':
                self.dds_beam = self.get_device('pump')
                self.beam_update_profile = 0
            elif self.beam_sweep_target == 'doppler_866':
                self.dds_beam = self.get_device('repump_cooling')
                self.beam_update_profile = 0
            elif self.beam_sweep_target == 'readout_397':
                self.dds_beam = self.get_device('pump')
                self.beam_update_profile = 1

            freq_beam_ftw_list = np.array([self.dds_beam.frequency_to_ftw(freq_mhz * MHz)
                                           for freq_mhz in self.freq_beam_mhz_list])
            ampl_beam_asf_list = np.array([self.dds_beam.amplitude_to_asf(ampl_pct / 100.)
                                           for ampl_pct in self.ampl_beam_pct_list])
        else:
            self.dds_beam = self.get_device('qubit')
            self.beam_update_profile = 6
            freq_beam_ftw_list = np.array([0x01], dtype=np.int32)
            ampl_beam_asf_list = np.array([0x01], dtype=np.int32)

        '''PREPARE SBC CONFIG/SCHEDULES'''
        self.num_modes = len(self.sideband_cooling_config_list)

        # create and fill SBC schedule
        self.sbc_config_base_list = np.zeros((self.num_modes, 3), dtype=np.int64)
        for i, params in enumerate(self.sideband_cooling_config_list.items()):
            self.sbc_config_base_list[i, 0] = self.qubit.frequency_to_ftw(params[0] * MHz)
            self.sbc_config_base_list[i, 1] = self.qubit.amplitude_to_asf(params[1][1] / 100.)
            self.sbc_config_base_list[i, 2] = np.int64(np.mean(time_sbc_mu_list) * params[1][0] / 100. / self.sideband_cycles_continuous)
        # create working copy for main loop
        self.sbc_config_update_list = np.copy(self.sbc_config_base_list)

        # create list of relative SBC times
        self.sbc_mode_time_frac_list = np.array([config_arr[0] / self.sideband_cycles_continuous
                                                 for config_arr in self.sideband_cooling_config_list.values()])

        '''
        CREATE EXPERIMENT CONFIG
        '''
        # create an array of values for the experiment to sweep
        self.config_experiment_list = np.array([
            list(vals)
            for vals in product(
                self.sidebandreadout_subsequence.freq_sideband_readout_ftw_list,
                time_sbc_mu_list, time_per_spinpol_mu_list,
                freq_sbc_scan_ftw_list, ampl_quench_scan_asf_list,
                freq_beam_ftw_list, ampl_beam_asf_list
            )
        ], dtype=np.int64)
        np.random.shuffle(self.config_experiment_list)

    def _prepare_argument_checks(self) -> TNone:
        """
        Check experiment arguments for validity.
        """
        pass
        # todo: actually do, and make it reasonable
        # # ensure a reasonable amount of modes (i.e. not too many)
        # if len(self.sideband_cooling_config_list) > 20:
        #     raise ValueError("Too many modes for SBC. Number of modes must be in [1, 20].")
        #
        # # ensure SBC config on all modes add up to 100%
        # mode_total_pct = np.sum([config_arr[0] for config_arr in self.sideband_cooling_config_list.values()])
        # if mode_total_pct > 100.:
        #     raise ValueError("Total sideband cooling mode percentages exceed 100%.")

        # todo: check SBC base config valid
        # todo: ensure SBC time enough for at least one mode
        # todo: check SBC time sweep isn't too extreme (incl. cycles)
        # todo: check time_per_spinpol isn't too extreme

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_experiment_list),
                8)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self) -> TNone:
        # record subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.sidebandreadout_subsequence.record_dma()
        self.readout_subsequence.record_dma()
        self.core.break_realtime()

        # update spinpol beam with SBC profile params
        self.probe.set_mu(self.probe.freq_spinpol_ftw, asf=self.probe.ampl_spinpol_asf, profile=self.profile_854_SBC,
                          phase_mode=ad9910.PHASE_MODE_CONTINUOUS)
        delay_mu(8000)

    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        for trial_num in range(self.repetitions):
            for config_vals in self.config_experiment_list:

                '''CONFIGURE & PREPARE'''
                # extract values from config list
                freq_readout_ftw =  np.int32(config_vals[0])
                time_sbc_mu =           config_vals[1]
                time_per_spinpol_mu =   config_vals[2]
                freq_sbc_scan_ftw =     config_vals[3]
                ampl_quench_scan_asf =  config_vals[4]
                freq_beam_ftw =     np.int32(config_vals[5])
                ampl_beam_asf =     np.int32(config_vals[6])

                # create & update SBC config w/target params
                self.sbc_config_update_list = self.sbc_config_base_list
                # update sbc config with new timing
                # self.sbc_config_update_list[:, 2] = np.int64(self.sbc_mode_time_frac_list * time_sbc_mu)
                for i in range(self.num_modes):
                    self.sbc_config_update_list[i, 2] = np.int64(time_sbc_mu * self.sbc_mode_time_frac_list[i])
                # update sbc config with target freq/quench ampl
                self.sbc_config_update_list[self.sbc_mode_target, 0] += freq_sbc_scan_ftw
                self.sbc_config_update_list[self.sbc_mode_target, 1] = ampl_quench_scan_asf

                # calculate penultimate mode & time
                time_counter_mu = np.int64(0)
                mode_counter = 0
                time_remainder_mu = np.int64(8)
                while time_counter_mu <= time_sbc_mu:
                    time_counter_mu += self.sbc_config_update_list[mode_counter % self.num_modes, 2]
                    mode_counter += 1
                if time_counter_mu > time_sbc_mu:
                    mode_counter -= 1
                    time_remainder_mu = np.int64(
                        time_sbc_mu -
                        (time_counter_mu - self.sbc_config_update_list[mode_counter % self.num_modes, 2])
                    )
                self.core.break_realtime()

                # print(self.sbc_config_update_list)
                # print(time_counter_mu)
                # print(mode_counter)
                # print(time_remainder_mu)
                # self.core.break_realtime()

                # prepare relevant beams
                self.qubit.set_mu(freq_readout_ftw, asf=self.sidebandreadout_subsequence.ampl_sideband_readout_asf,
                                  profile=self.profile_729_readout, phase_mode=ad9910.PHASE_MODE_CONTINUOUS)
                self.beam_update(freq_beam_ftw, ampl_beam_asf)
                # add heavy slack in case SBC config takes a while to schedule
                delay_mu(250000) # 250us


                '''INITIALIZE & SBC'''
                # doppler cool & initialize ion in S-1/2 state
                self.initialize_subsequence.run_dma()

                # synchronize spinpol w/SBC-ing
                time_start_mu = now_mu()
                self.spin_polarize(time_sbc_mu, time_per_spinpol_mu)
                at_mu(time_start_mu)
                self.sbc_schedule(mode_counter, time_remainder_mu)

                '''READOUT & SAVE RESULTS'''
                # sideband readout & detect fluorescence
                self.sidebandreadout_subsequence.run_dma()
                self.readout_subsequence.run_dma()

                # update dataset
                self.update_results(freq_readout_ftw, self.readout_subsequence.fetch_count(),
                                    time_sbc_mu, time_per_spinpol_mu, freq_sbc_scan_ftw, ampl_quench_scan_asf,
                                    freq_beam_ftw, ampl_beam_asf)
                self.core.break_realtime()

                # resuscitate ion
                self.rescue_subsequence.resuscitate()

            # rescue ion as needed
            self.rescue_subsequence.run(trial_num)

            # support graceful termination
            self.check_termination()
            self.core.break_realtime()


    @kernel(flags={"fast-math"})
    def spin_polarize(self, time_sbc_mu: TInt64=1, time_per_spinpol_mu: TInt64=10000) -> TNone:
        """
        Run spin polarization with variable timing.
        Arguments:
            num_spinpols: number of overall spin polarizations
            time_per_spinpol_mu: delay between successive spinpols
        """
        for _ in range(time_per_spinpol_mu, time_sbc_mu, time_per_spinpol_mu):
            delay_mu(time_per_spinpol_mu)

            self.probe.on()
            self.repump_cooling.on()
            delay_mu(self.initialize_subsequence.time_spinpol_mu)
            self.probe.off()
            self.repump_cooling.off()

    @kernel(flags={"fast-math"})
    def sbc_schedule(self, num_updates: TInt32=1, time_last_mu: TInt64=1000) -> TNone:
        """
        Run spin polarization with variable timing.
        Arguments:
            num_updates: number of modes to toggle through
            time_last_mu: amount of time to cool final mode
        """
        # prepare attenuation for SBC
        self.qubit.set_att_mu(self.att_sbc_mu)

        for i in range(num_updates):
            # update SBC beam parameters
            at_mu(now_mu() & ~7)
            with parallel:
                with sequential:
                    self.qubit.write64(ad9910._AD9910_REG_PROFILE0 + self.profile_729_SBC,
                                       self.qubit.ampl_qubit_asf << 16, # asf
                                       np.int32(self.sbc_config_update_list[i % self.num_modes, 0])) # ftw
                    delay_mu(np.int64(self.qubit.beam.sync_data.io_update_delay))
                    self.qubit.cpld.io_update.pulse_mu(8)

                with sequential:
                    self.repump_qubit.write64(ad9910._AD9910_REG_PROFILE0 + self.profile_854_SBC,
                                              np.int32(self.sbc_config_update_list[i % self.num_modes, 1] << 16), # asf
                                              self.repump_qubit.freq_repump_qubit_ftw) # ftw
                    delay_mu(np.int64(self.repump_qubit.beam.sync_data.io_update_delay))
                    self.repump_qubit.cpld.io_update.pulse_mu(8)
            # cool target mode
            delay_mu(self.sbc_config_update_list[i % self.num_modes, 2])

        # update SBC beam parameters - final mode
        at_mu(now_mu() & ~7)
        with parallel:
            with sequential:
                self.qubit.write64(ad9910._AD9910_REG_PROFILE0 + self.profile_729_SBC,
                                   self.qubit.ampl_qubit_asf << 16,  # asf
                                   np.int32(self.sbc_config_update_list[(num_updates + 1) % self.num_modes, 0]))  # ftw
                delay_mu(np.int64(self.qubit.beam.sync_data.io_update_delay))
                self.qubit.cpld.io_update.pulse_mu(8)

            with sequential:
                self.repump_qubit.write64(ad9910._AD9910_REG_PROFILE0 + self.profile_854_SBC,
                                          np.int32(self.sbc_config_update_list[(num_updates + 1) % self.num_modes, 1] << 16),  # asf
                                          self.repump_qubit.freq_repump_qubit_ftw)  # ftw
                delay_mu(np.int64(self.repump_qubit.beam.sync_data.io_update_delay))
                self.repump_qubit.cpld.io_update.pulse_mu(8)
        # cool target mode - final mode
        delay_mu(time_last_mu)

    @kernel(flags={"fast-math"})
    def beam_update(self, beam_freq_ftw: TInt32=0x01, beam_ampl_asf: TInt32=0x01) -> TNone:
        """
        Update target beam with parameters to enable scanning.
        Arguments:
            beam_freq_ftw: number of overall spin polarizations
            beam_ampl_asf: delay between successive spinpols
        """
        if self.enable_beam_sweep:
            self.dds_beam.set_mu(beam_freq_ftw, beam_ampl_asf, profile=self.beam_update_profile,
                                   phase_mode=ad9910.PHASE_MODE_CONTINUOUS)

