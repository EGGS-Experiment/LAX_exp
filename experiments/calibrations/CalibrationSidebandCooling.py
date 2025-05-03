import numpy as np
from artiq.experiment import *
from artiq.coredevice import ad9910

from LAX_exp.analysis import *
from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import (
    InitializeQubit, Readout, RescueIon,
    SidebandCoolContinuousRAM, SidebandCoolPulsed, SidebandReadout
)
# todo: sweep spinpol pct


class CalibrationSidebandCooling(LAXExperiment, Experiment):
    """
    Calibration: Sideband Cooling

    Sideband cooling but with scannable parameters for optimization.
    """
    name = 'Calibration Sideband Cooling'
    kernel_invariants = {
        # subsequences
        'initialize_subsequence', 'sidebandreadout_subsequence', 'readout_subsequence', 'rescue_subsequence',

        # configs
        'profile_729_readout', 'profile_729_SBC', 'config_experiment_list'
    }

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions",    NumberValue(default=40, precision=0, step=1, min=1, max=100000))
        self.setattr_argument("mode_target",    NumberValue(default=1, precision=0, step=1, min=1, max=10))

        # allocate profiles on 729nm for different subsequences
        self.profile_729_readout = 0
        self.profile_729_SBC = 1

        # base SBC
        self.setattr_argument("sideband_cooling_config_list", PYONValue({100.7555: [26., 5.], 100.455: [37., 5.], 100.315: [37., 5.]}),
                              tooltip="{freq_mode_mhz: [sbc_mode_pct_per_cycle, ampl_quench_mode_pct]}", group='SBC.base')
        self.setattr_argument("sideband_cycles_continuous", NumberValue(default=10, precision=0, step=1, min=1, max=10000),
                              tooltip="number of times to loop over the SBC configuration sequence", group='SBC.base')
        self.setattr_argument("time_per_spinpol_us",    NumberValue(default=600, precision=3, step=1, min=0.01, max=100000),
                              tooltip="time between spin polarization pulses (in us)", group='SBC.base')
        self.setattr_argument("att_sidebandcooling_continuous_db",  NumberValue(default=8, precision=1, step=0.5, min=8, max=31.5),
                              group='SBC_RAM.continuous')

        # SBC parameter scanning
        self.setattr_argument("time_sbc_us_list",   Scannable(
                                                        default=[
                                                            RangeScan(100, 2000, 100, randomize=True),
                                                            ExplicitScan([6.05]),
                                                            CenterScan(3.05, 5., 0.1, randomize=True),
                                                        ],
                                                        global_min=50, global_max=100000, global_step=1,
                                                        unit="us", scale=1, precision=5
                                                    ), group="SBC.sweep")
        self.setattr_argument("freq_sbc_scan_khz_list",  Scannable(
                                                        default=[
                                                            CenterScan(0, 20, 20, randomize=True),
                                                            ExplicitScan([0.]),
                                                            RangeScan(-10, 10, 20, randomize=True),
                                                        ],
                                                        global_min=-1000, global_max=1000, global_step=1.,
                                                        unit="kHz", scale=1, precision=3
                                                    ), group="SBC.sweep")
        self.setattr_argument("ampl_quench_pct_list",   Scannable(
                                                            default=[
                                                                RangeScan(1., 5., 15., randomize=True),
                                                                ExplicitScan([3.5]),
                                                                CenterScan(3.5, 4., 0.2, randomize=True),
                                                            ],
                                                            global_min=0.01, global_max=50., global_step=1,
                                                            unit="%", scale=1, precision=3
                                                        ), group="SBC.sweep")

        # Doppler cooling parameters
        self.setattr_argument("freq_397_doppler_mhz_list",  Scannable(
                                                        default=[
                                                            CenterScan(115, 6, 20, randomize=True),
                                                            ExplicitScan([115.]),
                                                            RangeScan(112, 118, 20, randomize=True),
                                                        ],
                                                        global_min=60, global_max=400, global_step=1.,
                                                        unit="MHz", scale=1, precision=6
                                                    ), group="Doppler.sweep")
        self.setattr_argument("ampl_397_doppler_pct_list",   Scannable(
                                                            default=[
                                                                RangeScan(9., 20., 22., randomize=True),
                                                                ExplicitScan([15]),
                                                                CenterScan(15, 6., 0.2, randomize=True),
                                                            ],
                                                            global_min=0.01, global_max=50., global_step=1,
                                                            unit="%", scale=1, precision=3
                                                        ), group="Doppler.sweep")

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
        # scan values
        freq_sbc_scan_ftw_list =    np.array([self.qubit.frequency_to_ftw(freq_khz * kHz)
                                              for freq_khz in self.freq_sbc_scan_khz_list])
        ampl_quench_asf_list =      np.array([self.repump_cooling.amplitude_to_asf(ampl_pct / 100.)
                                              for ampl_pct in self.ampl_quench_pct_list])
        time_sbc_mu_list =          np.array([self.core.seconds_to_mu(time_us * us)
                                              for time_us in self.time_sbc_us_list])
        freq_397_doppler_ftw_list = np.array([self.pump.frequency_to_ftw(freq_mhz * MHz)
                                              for freq_mhz in self.freq_397_doppler_mhz_list])
        ampl_397_doppler_asf_list = np.array([self.pump.amplitude_to_asf(ampl_pct / 100.)
                                              for ampl_pct in self.ampl_397_doppler_pct_list])

        # base SBC config
        self.att_sidebandcooling_mu =   att_to_mu(self.att_sidebandcooling_continuous_db * dB)

        # base sbc schedule
        mode_freqs_hz = np.array([
            freq_mhz * MHz
            for freq_mhz in self.sideband_cooling_config_list.keys()
        ])
        mode_time_pct = np.array([
            config_arr[0]
            for config_arr in self.sideband_cooling_config_list.values()
        ])
        mode_quench_ampls_frac = np.array([
            config_arr[1] / 100.
            for config_arr in self.sideband_cooling_config_list.values()
        ])


        sbc_schedule = np.zeros((len(self.sideband_cooling_config_list), 3))
        sbc_schedule[:, 0]



        '''
        CREATE EXPERIMENT CONFIG
        '''
        # create an array of values for the experiment to sweep
        # (i.e. heating time & readout FTW)
        self.config_experiment_list = np.stack(np.meshgrid(
            time_readout_mu_list,
            freq_397_readout_ftw_list, ampl_397_readout_asf_list,
            freq_866_readout_ftw_list, ampl_866_readout_asf_list
        ), -1).reshape(-1, 5)
        self.config_experiment_list = np.array(self.config_experiment_list, dtype=np.int64)
        np.random.shuffle(self.config_experiment_list)

    def _prepare_argument_checks(self) -> TNone:
        """
        Check experiment arguments for validity.
        """
        # ensure a reasonable amount of modes (i.e. not too many)
        if len(self.sideband_cooling_config_list) > 20:
            raise ValueError("Too many modes for SBC. Number of modes must be in [1, 20].")

        # ensure SBC config on all modes add up to 100%
        mode_total_pct = np.sum([config_arr[0] for config_arr in self.sideband_cooling_config_list.values()])
        if mode_total_pct > 100.:
            raise ValueError("Total sideband cooling mode percentages exceed 100%.")

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_experiment_list),
                6)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self) -> TNone:
        # record subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.sidebandreadout_subsequence.record_dma()
        self.readout_subsequence.record_dma()
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        self.core.break_realtime()

        for trial_num in range(self.repetitions):

            # scan over sideband readout frequencies
            for config_vals in self.config_experiment_list:
                self.core.break_realtime()

                '''CONFIGURE'''
                # extract values from config list
                time_readout_mu = config_vals[0]
                freq_397_readout_ftw = np.int32(config_vals[1])
                ampl_397_readout_asf = np.int32(config_vals[2])
                freq_866_readout_ftw = np.int32(config_vals[3])
                ampl_866_readout_asf = np.int32(config_vals[4])

                # set beam frequencies
                self.qubit.set_mu(freq_ftw, asf=self.sidebandreadout_subsequence.ampl_sideband_readout_asf,
                                  profile=self.profile_729_readout, phase_mode=ad9910.PHASE_MODE_CONTINUOUS)
                self.pump.set_mu()
                delay_mu(50000)

                # set cooling waveform
                self.pump.cooling(freq_doppler_ftw, asf=ampl_doppler_asf, profile=0, phase_mode=)

                # turn repumps on (and ensure spinpol off)
                self.probe.off()
                self.repump_qubit.on()
                self.repump_cooling.on()

                # doppler cooling
                self.pump.on()
                delay_mu(self.time_doppler_cooling_mu)
                self.pump.off()

                '''INITIALIZE & SBC'''
                # doppler cool & initialize ion in S-1/2 state
                self.initialize_subsequence.run_dma()

                # run sideband cooling sequence
                self.sidebandcool_subsequence.run_dma()

                '''READOUT & SAVE RESULTS'''
                # sideband readout & detect fluorescence
                self.sidebandreadout_subsequence.run_dma()
                self.readout_subsequence.run_dma()

                # update dataset
                self.update_results(freq_ftw, self.readout_subsequence.fetch_count())
                self.core.break_realtime()

                # resuscitate ion
                self.rescue_subsequence.resuscitate()

            # rescue ion as needed
            self.rescue_subsequence.run(trial_num)

            # support graceful termination
            self.check_termination()
            self.core.break_realtime()

