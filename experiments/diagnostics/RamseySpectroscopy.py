import numpy as np
from artiq.experiment import *
from artiq.coredevice.ad9910 import PHASE_MODE_CONTINUOUS, PHASE_MODE_ABSOLUTE

from LAX_exp.language import *
from LAX_exp.system.subsequences import (
    InitializeQubit, Readout, RescueIon, NoOperation, SidebandCoolContinuousRAM
)
# todo: enable pi/2 pulse for motional stuff


class RamseySpectroscopy(LAXExperiment, Experiment):
    """
    Experiment: Ramsey Spectroscopy

    Measures ion fluorescence after conducting a Ramsey Spectroscopy sequence.
    """
    name = 'Ramsey Spectroscopy'
    kernel_invariants = {
        # hardware values
        'time_pulse_mu', 'ampl_ramsey_asf', 'att_ramsey_mu',
        'time_delay_mu_list', 'freq_ramsey_ftw_list', 'phase_ramsey_pow_list',

        # subsequences
        'cooling_subsequence', 'initialize_subsequence', 'readout_subsequence', 'rescue_subsequence',
        'time_linetrig_holdoff_mu_list',

        # configs
        'profile_729_ramsey', 'profile_729_ramsey', 'config_experiment_list'
    }

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions",    NumberValue(default=40, precision=0, step=1, min=1, max=100000))
        self.setattr_argument("cooling_type",   EnumerationValue(["Doppler", "SBC - Continuous"],
                                                                 default="Doppler"))

        # allocate profiles on 729nm for different subsequences
        self.profile_729_SBC =      4
        self.profile_729_ramsey =   7

        # linetrigger
        self.setattr_argument("enable_linetrigger", BooleanValue(default=False), group='linetrigger')
        self.setattr_argument("time_linetrig_holdoff_ms_list",   Scannable(
                                                                default=[
                                                                    ExplicitScan([0.1]),
                                                                    RangeScan(1, 3, 3, randomize=True),
                                                                ],
                                                                global_min=0.01, global_max=1000, global_step=1,
                                                                unit="ms", scale=1, precision=3
                                                            ), group='linetrigger')

        # ram-based continuous sideband cooling
        self.sidebandcool_continuous_subsequence =  SidebandCoolContinuousRAM(
            self, profile_729=self.profile_729_SBC, profile_854=3,
            ram_addr_start_729=0, ram_addr_start_854=0, num_samples=200
        )

        # ramsey pulse parameters
        self.setattr_argument('time_pulse_us', NumberValue(default=10., precision=2, step=0.5, min=1., max=10000., scale=1., unit="us"),
                              group="ramsey.pulse")
        self.setattr_argument('ampl_ramsey_pct', NumberValue(default=50., precision=2, step=5., min=0., max=50., scale=1., unit="%"),
                              group="ramsey.pulse")
        self.setattr_argument('att_ramsey_db', NumberValue(default=8., precision=1, step=0.5, min=8., max=31.5, scale=1., unit="dB"),
                              group="ramsey.pulse")

        # ramsey sweep parameters
        self.setattr_argument("time_delay_us_list", Scannable(
                                                        default=[
                                                            ExplicitScan([100.]),
                                                            RangeScan(20., 70., 10, randomize=True)
                                                        ],
                                                        global_min=3, global_max=10000000, global_step=100,
                                                        unit="us", scale=1, precision=2
                                                    ), group="ramsey.sweep")
        self.setattr_argument("freq_ramsey_mhz_list",   Scannable(
                                                            default=[
                                                                CenterScan(102.2152, 0.01, 0.0001, randomize=True),
                                                                ExplicitScan([104.335]),
                                                            ],
                                                            global_min=30, global_max=200, global_step=1,
                                                            unit="MHz", scale=1, precision=6
                                                        ), group="ramsey.sweep")
        self.setattr_argument("phase_ramsey_turns_list",    Scannable(
                                                            default=[
                                                                ExplicitScan([0.5]),
                                                                RangeScan(0., 1., 21, randomize=True),
                                                            ],
                                                            global_min=-1.1, global_max=1.1, global_step=0.1,
                                                            unit="turns", scale=1, precision=4
                                                        ), group="ramsey.sweep")
        # get devices
        self.setattr_device('qubit')
        self.setattr_device('trigger_line')

        # tmp remove
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')
        self.setattr_device('repump_qubit')
        # tmp remove

        # prepare sequences
        self.initialize_subsequence =   InitializeQubit(self)
        self.doppler_subsequence =      NoOperation(self)
        self.readout_subsequence =      Readout(self)
        self.rescue_subsequence =       RescueIon(self)

    def prepare_experiment(self):
        """
        Prepare & precompute experimental values.
        """
        '''
        CHOOSE COOLING SEQUENCE
        '''
        # choose correct cooling subsequence
        if self.cooling_type == "Doppler":              self.cooling_subsequence =  self.doppler_subsequence
        elif self.cooling_type == "SBC - Continuous":   self.cooling_subsequence =  self.sidebandcool_continuous_subsequence

        '''
        CONVERT VALUES TO MACHINE UNITS
        '''
        # ramsey pulse parameters
        self.time_pulse_mu =    self.core.seconds_to_mu(self.time_pulse_us * us)
        self.ampl_ramsey_asf =  self.qubit.amplitude_to_asf(self.ampl_ramsey_pct / 100.)
        self.att_ramsey_mu =    self.qubit.cpld.att_to_mu(self.att_ramsey_db * dB)

        self.time_delay_mu_list =       np.array([self.core.seconds_to_mu(time_us * us)
                                                  for time_us in list(self.time_delay_us_list)])
        self.freq_ramsey_ftw_list =     np.array([self.qubit.frequency_to_ftw(freq_mhz * MHz)
                                                  for freq_mhz in list(self.freq_ramsey_mhz_list)])
        self.phase_ramsey_pow_list =    np.array([self.qubit.turns_to_pow(phas_turn)
                                                  for phas_turn in list(self.phase_ramsey_turns_list)])

        # linetrigger parameters
        self.time_linetrig_holdoff_mu_list =    np.array([self.core.seconds_to_mu(time_ms * ms)
                                                          for time_ms in self.time_linetrig_holdoff_ms_list])

        '''
        CREATE EXPERIMENT CONFIG
        '''
        # create an array of values for the experiment to sweep
        # (i.e. heating time & readout FTW)
        self.config_experiment_list = np.stack(np.meshgrid(self.time_delay_mu_list,
                                                       self.freq_ramsey_ftw_list,
                                                       self.phase_ramsey_pow_list,
                                                       self.time_linetrig_holdoff_mu_list),
                                                 -1).reshape(-1, 4)
        self.config_experiment_list = np.array(self.config_experiment_list, dtype=np.int64)
        np.random.shuffle(self.config_experiment_list)

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_experiment_list),
                5)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def initialize_experiment(self) -> TNone:
        # record subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.cooling_subsequence.record_dma()
        self.readout_subsequence.record_dma()

    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        for trial_num in range(self.repetitions):

            # sweep exp config
            for config_vals in self.config_experiment_list:

                # tmp remove
                self.pump.rescue()
                self.repump_cooling.on()
                self.repump_qubit.on()
                self.pump.on()
                # tmp remove

                # extract values from config list
                time_delay_mu =     config_vals[0]
                freq_ramsey_ftw =   np.int32(config_vals[1])
                phase_ramsey_pow =  np.int32(config_vals[2])
                time_holdoff_mu =   config_vals[3]
                self.core.break_realtime()

                # wait for linetrigger
                if self.enable_linetrigger:
                    self.trigger_line.trigger(self.trigger_line.time_timeout_mu, time_holdoff_mu)

                # initialize ion in S-1/2 state & sideband cool
                self.initialize_subsequence.run_dma()
                self.cooling_subsequence.run_dma()

                # do ramsey sequence
                self.run_ramsey(time_delay_mu, freq_ramsey_ftw, phase_ramsey_pow)

                # do readout
                self.readout_subsequence.run_dma()

                # update dataset
                self.update_results(
                    freq_ramsey_ftw,
                    self.readout_subsequence.fetch_count(),
                    time_delay_mu,
                    phase_ramsey_pow,
                    time_holdoff_mu
                )
                self.core.break_realtime()

                # resuscitate ion
                self.rescue_subsequence.resuscitate()

            # rescue ion as needed
            self.rescue_subsequence.run(trial_num)

            # support graceful termination
            self.check_termination()
            self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def run_ramsey(self, time_mu: TInt64, freq_ftw: TInt32, phase_pow: TInt32) -> TNone:
        """
        Run ramsey pulse sequence.
        Arguments:
            time_mu:    the ramsey delay time (in machine units).
            freq_ftw:   the frequency of the ramsey pulse (in ftw).
            phase_pow:  the phase of the reverse pulse (relative to the forward) in pow.
        """
        # prepare ramsey beam
        self.qubit.set_att_mu(self.att_ramsey_mu)
        self.qubit.set_profile(self.profile_729_ramsey)
        self.qubit.io_update()

        # first pulse
        self.qubit.set_mu(freq_ftw, asf=self.ampl_ramsey_asf, pow_=0x0,
                          profile=self.profile_729_ramsey, phase_mode=PHASE_MODE_ABSOLUTE)
        self.qubit.on()
        delay_mu(self.time_pulse_mu)
        self.qubit.off()

        # ramsey delay
        delay_mu(time_mu)

        # second pulse
        self.qubit.set_mu(freq_ftw, asf=self.ampl_ramsey_asf, pow_=phase_pow,
                          profile=self.profile_729_ramsey, phase_mode=PHASE_MODE_CONTINUOUS)
        self.qubit.on()
        delay_mu(self.time_pulse_mu)
        self.qubit.off()

