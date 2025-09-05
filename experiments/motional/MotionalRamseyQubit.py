from artiq.experiment import *
from artiq.coredevice import ad9910
from numpy import array, int32, int64

from LAX_exp.language import *
from LAX_exp.system.subsequences import (
    InitializeQubit, Readout, RescueIon, NoOperation, SidebandCoolContinuousRAM
)


class MotionalRamseyQubit(LAXExperiment, Experiment):
    """
    Experiment: Motional Ramsey Qubit

    Conduct motional Ramsey spectroscopy using the 729nm laser.

    """
    name = 'Motional Ramsey Qubit'
    kernel_invariants = {
        # hardware values
        'time_pulse_mu', 'ampl_ramsey_asf', 'att_ramsey_mu',
        'freq_carrier_ftw', 'ampl_carrier_asf', 'phase_carrier_rev_pow', 'time_carrier_mu',

        # subsequences
        'cooling_subsequence', 'initialize_subsequence', 'readout_subsequence', 'rescue_subsequence',

        # configs
        'profile_729_ramsey', 'profile_729_ramsey', 'config_experiment_list'
    }

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions",    NumberValue(default=40, precision=0, step=1, min=1, max=100000))
        self.setattr_argument("cooling_type",   EnumerationValue(["Doppler", "SBC"],
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
        self.sbc_subsequence =  SidebandCoolContinuousRAM(
            self, profile_729=self.profile_729_SBC, profile_854=3,
            ram_addr_start_729=0, ram_addr_start_854=0, num_samples=200
        )

        # ramsey pulse parameters (carrier)
        self.setattr_argument("freq_carrier_mhz",   NumberValue(default=101.0581, precision=6, step=1e-2, min=60, max=400, unit="MHz", scale=1.),
                              group="ramsey.carr")
        self.setattr_argument('ampl_carrier_pct',   NumberValue(default=50., precision=2, step=5., min=0.01, max=50., scale=1., unit="%"),
                              group="ramsey.carr")
        self.setattr_argument('time_carrier_us',    NumberValue(default=1.5, precision=3, step=0.5, min=0.01, max=10000., scale=1., unit="us"),
                              group="ramsey.carr")
        self.setattr_argument('phase_carrier_rev_turns', NumberValue(default=0.5, precision=4, step=0.1, min=-1., max=1., scale=1., unit="turns"),
                              group="ramsey.carr")

        # ramsey pulse parameters (RSB)
        self.setattr_argument('time_pulse_us', NumberValue(default=10., precision=2, step=0.5, min=1., max=10000., scale=1., unit="us"),
                              group="ramsey.RSB")
        self.setattr_argument('ampl_ramsey_pct', NumberValue(default=50., precision=2, step=5., min=0., max=50., scale=1., unit="%"),
                              group="ramsey.RSB")
        self.setattr_argument('att_ramsey_db', NumberValue(default=8., precision=1, step=0.5, min=8., max=31.5, scale=1., unit="dB"),
                              group="ramsey.RSB")


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
        if self.cooling_type == "Doppler":  self.cooling_subsequence =  self.doppler_subsequence
        elif self.cooling_type == "SBC":    self.cooling_subsequence =  self.sbc_subsequence

        '''
        CONVERT VALUES TO MACHINE UNITS
        '''
        # ramsey carrier parameters
        self.freq_carrier_ftw =  self.qubit.frequency_to_ftw(self.freq_carrier_mhz * MHz)
        self.ampl_carrier_asf =  self.qubit.amplitude_to_asf(self.ampl_carrier_pct / 100.)
        self.phase_carrier_rev_pow = self.qubit.phase_to_pow(self.phase_carrier_rev_turns)
        self.time_carrier_mu =    self.core.seconds_to_mu(self.time_carrier_us * us)

        # ramsey RSB parameters
        self.time_pulse_mu =    self.core.seconds_to_mu(self.time_pulse_us * us)
        self.ampl_ramsey_asf =  self.qubit.amplitude_to_asf(self.ampl_ramsey_pct / 100.)
        self.att_ramsey_mu =    self.qubit.cpld.att_to_mu(self.att_ramsey_db * dB)

        time_delay_mu_list =    array([self.core.seconds_to_mu(time_us * us)
                                       for time_us in list(self.time_delay_us_list)])
        freq_ramsey_ftw_list =  array([self.qubit.frequency_to_ftw(freq_mhz * MHz)
                                       for freq_mhz in list(self.freq_ramsey_mhz_list)])
        phase_ramsey_pow_list = array([self.qubit.turns_to_pow(phas_turn)
                                       for phas_turn in list(self.phase_ramsey_turns_list)])

        # linetrigger parameters
        time_linetrig_holdoff_mu_list = array([self.core.seconds_to_mu(time_ms * ms)
                                               for time_ms in self.time_linetrig_holdoff_ms_list])

        '''
        CREATE EXPERIMENT CONFIG
        '''
        self.config_experiment_list = create_experiment_config(
            time_delay_mu_list,
            freq_ramsey_ftw_list,
            phase_ramsey_pow_list,
            time_linetrig_holdoff_mu_list,
            shuffle_config=True,
            config_type=int64
        )

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_experiment_list),
                5)


    '''
    MAIN SEQUENCE
    '''
    @kernel(flags={"fast-math"})
    def initialize_experiment(self) -> TNone:
        # record subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.cooling_subsequence.record_dma()
        self.readout_subsequence.record_dma()

    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        for trial_num in range(self.repetitions):
            for config_vals in self.config_experiment_list:

                # tmp remove
                self.core.break_realtime()
                self.pump.rescue()
                self.repump_cooling.on()
                self.repump_qubit.on()
                self.pump.on()
                # tmp remove

                # extract values from config list
                time_delay_mu =     config_vals[0]
                freq_ramsey_ftw =   int32(config_vals[1])
                phase_ramsey_pow =  int32(config_vals[2])
                time_holdoff_mu =   config_vals[3]

                # wait for linetrigger
                if self.enable_linetrigger:
                    self.trigger_line.trigger(self.trigger_line.time_timeout_mu, time_holdoff_mu)

                # initialize ion in S-1/2 state & sideband cool
                self.initialize_subsequence.run_dma()
                self.cooling_subsequence.run_dma()

                # do ramsey sequence
                self.run_ramsey(time_delay_mu, freq_ramsey_ftw, phase_ramsey_pow)

                # read out, clean up, and update dataset
                self.readout_subsequence.run_dma()
                self.rescue_subsequence.resuscitate()
                self.update_results(
                    freq_ramsey_ftw,
                    self.readout_subsequence.fetch_count(),
                    time_delay_mu,
                    phase_ramsey_pow,
                    time_holdoff_mu
                )

            # rescue ion as needed & support graceful termination
            self.rescue_subsequence.run(trial_num)
            self.check_termination()

    @kernel(flags={"fast-math"})
    def run_ramsey(self, time_mu: TInt64, freq_ftw: TInt32, phase_pow: TInt32) -> TNone:
        """
        Run ramsey pulse sequence.
        :param time_mu: the ramsey delay time (in machine units).
        :param freq_ftw: frequency of the ramsey RSB pulse (in ftw).
        :param phase_pow: phase of the reverse RSB pulse (relative to the forward) in pow.
        """
        # prepare ramsey beam
        self.qubit.set_att_mu(self.att_ramsey_mu)
        self.qubit.set_profile(self.profile_729_ramsey)
        self.qubit.io_update()

        # get starting reference time
        t_start_mu = now_mu()

        # forward - sigma_x
        self.qubit.set_mu(self.freq_carrier_ftw, asf=self.ampl_carrier_asf, pow_=0x0,
                          profile=self.profile_729_ramsey,
                          phase_mode=ad9910.PHASE_MODE_TRACKING,
                          ref_time_mu=t_start_mu)
        self.qubit.on()
        delay_mu(self.time_carrier_mu)
        self.qubit.off()
        # forward - RSB
        self.qubit.set_mu(freq_ftw, asf=self.ampl_ramsey_asf, pow_=0x0,
                          profile=self.profile_729_ramsey,
                          phase_mode=ad9910.PHASE_MODE_TRACKING,
                          ref_time_mu=t_start_mu)
        self.qubit.on()
        delay_mu(self.time_pulse_mu)
        self.qubit.off()

        # ramsey delay
        delay_mu(time_mu)

        # reverse - RSB
        self.qubit.set_mu(freq_ftw, asf=self.ampl_ramsey_asf, pow_=phase_pow,
                          profile=self.profile_729_ramsey,
                          phase_mode=ad9910.PHASE_MODE_TRACKING,
                          ref_time_mu=t_start_mu)
        self.qubit.on()
        delay_mu(self.time_pulse_mu)
        self.qubit.off()
        # forward - sigma_x
        self.qubit.set_mu(self.freq_carrier_ftw, asf=self.ampl_carrier_asf, pow_=self.phase_carrier_rev_pow,
                          profile=self.profile_729_ramsey,
                          phase_mode=ad9910.PHASE_MODE_TRACKING,
                          ref_time_mu=t_start_mu)
        self.qubit.on()
        delay_mu(self.time_carrier_mu)
        self.qubit.off()

