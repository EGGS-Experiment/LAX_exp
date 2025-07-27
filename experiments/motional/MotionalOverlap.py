from numpy import int32
from artiq.experiment import *
from artiq.coredevice.ad9910 import PHASE_MODE_CONTINUOUS

from LAX_exp.language import *
from LAX_exp.system.subsequences import (
    InitializeQubit, Readout, RescueIon,
    SidebandCoolContinuousRAM, AgilePulseGenerator, OverlapReadout
)


class MotionalOverlap(LAXExperiment, Experiment):
    """
    Experiment: Motional Overlap
    Conduct metrology using the motional fock overlap technique from F. Wolf (P.O. Schmidt group).
    https://www.nature.com/articles/s41467-019-10576-4
    """
    name = 'Motional Overlap'
    kernel_invariants = {
        # hardware parameters
        'ampl_rabi_asf', 'att_rabi_mu',

        # subsequences
        'initialize_subsequence', 'sidebandcool_subsequence', 'readout_subsequence', 'rescue_subsequence',

        # configs
        'config_experiment_list', 'profile_729_rabi', 'profile_729_SBC', 'profile_729_fock', 'profile_729_RAP',
        'profile_729_overlap'
    }

    def build_experiment(self):
        # core arguments
        self.setattr_argument("repetitions", NumberValue(default=40, precision=0, step=1, min=1, max=10000))
    
        # allocate relevant beam profiles
        self.profile_729_rabi =     0
        self.profile_729_SBC =      1
        self.profile_729_fock =     3
        self.profile_729_RAP =      4
        self.profile_729_overlap =  5

        # prepare argument-based subsequences
        self.sidebandcool_subsequence = SidebandCoolContinuousRAM(
            self, profile_729=self.profile_729_SBC, profile_854=3,
            ram_addr_start_729=0, ram_addr_start_854=0, num_samples=100
        )

        # fock state generation: arguments
        self.setattr_argument("enable_fock_prep", BooleanValue(default=False), group='fock_prep')
        self.setattr_argument("fock_pulse_config", PYONValue([[100.7528, 26., 3.], [100.7528, 26., 3.], [100.7528, 26., 3.]]),
                              group='fock_prep', tooltip="List of [freq_mhz, ampl_pct, time_us].")
        self.setattr_argument("att_fock_db", NumberValue(default=31.5, precision=1, step=0.5, min=8, max=31.5, scale=1., unit='dB'),
                              group='fock_prep')

        # spectroscopy: arguments
        self.setattr_argument("enable_rabi", BooleanValue(default=False), group='Spectroscopy')
        self.setattr_argument("time_rabi_us_list",  Scannable(
                                                        default=[
                                                            RangeScan(1, 150, 100, randomize=True),
                                                            ExplicitScan([6.05]),
                                                            CenterScan(3.05, 5., 0.1, randomize=True),
                                                        ],
                                                        global_min=1, global_max=100000, global_step=1,
                                                        unit="us", scale=1, precision=5
                                                    ), group="Spectroscopy")
        self.setattr_argument("freq_rabi_mhz_list", Scannable(
                                                        default=[
                                                            CenterScan(100.3172, 0.02, 0.00025, randomize=True),
                                                            ExplicitScan([100.7044]),
                                                        ],
                                                        global_min=30, global_max=200, global_step=1,
                                                        unit="MHz", scale=1, precision=5
                                                    ), group="Spectroscopy")
        self.setattr_argument("ampl_rabi_pct",  NumberValue(default=50, precision=3, step=5, min=1, max=50, scale=1., unit="%"),
                              group="Spectroscopy")
        self.setattr_argument("att_rabi_db",    NumberValue(default=8, precision=1, step=0.5, min=8, max=31.5, scale=1., unit="dB"),
                              group="Spectroscopy")

        # instantiate motional overlap subsequences
        self.overlap_subsequence = OverlapReadout(
            self, ram_profile=self.profile_729_rap, profile_shelve=self.profile_729_overlap,
            ram_addr_start=0, num_samples=500,
            ampl_max_pct=self.ampl_rap_pct, pulse_shape="blackman"
        )
        self.fock_subsequence = AgilePulseGenerator(
            self, profile_agile=self.profile_729_fock, att_pulse_db=self.att_fock_db,
            pulse_config=self.fock_pulse_config
        )

        # prepare other sequences
        self.initialize_subsequence =   InitializeQubit(self)
        self.readout_subsequence =      Readout(self)
        self.rescue_subsequence =       RescueIon(self)

        # get devices
        self.setattr_device('qubit')

    def prepare_experiment(self):
        """
        Prepare values for speedy evaluation.
        """
        # prepare arguments: rabi/spectroscopy
        freq_rabi_ftw_list = [self.qubit.frequency_to_ftw(freq_mhz * MHz) for freq_mhz in self.freq_rabi_mhz_list]
        time_rabi_mu_list = [self.core.seconds_to_mu(time_us * us) for time_us in self.time_rabi_us_list]
        self.ampl_rabi_asf = self.qubit.amplitude_to_asf(self.ampl_rabi_pct / 100.)
        self.att_rabi_mu = att_to_mu(self.att_rabi_db * dB)

        # create an array of values for the experiment to sweep
        self.config_experiment_list = create_experiment_config(
            freq_rabi_ftw_list, time_rabi_mu_list,
            shuffle_config=True
        )

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_experiment_list),
                3)


    '''
    MAIN SEQUENCE
    '''
    @kernel(flags={"fast-math"})
    def initialize_experiment(self) -> TNone:
        # store qubit CPLD attenuations
        self.qubit.cpld.get_att_mu()
        self.core.break_realtime()

        # record subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.sidebandcool_subsequence.record_dma()
        self.fock_subsequence.record_dma()
        self.core.break_realtime()

        self.overlap_subsequence.record_dma()
        self.readout_subsequence.record_dma()
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        # MAIN LOOP
        for trial_num in range(self.repetitions):

            for config_vals in self.config_experiment_list:
                
                '''PREPARE & CONFIGURE'''
                # extract values from config list
                freq_rabiflop_ftw = int32(config_vals[0])
                time_rabiflop_mu =  config_vals[1]
                self.core.break_realtime()


                '''STATE PREPARATION'''
                # initialize ion in S-1/2 state & sideband cool
                self.initialize_subsequence.run_dma()
                self.cooling_subsequence.run_dma()
                if self.enable_fock_prep:
                    self.fock_subsequence.run_dma()


                '''SPECTROSCOPY'''
                # run rabi flopping
                if self.enable_rabi:
                    self.qubit.set_att_mu(self.att_rabi_mu)
                    self.qubit.set_profile(self.profile_729_rabi)
                    self.qubit.set_mu(freq_rabiflop_ftw, asf=self.ampl_rabi_asf,
                                      profile=self.profile_729_rabi,
                                      phase_mode=PHASE_MODE_CONTINUOUS)
                    self.qubit.on()
                    delay_mu(time_rabiflop_mu)
                    self.qubit.off()


                '''READOUT'''
                # run motional overlap subsequence,
                # then read out via state-selective fluorescence
                if self.enable_overlap:
                    self.overlap_subsequence.run_dma()
                self.readout_subsequence.run_dma()

                # update dataset & clean up
                self.update_results(freq_rabiflop_ftw, self.readout_subsequence.fetch_count(), time_rabiflop_mu)
                self.rescue_subsequence.resuscitate()

            # rescue ion as needed
            self.rescue_subsequence.run(trial_num)

            # support graceful termination
            self.check_termination()
            self.core.break_realtime()


    '''
    ANALYSIS
    '''
    def analyze_experiment(self):
        # todo
        pass
