from LAX_exp.base.base_experiment import LAXExperiment
from artiq.experiment import *
from LAX_exp.language import *
from LAX_exp.system.subsequences import InitializeQubit, Readout


class OpticalMassSpectrometry(LAXExperiment):

    def build_experiment(self):

        self.profile_CPT = 4

        # set pump parameters
        self.setattr_argument('pump_freq_mhz',
                            NumberValue(default=120., min=100, max=140,
                                        step=0.001, precision=5, unit='MHz'),
                              group ='pump_parameters')
        self.setattr_argument('pump_ampl_pct',
                              NumberValue(default=50, min=0, max=50,
                                          step=0.01, precision=2, unit='%' ),
                              group ='pump_parameters')
        self.setattr_argument('pump_att_dB',
                              NumberValue(default=14., min=14., max=31.5,
                                          step=0.5, precision=1, unit='dB'),
                              group ='pump_parameters')

        # set repump parameters
        self.setattr_argument('repump_freq_mhz',
                              NumberValue(default=120., min=100, max=140,
                                          step=0.001, precision=5, unit='MHz'),
                              group='repump_parameters')

        self.setattr_argument('repump_ampl_pct',
                              NumberValue(default=50, min=0, max=50,
                                          step=0.01, precision=2, unit='%'),
                              group='repump_parameters')

        self.setattr_argument('repump_att_dB',
                              NumberValue(default=14., min=14., max=31.5,
                                          step=0.5, precision=1, unit='dB'),
                              group='repump_parameters')


        # set number of counts to get
        self.setattr_argument('num_pmt_counts',
                             NumberValue(default=1e6, min=1, max=1e8,
                                         precision=0, step=1),
                              group ='pmt_parameters')
        self.setattr_argument('max_wait_time_per_count_us',
                              NumberValue(default=10., min=1., max=1000.,
                                          precision=2, step=1, unit='us'),
                              group ='pmt_parameters')


        # get devices
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')
        self.setattr_device('pmt')

        # get subsequences
        self.initialize_qubit =         InitializeQubit(self)
        self.rescue_subsequence =       RescueIon(self)


    def prepare_experiment(self):

        self.timestamps = np.zeros(self_num_pmt_counts)

        # convert pump units
        self.pump_freq_ftw = mhz_to_ftw(self.pump_freq_mhz)
        self.pump_ampl_asf = pct_to_asf(self.pump_ampl_asf)
        self.pump_att_mu = att_to_mu(self.pump_att_dB)

        # convert repump units
        self.repump_freq_ftw = mhz_to_ftw(self.repump_freq_mhz)
        self.repump_ampl_asf = pct_to_asf(self.repump_ampl_asf)
        self.repump_att_mu = att_to_mu(self.repump_att_dB)

        # convert max wait time
        self.max_wait_time_mu = us_to_mu(self.max_wait_time_per_count_us)

    def results_shape(self):
        return (self.num_pmt_counts)

    @kernel(flags={'fast-math'})
    def initialize_experiment(self) -> TNone:
        self.core.break_realtime()

        self.initialize_qubit.record_dma()

    @kernel(flags={'fast-math'})
    def run_main(self) -> None:

        # get slack
        self.core.break_realtime()
        delay_mu(100000)

        # Initialize Calcium to S1/2 mj=-1/2
        self.initialize_qubit.run_dma()

        # set pump (397nm cooling) parameters
        self.pump.set_mu(self.pump_freq_ftw, asf = self.pump_ampl_asf, profile = self.profile_CPT)
        self.pump.set_att_mu(self.pump_att_mu)
        self.pump.set_profile(self.profile_CPT)

        # set rempump (866nm) parameters
        self.repump_cooling.set_mu(self.repump_freq_ftw, asf=self.repump_ampl_asf, profile = self.profile_CPT)
        self.repump_cooling.set_att_mu(self.repump_att_mu)
        self.repump_cooling.set_profile(self.profile_CPT)

        # turn beams on
        self.pump.on()
        self.repump_cooling.on()

        # todo: make sure this works
        # fill timestamp counts
        self.pmt.timestamp_counts(self.timestamps, self.max_wait_time_mu)

        # turn beams off
        self.pump.off()
        self.repump_cooling.off()

        # store results
        self.set_dataset('results', self.timestamps)


    def analyze(self):
        pass



