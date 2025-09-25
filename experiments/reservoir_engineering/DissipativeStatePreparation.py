from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.language import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import InitializeQubit, Readout, RabiflopReadout, SpinPolarizationRE
from numpy import tanh, sqrt
from artiq.coredevice import ad9910


class DissipativeStatePreparation(LAXExperiment, Experiment):
    """
    Experiment: Dissipative State Preparation

    From a Doppler-cooled ion apply multi-chromatic beams to reservoir engineer a specified motional state
    """
    name = 'Dissipative State Preparation'
    kernel_invariants = {

        'repetitions',

        # 729 multi-chromatic
        'dsp_phase_bsb_turns', 'dsp_freq_carrier_MHz', 'dsp_freq_secular_kHz', 'dsp_chromatic_att_dB',
        'dsp_ampl_rsb_pct', 'dsp_freq_rsb_MHz', 'dsp_freq_bsb_MHz', 'dsp_ampl_bsb_pct', 'dsp_squeeze_r',

        # 729 multi-chromatic machine unit params
        'dsp_freq_carrier_ftw', 'dsp_freq_secular_ftw', 'dsp_freq_rsb_ftw', 'dsp_freq_bsb_ftw', 'dsp_chromatic_att_mu',
        'dsp_phase_bsb_pow', 'dsp_phase_rsb_pow', 'dsp_ampl_bsb_asf', 'dsp_ampl_rsb_asf',

        # 729 multi-chromatic default values
        'singlepass0_default_amp_asf', 'singlepass1_default_amp_asf','singlepass0_default_att_mu',
        'singlepass1_default_att_mu', 'singlepass0_default_freq_MHz', 'singlepass1_default_freq_MHz',

        # timing
        'dsp_time_mu', 'rabiflop_readout_times_us_list', 'rabiflop_readout_times_mu_list',

        # 729 double pass
        'dsp_douplepass_carrier_att_mu', 'dsp_douplepass_carrier_asf',

        # subsequences
        'readout_subsequence', 'initialize_subsequence', 'rabiflop_readout_subsequence',
        'spin_polarization_re_subsequence',

        # devices
        'singlepass0', 'singlepass1',
    }

    def build_experiment(self):

        # get arguments
        self.setattr_argument('repetitions',
                              NumberValue(default=50., precision=0, step=1, min=0.,
                                          max=1e6), group=self.name,
                              tooltip='parameter determining amount of squeeze after reservoir engineering')

        # dissipative state preparation parameters
        self.setattr_argument('dsp_squeeze_r',
                              NumberValue(default=1., precision=3, step=0.01, min=0.,
                                          max=5.), group=self.name,
                              tooltip='parameter determining amount of squeeze after reservoir engineering')

        self.setattr_argument('dsp_time_us', NumberValue(default=20000, precision=3, step=10, min=1,
                                                         max=1000000), group=self.name,
                              tooltip='How long to apply bichromatic beams for reservoir engineering')

        # AOM frequenices
        self.setattr_argument('dsp_freq_carrier_MHz', NumberValue(default=101.0928, precision=3, step=1e-3, max=120.,
                                                                  min=80., unit="MHz"), group=self.name,
                              tooltip='frequency of the urukul channel used to drive the red sideband')

        self.setattr_argument('dsp_freq_secular_kHz', NumberValue(default=701.6, precision=3, step=1e-3, max=120.,
                                                                  min=80., unit="kHz"), group=self.name,
                              tooltip='frequency of the urukul channel used to drive the blue sideband')

        # AOM attenuations
        self.setattr_argument('dsp_chromatic_att_dB', NumberValue(default=13., precision=1, step=1e-3, max=31.5,
                                                                  min=13., unit="dB"), group=self.name,
                              tooltip='attenuation to apply to both of the urukuls channels (i.e the urukul channels used to '
                                      'drive the red and blue sideband are set to the same attenuation')

        # AOM phases
        self.setattr_argument('dsp_phase_bsb_turns', NumberValue(default=0., precision=2, step=0.01, max=1.,
                                                                 min=0.), group=self.name,
                              tooltip='phase offset between bsb and rsb tones')

        # AOM amplitudes
        self.setattr_argument('dsp_ampl_rsb_pct',
                              NumberValue(default=40., step=5., min=0.01, max=50., precision=1,  unit="%", scale=1),
                              group=self.name,
                              tooltip='amplitude of drive of urukul channel driving red sideband. Blue sideband amplitude will be '
                                      'calculated from this value')

        # rabi flopping readout timing
        self.setattr_argument("rabiflop_readout_times_us_list", Scannable(
            default=[
                RangeScan(1, 500, 250, randomize=True),
                ExplicitScan([6.05]),
                CenterScan(3.05, 5., 0.1, randomize=True),
            ],
            global_min=1, global_max=100000, global_step=1,
            unit="us", scale=1, precision=5
        ), group=self.name)

        # use as placeholder
        self.dsp_ampl_bsb = 0.

        self.setattr_device('qubit')
        self.setattr_device('repump_qubit')
        self.setattr_device('repump_cooling')
        self.setattr_device('probe')
        self.singlepass0 = self.get_device('urukul0_ch1')
        self.singlepass1 = self.get_device('urukul0_ch2')

        self.readout_subsequence = Readout(self)
        self.initialize_subsequence = InitializeQubit(self)
        self.rabiflop_readout_subsequence = RabiflopReadout(self)
        self.spin_polarization_re_subsequence = SpinPolarizationRE(self)

    def prepare_experiment(self):
        #
        self.get_singlepass_default_values()

        # get red and blue sideband values
        self.dsp_freq_rsb_MHz = (self.singlepass0_default_freq_MHz * MHz - self.dsp_freq_secular_kHz * kHz) / MHz
        self.dsp_freq_bsb_MHz = (self.singlepass1_default_freq_MHz * MHz + self.dsp_freq_secular_kHz * kHz) / MHz

        # convert parameters to machine units
        self.dsp_freq_carrier_ftw = self.qubit.frequency_to_ftw(self.dsp_freq_carrier_MHz * MHz)
        self.dsp_freq_secular_ftw = self.singlepass1.frequency_to_ftw(self.dsp_freq_secular_kHz * kHz)
        self.dsp_freq_rsb_ftw = self.singlepass0.frequency_to_ftw(self.dsp_freq_rsb_MHz * MHz)
        self.dsp_freq_bsb_ftw = self.singlepass1.frequency_to_ftw(self.dsp_freq_bsb_MHz * MHz)

        self.readout_freq_bsb_ftw = self.qubit.frequency_to_ftw(101.423*MHz)

        self.dsp_time_mu = self.core.seconds_to_mu(self.dsp_time_us * us)

        self.dsp_chromatic_att_mu = att_to_mu(self.dsp_chromatic_att_dB * dB)
        # todo: see if can use a parameter???
        self.dsp_douplepass_carrier_att_mu = att_to_mu(8. * dB)

        # todo: see if can use a parameter???
        self.dsp_phase_rsb_pow = self.singlepass0.turns_to_pow(0.)
        self.dsp_phase_bsb_pow = self.singlepass1.turns_to_pow(self.dsp_phase_bsb_turns)

        # todo: implement interpolation functions to scale frequency effects coupling to fiber
        # todo: implement interpolation functions to scale difference in power from channels

        # find bsb amplitude from rsb amplitude
        self.dsp_ampl_bsb_pct = sqrt(tanh(self.dsp_squeeze_r)) * self.dsp_ampl_rsb_pct

        self.dsp_ampl_rsb_asf = pct_to_asf(self.dsp_ampl_rsb_pct)
        self.dsp_ampl_bsb_asf = pct_to_asf(self.dsp_ampl_bsb_pct)
        # todo: see if can use a parameter???
        self.dsp_douplepass_carrier_asf = pct_to_asf(50.)

        # get time values
        self.rabiflop_readout_times_mu_list = [self.core.seconds_to_mu(rabiflop_readout_time_us*us)
                                               for rabiflop_readout_time_us in self.rabiflop_readout_times_us_list]


    @property
    def results_shape(self):
        return (self.repetitions * len(self.rabiflop_readout_times_mu_list),
                4)

    @kernel(flags={"fast-math"})
    def initialize_experiment(self) -> TNone:

        self.core.break_realtime()

        self.qubit.set_profile(0)
        self.qubit.set_att_mu(self.dsp_douplepass_carrier_att_mu)
        self.qubit.set_mu(ftw=self.dsp_freq_carrier_ftw, asf=self.dsp_douplepass_carrier_asf, profile =0)
        self.qubit.set_att_mu(self.dsp_douplepass_carrier_att_mu)
        self.qubit.off()

        self.singlepass0.set_att_mu(self.dsp_chromatic_att_mu)
        self.singlepass1.set_att_mu(self.dsp_chromatic_att_mu)

        self.initialize_subsequence.record_dma()
        self.readout_subsequence.record_dma()
        self.spin_polarization_re_subsequence.record_dma()
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        self.core.break_realtime()

        for trial_num in range(self.repetitions):

            for rabiflop_readout_time_mu in self.rabiflop_readout_times_mu_list:
                self.core.break_realtime()

                # Doppler cool and optically pump 40Ca+ to |S_1/2, mj=-1/2>
                self.initialize_subsequence.run_dma()

                # get a reference time so all subsequent pulses are phase coherent
                ref_time_mu = now_mu() & ~0x7
                # # configure the single pass values
                self.configure_singlepass(ref_time_mu)

                # # execute the dissipative state preparation
                self.pulse_dsp()

                # # ensure everything is in ground state by optical pumping (turn on 397 sigma+, 866, 854)
                self.spin_polarization_re_subsequence.run_dma()

                # # rabi flop
                self.set_default_singlepass_values()
                self.singlepass0.sw.on()
                self.singlepass1.sw.on()
                self.qubit.set_profile(0)
                self.qubit.set_att_mu(self.dsp_douplepass_carrier_att_mu)
                self.qubit.set_mu(ftw=self.readout_freq_bsb_ftw, asf=self.dsp_douplepass_carrier_asf, profile=0)
                self.qubit.on()
                delay_mu(rabiflop_readout_time_mu)
                self.qubit.off()
                self.qubit.set_mu(ftw=self.dsp_freq_carrier_ftw, asf=self.dsp_douplepass_carrier_asf, profile=0)

                # # scatter photons off |S_1/2> -> |P_1/2>  transition
                self.readout_subsequence.run_dma()
                counts = self.readout_subsequence.fetch_count()

                # self.check_termination()  # check termination b/c we haven't in a while
                self.update_results(rabiflop_readout_time_mu,
                                    counts,
                                    0,
                                    0)

                # self.core.break_realtime()
                # self.check_termination()
    @kernel(flags={"fast-math"})
    def cleanup_experiment(self) -> TNone:

        # reset singlepass values
        self.singlepass0.set_mu(ftw=self.singlepass0_default_freq_ftw,
                                asf=self.singlepass0_default_amp_asf,
                                )
        self.singlepass0.set_att_mu(self.singlepass0_default_att_mu)

        self.singlepass1.set_mu(ftw=self.singlepass1_default_freq_ftw,
                                asf=self.singlepass1_default_amp_asf,
                                )
        self.singlepass1.set_att_mu(self.singlepass1_default_att_mu)

        # turn off lasers
        self.qubit.off()
        self.probe.off()
        self.repump_cooling.off()
        self.repump_qubit.off()

    @kernel(flags={"fast-math"})
    def configure_singlepass(self, ref_time_mu: TInt64):
        self.singlepass0.set_cfr1()
        self.singlepass1.set_cfr1()
        self.qubit.cpld.io_update.pulse_mu(8)
        self.singlepass0.set_mu(ftw=self.dsp_freq_rsb_ftw,
                                pow_=self.dsp_phase_rsb_pow,
                                asf=self.dsp_ampl_rsb_asf,
                                phase_mode=ad9910.PHASE_MODE_TRACKING,
                                ref_time_mu=ref_time_mu
                                )

        self.singlepass1.set_mu(ftw=self.dsp_freq_bsb_ftw,
                                pow_=self.dsp_phase_bsb_pow,
                                asf=self.dsp_ampl_bsb_asf,
                                phase_mode=ad9910.PHASE_MODE_TRACKING,
                                ref_time_mu=ref_time_mu
                                )

    @kernel(flags={"fast-math"})
    def set_default_singlepass_values(self):
        self.singlepass0.set_mu(ftw=self.singlepass0_default_freq_ftw ,
                                asf=self.singlepass0_default_amp_asf
                                )

        self.singlepass1.set_mu(ftw=self.singlepass0_default_freq_ftw ,
                                asf= self.singlepass1_default_amp_asf,
                                )

        self.singlepass0.set_att_mu(self.singlepass0_default_att_mu)
        self.singlepass1.set_att_mu(self.singlepass1_default_att_mu)


    @kernel(flags={"fast-math"})
    def pulse_dsp(self) -> TNone:
        # ensure doublepass is on
        self.qubit.on()

        # turn on dissipation lasers
        self.repump_cooling.on()
        self.repump_qubit.on()
        self.probe.on()

        # turn on single pass
        self.singlepass0.sw.on()
        self.singlepass1.sw.on()

        delay_mu(self.dsp_time_mu)

        # turn off singlepass
        self.singlepass0.sw.off()
        self.singlepass1.sw.off()

        # turn off all the lasers
        self.qubit.off()
        self.repump_cooling.off()
        self.repump_qubit.off()
        self.probe.off()

    def get_singlepass_default_values(self) -> TNone:
            # get default values for the 729 single pass
            self.singlepass0_default_amp_asf = self.get_parameter('ampl_729_singlepass0_pct',
                                                                  group='beams.ampl_pct', override=True,
                                                                  conversion_function=pct_to_asf)

            self.singlepass1_default_amp_asf = self.get_parameter('ampl_729_singlepass1_pct',
                                                                  group='beams.ampl_pct', override=True,
                                                                  conversion_function=pct_to_asf)

            self.singlepass0_default_att_mu = self.get_parameter('att_729_singlepass0_db',
                                                                 group='beams.att_db', override=True,
                                                                 conversion_function=att_to_mu)

            self.singlepass1_default_att_mu = self.get_parameter('att_729_singlepass1_db',
                                                                 group='beams.att_db', override=True,
                                                                 conversion_function=att_to_mu)

            self.singlepass0_default_freq_MHz = self.get_parameter('freq_729_singlepass0_mhz',
                                                                   group='beams.freq_mhz', override=True,
                                                                   )

            self.singlepass1_default_freq_MHz = self.get_parameter('freq_729_singlepass1_mhz',
                                                                   group='beams.freq_mhz', override=True,
                                                                   )

            # self.singlepass0_default_freq_MHz = 120.339
            # self.singlepass1_default_freq_MHz = 120.339
            self.singlepass0_default_freq_ftw = self.singlepass0.frequency_to_ftw(self.singlepass0_default_freq_MHz*MHz)
            self.singlepass1_default_freq_ftw = self.singlepass1.frequency_to_ftw(self.singlepass1_default_freq_MHz*MHz)
            # self.singlepass0_default_amp_asf = pct_to_asf(50.)
            # self.singlepass1_default_amp_asf = pct_to_asf(0.01)
            # self.singlepass0_default_att_mu = att_to_mu(7.*dB)
            # self.singlepass1_default_att_mu = att_to_mu(31.5*dB)
