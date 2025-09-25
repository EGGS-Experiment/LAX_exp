from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import InitializeQubit, Readout, RabiflopReadout
import numpy as np
from artiq.coredevice import ad9910


class DissipativeStatePreparation(LAXExperiment):
    """
    Experiment: Dissipative State Preparation

    From a Doppler-cooled ion apply birchromatic beams to reservoir engineer a motional state
    """
    name = 'Dissipative State Preparation'
    kernel_invariants = {
        # 729 bichromatic
        'dsp_phase_bsb_turns', 'dsp_freq_carrier_MHz', 'dsp_freq_secular_kHz', 'dsp_chromatic_att_dB', 'dsp_ampl_rsb',
        'dsp_freq_rsb_MHz', 'dsp_freq_bsb_MHz', 'dsp_ampl_bsb_pct', 'squeeze_r',

        # 729 bichromatic machine unit params
        'dsp_freq_carrier_ftw', 'dsp_freq_secular_ftw', 'dsp_freq_rsb_ftw', 'dsp_freq_bsb_ftw', 'dsp_chromatic_att_mu',
        'dsp_phase_bsb_pow', 'dsp_ampl_bsb_asf', 'dsp_ampl_rsb_asf'

        # 729 bichromatic default values
                                                 'singlepass0_default_amp_asf', 'singlepass1_default_amp_asf',
        'singlepass0_default_att_mu',
        'singlepass1_default_att_mu', 'singlepass0_default_freq_mhz', 'singlepass1_default_freq_mhz'

        # timing
                                                                      'dsp_time_mu', 'rabiflop_readout_times_us_list',
        'rabiflop_readout_times_mu_list',

        # 729 double pass
        'dsp_douplepass_carrier_att_mu', 'dsp_douplepass_carrier_asf'

        # subsequences
                                         'readout_subsequence', 'initialize_subsequence',
        'rabiflop_readout_subsequence',
        'spin_polarization_re_subsequence'

        # devices
        'singlepass0', 'singlepass1',
    }

    def build_experiment(self):
        # get devices

        # get arguments
        self.setattr_argument('squeeze_r',
                              NumberValue(default=1., precision=3, step=0.01, min=0.,
                                          max=5.), group=self.name,
                              tooltip='parameter determining amount of squeeze after reservoir engineering')

        self.setattr_argument('dsp_time_us', NumberValue(default=125, precision=3, step=10, min=1,
                                                         max=1000000), group=self.name,
                              tooltip='How long to apply bichromatic beams for reservoir engineering')

        # set frequencies
        self.setattr_argument('dsp_freq_carrier_MHz', NumberValue(default=101.09, precision=3, step=1e-3, max=120.,
                                                                  min=80., unit="MHz"), group=self.name,
                              tooltip='frequency of the urukul channel used to drive rthe red sideband')

        self.setattr_argument('dsp_freq_secular_kHz', NumberValue(default=702.0, precision=3, step=1e-3, max=120.,
                                                                  min=80., unit="kHz"), group=self.name,
                              tooltip='frequency of the urukul channel used to drive the blue sideband')

        # set attenuations
        self.setattr_argument('dsp_chromatic_att_dB', NumberValue(default=12, precision=1, step=1e-3, max=31.5,
                                                                  min=13., unit="dB"), group=self.name,
                              tooltip='attenuation to apply to both of the urukuls channels (i.e the urukul channels used to '
                                      'drive the red and blue sideband are set to the same attenuation')

        # set phases
        self.setattr_argument('dsp_phase_bsb_turns', NumberValue(default=0., precision=2, step=0.01, max=1.,
                                                                 min=0., unit="turns"), group=self.name,
                              tooltip='phase offset between bsb and rsb tones')

        # set amplitudes
        self.setattr_argument('dsp_ampl_rsb',
                              NumberValue(default=50., step=5., min=0.01, max=50., precision=1, unit='%'),
                              group=self.name,
                              tooltip='amplitude of drive of urukul channel driving red sideband. Blue sideband amplitude will be '
                                      'calculated from this value')

        # set timing
        # rabi flopping arguments
        self.setattr_argument("rabiflop_readout_times_us_list", Scannable(
            default=[
                RangeScan(1, 100, 100, randomize=True),
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

        self.get_default_values()

        # get red and blue sideband values
        self.dsp_freq_rsb_MHz = (self.singlepass0_default_freq_ftw * MHz - self.freq_secular_kHz * kHz) / MHz
        self.dsp_freq_bsb_MHz = (self.singlepass1_default_freq_ftw * MHz + self.freq_secular_kHz * kHz) / MHz
        # convert parameters to machine units
        self.dsp_freq_carrier_ftw = self.singlepass0.frequency_to_ftw(self.dsp_freq_carrier_MHz * MHz)
        self.dsp_freq_secular_ftw = self.singlepass1.frequency_to_ftw(self.dsp_freq_secular_MHz * MHz)
        self.dsp_freq_rsb_ftw = self.singlepass0.frequency_to_ftw(self.dsp_freq_rsb_MHz * MHz)
        self.dsp_freq_bsb_ftw = self.singlepass1.frequency_to_ftw(self.dsp_freq_bsb_MHz * MHz)

        self.dsp_time_mu = self.core.seconds_to_mu(self.dsp_time_us * us)

        self.dsp_chromatic_att_mu = att_to_mu(self.dsp_att_dB * dB)
        self.dsp_douplepass_carrier_att_mu = att_to_mu(8. * dB)

        self.dsp_phase_bsb_pow = self.singlepass1.turns_to_pow(self.dsp_phase_bsb_turns)

        # find bsb amplitude from rsb amplitude
        self.dsp_ampl_bsb_pct = np.sqrt(np.tanh(self.dsp_r)) * self.dsp_ampl_rsb_pct

        self.dsp_ampl_rsb_asf = pct_to_asf(self.dsp_ampl_rsb_pct)
        self.dsp_ampl_bsb_asf = pct_to_asf(self.dsp_ampl_bsb_pct)
        self.dsp_douplepass_carrier_asf = pct_to_asf(50.)

        # get time values
        self.rabiflop_readout_times_mu_list = [self.core.seconds_to_mu(self.rabiflop_readout_time_us * us)
                                               for rabiflop_readout_time_us in rabiflop_readout_times_list_us]

    @property
    def results_shape(self):
        return (self.repetitions * self.sub_repetitions * len(self.rabiflop_time_mu_list),
                4)

    def initialize_experiment(self) -> TNone:

        self.core.break_realtime()

        self.qubit.set_mu(ftw=self.dsp_freq_carrier_ftw, ampl=self.dsp_douplepass_carrier_asf)
        self.qubit.set_att_mu(self.dsp_douplepass_carrier_att_mu)

        self.singlepass0.set_att_mu(self.dsp_chromatic_att_mu)
        self.singlepass1.set_att_mu(self.dsp_chromatic_att_mu)

        self.initialize_subsequence.record_dma()
        self.readout_subsequence.record_dma()
        self.spin_polarization_re_subsequence.record_dma()
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def run(self, start_time) -> TNone:
        self.core.break_realtime()

        for trial_num in range(self.repetitions):

            for rabiflop_readout_time_mu in self.rabiflop_readout_times_mu_list:
                self.core.break_realtime()

                self.initialize_subsequence.run_dma()

                ref_time_mu = now_mu() & ~0x7
                self.configure_bichromatic(ref_time_mu)

                self.pulse_dsp()

                # ensure everything is in ground state by optical pumping
                self.spin_polarization_re_subsequence.run_dma()

                self.rabiflop_readout_subsequence.run_time(rabiflop_readout_time_mu)
                self.readout_subsequence.run_dma()

    @kernel(flags={"fast-math"})
    def cleanup_experiment(self) -> TNone:

        # reset singlepass values
        self.singlepass0.set_mu(ftw=self.singlepass0_default_freq_ftw, pow=0.,
                                asf=self.singlepass0_default_amp_asf,
                                )
        self.singlepass0.set_att_mu(self.singlepass0_default_att_mu)

        self.singlepass1.set_mu(ftw=self.singlepass1_default_freq_ftw, pow=0.,
                                asf=self.singlepass1_default_amp_asf,
                                )
        self.singlepass1.set_att_mu(self.singlepass1_default_att_mu)

        # turn off lasers
        self.qubit.off()
        self.probe.off()
        self.repump_cooling.off()
        self.repump_qubit.off()

    @kernel(flags={"fast-math"})
    def configure_bichromatic(self, ref_time_mu: TInt64):
        self.singlepass0.set_cfr1()
        self.singlepass1.set_cfr1()
        self.qubit.cpld.io_update.pulse_mu()
        self.singlepass0.set_mu(ftw=self.dsp_freq_rsb_ftw, pow=0.,
                                asf=self.dsp_ampl_rsb_asf,
                                phase_mode=ad9910.PHASE_MODE_TRACKING,
                                ref_time_mu=ref_time_mu
                                )

        self.singlepass1.set_mu(ftw=self.dsp_freq_bsb_ftw,
                                pow=self.dsp_phase_bsb_pow,
                                asf=self.dsp_ampl_bsb_asf,
                                phase_mode=ad9910.PHASE_MODE_TRACKING,
                                ref_time_mu=ref_time_mu
                                )

    @kernel(flags={"fast-math"})
    def pulse_dsp(self) -> TNone:
        # ensure carrier is on
        self.qubit.on()

        # turn on dissipation lasers
        self.repump_cooling.on()
        self.repump_qubit.on()
        self.probe.on()

        # turn on 729
        self.singlepass0.sw.on()
        self.singlepass1.sw.on()

        delay_mu(self.dsp_time_mu)
        self.singlepass0.sw.off()
        self.singlepass1.sw.off()

        # turn off all the lasers
        self.qubit.of()
        self.repump_cooling.off()
        self.repump_qubit.off()
        self.probe.off()

    def get_default_values(self) -> TNone:
        # get default values for the 729 single pass
        self.singlepass0_default_amp_asf = self.get_parameter('ampl_bichromatic_0_pct',
                                                              group='beams.ampl_pct', override=True,
                                                              conversion_function=pct_to_asf)

        self.singlepass1_default_amp_asf = self.get_parameter('ampl_bichromatic_1_pct',
                                                              group='beams.ampl_pct', override=True,
                                                              conversion_function=pct_to_asf)

        self.singlepass0_default_att_mu = self.get_parameter('att_bichromatic_0_dB',
                                                             group='beams.att_dB', override=True,
                                                             conversion_function=att_to_mu)

        self.singlepass1_default_att_mu = self.get_parameter('att_bichromatic_1_dB',
                                                             group='beams.att_dB', override=True,
                                                             conversion_function=att_to_mu)

        self.singlepass0_default_freq_mhz = self.get_parameter('freq_bichromatic_0_ftw',
                                                               group='beams.freq_MHz', override=True,
                                                               )

        self.singlepass1_default_freq_mhz = self.get_parameter('freq_bichromatic_1_ftw',
                                                               group='beams.freq_MHz', override=True,
                                                               )
