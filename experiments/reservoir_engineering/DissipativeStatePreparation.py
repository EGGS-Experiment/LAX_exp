from artiq.coredevice.ad9914 import PHASE_MODE_CONTINUOUS, PHASE_MODE_TRACKING
from artiq.experiment import *

from LAX_exp.extensions import *
from LAX_exp.language import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import InitializeQubit, Readout, RabiflopReadout, SpinPolarizationRE, SidebandReadout
from LAX_exp.system.subsequences import SidebandCoolContinuousRAM

import numpy as np
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
        'dsp_phase_bsb_turns', 'dsp_freq_carrier_MHz', 'dsp_freq_secular_kHz_list', 'dsp_chromatic_att_dB',
        'dsp_ampl_rsb_pct', 'dsp_ampl_bsb_pct', 'dsp_squeeze_r',

        # 729 multi-chromatic machine unit params
        # 'dsp_freq_carrier_ftw', 'dsp_freq_secular_ftw_list',
        'dsp_chromatic_att_mu',
        'dsp_phase_bsb_pow', 'dsp_phase_rsb_pow', 'dsp_ampl_bsb_asf', 'dsp_ampl_rsb_asf',

        # 729 multi-chromatic default values
        'singlepass0_default_amp_asf', 'singlepass1_default_amp_asf', 'singlepass0_default_att_mu',
        'singlepass1_default_att_mu', 'singlepass0_default_freq_MHz', 'singlepass1_default_freq_MHz',

        # timing
        'dsp_time_mu', 'rabiflop_readout_times_us_list', 'rabiflop_readout_times_mu_list',

        # 729 double pass
        'dsp_doublepass_att_mu', 'dsp_doublepass_asf',

        # subsequences
        'readout_subsequence', 'initialize_subsequence', 'rabiflop_readout_subsequence',
        'spin_polarization_re_subsequence',

        # devices
        'singlepass0', 'singlepass1',
    }

    def build_experiment(self):

        # allocate profiles
        self.profile_readout = 0
        self.profile_SBC = 1
        self.profile_dsp = 2

        # get arguments
        self.setattr_argument('repetitions',
                              NumberValue(default=50., precision=0, step=1, min=1.,
                                          max=1e6), group=self.name,
                              tooltip='parameter determining level of squeezing after reservoir engineering')

        # dissipative state preparation parameters
        self.setattr_argument('dsp_squeeze_r',
                              NumberValue(default=1., precision=3, step=0.01, min=0.,
                                          max=5.), group=self.name,
                              tooltip='parameter determining amount of squeeze after reservoir engineering')

        self.setattr_argument('dsp_time_us', NumberValue(default=8000, precision=3, step=10, min=1,
                                                         max=1000000), group=self.name,
                              tooltip='How long to apply multi-chromatic beams for reservoir engineering')

        # AOM frequenices
        self.setattr_argument('dsp_freq_carrier_MHz', NumberValue(default=101.092, precision=3, step=1e-3, max=120.,
                                                                  min=80., unit="MHz"), group=self.name,
                              tooltip='Frequency of the carrier transition in the calcium ion')

        self.setattr_argument('dsp_freq_secular_kHz_list', Scannable(
            default=[
                # CenterScan(701., 60., 0.25, randomize=True),
                # RangeScan(, 200, 200, randomize=True),
                ExplicitScan([702]),
            ],
            global_min=0.1, global_max=5000, global_step=1e-3,
            unit=kHz, scale=1, precision=5
        ), tooltip='difference in frequency between the carrier and the sidebands',
                              group=self.name)

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
                              NumberValue(default=50., step=0.1, min=0.01, max=50., precision=3, unit="%", scale=1),
                              group=self.name,
                              tooltip='amplitude of drive of urukul channel driving red sideband. Blue sideband amplitude will be '
                                      'calculated from this value')

        self.dsp_doublepass_amp_pct = self.setattr_argument('dsp_doublepass_amp_pct',
                                                            NumberValue(default=50., step=0.5, min=0.01, max=50.,
                                                                        precision=3, unit="%", scale=1),
                                                            group=self.name)

        self.dsp_quench_amp_pct = self.setattr_argument('dsp_quench_amp_pct',
                                                            NumberValue(default=12.28, step=0.01, min=0.01, max=50.,
                                                                        precision=3, unit="%", scale=1),
                                                            group=self.name)

        self.dsp_doublepass_carrier_att_dB = self.setattr_argument('dsp_doublepass_att_dB', NumberValue(8.,
                                                                                                        min=8.,
                                                                                                        max=31.5,
                                                                                                        step=0.5,
                                                                                                        precision=1,
                                                                                                        unit='dB'),
                                                                   group=self.name)

        self.num_spinpol_cycles = self.setattr_argument('num_spinpol_cycles', NumberValue(20, min=0, max=100,
                                                                                          step=1, precision=0),
                                                        group=self.name,
                                                        tooltip='how many times to perform spin polarization during dissipative state prep')

        # sideband readout
        self.setattr_argument('sidebandreadout_time_us', NumberValue(default=27.4, precision=3, step=10, min=1,
                                                                     max=1000), group='sideband_readout',
                              tooltip='how long to readout readout sidebands')

        # rabi flopp readout
        self.setattr_argument('enable_rabiflop_readout', BooleanValue(True),
                              tooltip='true/false value indicating if rabi flopping should be down after reservoir engineering',
                              group='rabiflop_readout')
        self.setattr_argument("rabiflop_readout_times_us_list", Scannable(
            default=[
                RangeScan(1, 1500, 500, randomize=True),
                # ExplicitScan([1]),
                # CenterScan(3.05, 5., 0.1, randomize=True),
            ],
            global_min=1, global_max=100000, global_step=1,
            unit="us", scale=1, precision=5
        ), group='rabiflop_readout', tooltip='times to rabiflop')

        self.setattr_argument('rabiflop_readout_freq_MHz', NumberValue(default=101.431, min=80., max=120.,
                                                                       precision=5, step=0.0001, unit='MHz'),
                              group='rabiflop_readout',
                              tooltip='AOM frequency to rabi flop at')

        self.default_doublepass_carrier_att_dB = self.setattr_argument('default_doublepass_att_dB', NumberValue(8.,
                                                                                                                min=8.,
                                                                                                                max=31.5,
                                                                                                                step=0.5,
                                                                                                                precision=1,
                                                                                                                unit='dB'),
                                                                       group='defaults')

        self.rabiflop_readout_doublepass_amp_pct = self.setattr_argument('rabiflop_readout_doublepass_amp_pct',
                                                            NumberValue(default=50., step=0.5, min=0.01, max=50.,
                                                                        precision=3, unit="%", scale=1),
                                                            group='rabiflop_readout')

        # use as placeholder
        self.dsp_ampl_bsb_pct = 0.

        # get all necessary devices
        self.setattr_device('qubit')
        self.setattr_device('repump_qubit')
        self.setattr_device('repump_cooling')
        self.setattr_device('probe')
        self.singlepass0 = self.get_device('urukul0_ch1')
        self.singlepass1 = self.get_device('urukul0_ch2')

        # get all necessary subsequences
        self.sidebandreadout_subsequence = SidebandReadout(self, profile_dds=self.profile_readout)
        self.readout_subsequence = Readout(self)
        self.initialize_subsequence = InitializeQubit(self)
        self.rabiflop_readout_subsequence = RabiflopReadout(self)
        self.spin_polarization_re_subsequence = SpinPolarizationRE(self)

        self.sidebandcool_subsequence = SidebandCoolContinuousRAM(
            self, profile_729=self.profile_SBC, profile_854=self.profile_SBC,
            ram_addr_start_729=0, ram_addr_start_854=0, num_samples=200
        )

    def prepare_experiment(self):

        # retrieve parameters that we typically use for the single pass
        self.get_singlepass_default_values()

        # convert parameters to machine units
        self.dsp_freq_carrier_ftw = self.qubit.frequency_to_ftw(self.dsp_freq_carrier_MHz * MHz)
        self.dsp_freq_secular_ftw_list = np.array([self.singlepass1.frequency_to_ftw(dsp_freq_secular_kHz * kHz)
                                                   for dsp_freq_secular_kHz in self.dsp_freq_secular_kHz_list])

        self.rabiflop_readout_freq_ftw = self.qubit.frequency_to_ftw(self.rabiflop_readout_freq_MHz * MHz)

        self.dsp_time_mu = self.core.seconds_to_mu(self.dsp_time_us * us)

        self.dsp_chromatic_att_mu = att_to_mu(self.dsp_chromatic_att_dB * dB)

        self.dsp_doublepass_att_mu = att_to_mu(self.dsp_doublepass_att_dB)
        self.default_doublepass_att_mu = att_to_mu(self.default_doublepass_att_dB)

        self.dsp_phase_rsb_pow = self.singlepass0.turns_to_pow(0.)
        self.dsp_phase_bsb_pow = self.singlepass1.turns_to_pow(self.dsp_phase_bsb_turns)

        # todo: implement interpolation functions to scale frequency effects coupling to fiber
        # todo: implement interpolation functions to scale difference in power from channels

        # find bsb amplitude from rsb amplitude
        self.dsp_ampl_bsb_pct = np.tanh(self.dsp_squeeze_r) * self.dsp_ampl_rsb_pct

        # convert amplitudes into machine units
        self.dsp_ampl_rsb_asf = pct_to_asf(self.dsp_ampl_rsb_pct)
        self.dsp_ampl_bsb_asf = pct_to_asf(self.dsp_ampl_bsb_pct)

        self.dsp_doublepass_asf = pct_to_asf(self.dsp_doublepass_amp_pct)

        # get rabi flop time values
        self.rabiflop_readout_times_mu_list = np.array([self.core.seconds_to_mu(rabiflop_readout_time_us * us)
                                                        for rabiflop_readout_time_us in
                                                        self.rabiflop_readout_times_us_list])

        self.rabiflop_readout_doublepass_asf = pct_to_asf(self.rabiflop_readout_doublepass_amp_pct)

        # spin pol
        self.cycle_time_us = self.dsp_time_us / self.num_spinpol_cycles

        # sideband readout
        self.sidebandreadout_time_mu = self.core.seconds_to_mu(self.sidebandreadout_time_us * us)

        # other lasers
        self.cycle_time_mu = self.core.seconds_to_mu(self.cycle_time_us * us)
        self.ampl_quench_asf = pct_to_asf(self.dsp_quench_amp_pct)
        self.dsp_doublepass_asf = pct_to_asf(self.dsp_doublepass_amp_pct)

        # gather parameters
        self.time_spinpol_mu = self.get_parameter('time_spinpol_us', group='timing', override=True,
                                                  conversion_function=seconds_to_mu, units=us)
        self.freq_repump_qubit_ftw = self.get_parameter('freq_repump_qubit_mhz', group='beams.freq_mhz', override=False,
                                                        conversion_function=hz_to_ftw, units=MHz)

        self.freq_sideband_readout_ftw_list = self.sidebandreadout_subsequence.freq_sideband_readout_ftw_list

        if not self.enable_rabiflop_readout:
            self.rabiflop_readout_times_mu_list = np.array([0])
        else:
            self.freq_sideband_readout_ftw_list = np.array([0])

        self.config_experiment_list = create_experiment_config(self.freq_sideband_readout_ftw_list,
                                                               self.rabiflop_readout_times_mu_list,
                                                               self.dsp_freq_secular_ftw_list,
                                                               config_type=np.int32)

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_experiment_list), 4)

    @kernel(flags={"fast-math"})
    def initialize_experiment(self) -> TNone:

        self.core.break_realtime()

        # ensure we are not autoclearing the phase accumulator
        self.qubit.set_cfr1()
        self.singlepass0.set_cfr1()
        self.singlepass1.set_cfr1()
        self.qubit.cpld.io_update.pulse_mu(8)
        delay_mu(25000)

        # set up singlepass AOMs to default values (b/c AOM thermal drift) on ALL profiles
        for i in range(8):
            self.singlepass0.set_mu(ftw=self.singlepass0_default_freq_ftw,
                                    asf=self.singlepass0_default_amp_asf,
                                    profile=i)

            self.singlepass1.set_mu(ftw=self.singlepass1_default_freq_ftw,
                                    asf=self.singlepass1_default_amp_asf,
                                    profile=i)

            delay_mu(8000)

        self.singlepass0.set_att_mu(self.singlepass0_default_att_mu)
        self.singlepass1.set_att_mu(self.singlepass1_default_att_mu)

        self.singlepass0.sw.on()
        self.singlepass1.sw.off()
        delay_mu(25000)

        self.repump_qubit.set_mu(self.freq_repump_qubit_ftw, asf=self.ampl_quench_asf, phase_mode=PHASE_MODE_CONTINUOUS,
                                 profile=self.profile_dsp)

        # record subsequences onto dma
        self.sidebandcool_subsequence.record_dma()
        self.initialize_subsequence.record_dma()
        self.readout_subsequence.record_dma()
        self.spin_polarization_re_subsequence.record_dma()
        self.sidebandreadout_subsequence.record_dma()
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def run_main(self) -> TNone:
        self.core.break_realtime()

        for trial_num in range(self.repetitions):

            for config_vals in self.config_experiment_list:
                freq_readout_ftw = config_vals[0]
                rabiflop_readout_time_mu = np.int64(config_vals[1])
                secular_freq_ftw = config_vals[2]

                # Doppler cool and optically pump 40Ca+ to |S_1/2, mj=-1/2>
                self.initialize_subsequence.run_dma()
                self.sidebandcool_subsequence.run_dma()

                # get a reference time so all subsequent pulses are phase coherent
                ref_time_mu = now_mu() & ~0x7

                # # configure the single pass values
                self.configure_singlepass(ref_time_mu, secular_freq_ftw)

                # # # execute the dissipative state preparation
                self.pulse_dsp()

                # # ensure everything is in ground state by optical pumping (turn on 397 sigma+, 866, 854)
                self.spin_polarization_re_subsequence.run_dma()

                # # rabi flop
                self.set_default_singlepass_values()
                self.singlepass0.sw.on()
                self.singlepass1.sw.off()
                if self.enable_rabiflop_readout:
                    self.qubit.set_mu(self.rabiflop_readout_freq_ftw, asf=self.rabiflop_readout_doublepass_asf,
                                      profile=self.profile_readout)
                    self.rabiflop_readout_subsequence.run_time(rabiflop_readout_time_mu)

                else:
                    # todo: do we need these delays
                    delay_mu(25000)
                    self.qubit.set_mu(freq_readout_ftw, asf=self.sidebandreadout_subsequence.ampl_sideband_readout_asf,
                                      profile=self.profile_readout)
                    delay_mu(8000)
                    self.sidebandreadout_subsequence.run_time(self.sidebandreadout_time_mu)

                # # # scatter photons off |S_1/2> -> |P_1/2>  transition
                self.readout_subsequence.run_dma()
                counts = self.readout_subsequence.fetch_count()

                # # # store results
                self.update_results(rabiflop_readout_time_mu,
                                    counts,
                                    freq_readout_ftw,
                                    secular_freq_ftw)
                self.check_termination()  # check termination b/c we haven't in a while
                self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def cleanup_experiment(self) -> TNone:

        for i in range(8):
            # reset singlepass values
            self.singlepass0.set_mu(ftw=self.singlepass0_default_freq_ftw,
                                    asf=self.singlepass0_default_amp_asf,
                                    profile=i)

            self.singlepass1.set_mu(ftw=self.singlepass1_default_freq_ftw,
                                    asf=self.singlepass1_default_amp_asf,
                                    profile=i)

            delay_mu(8000)

        # reset singlepass to normal configuration
        self.singlepass0.set_att_mu(self.singlepass0_default_att_mu)
        self.singlepass1.set_att_mu(self.singlepass1_default_att_mu)
        self.singlepass0.sw.on()
        self.singlepass1.sw.off()
        delay_mu(25000)

        # reset lasers
        self.qubit.off()
        self.probe.on()
        self.repump_cooling.on()
        self.repump_qubit.on()

    @kernel(flags={"fast-math"})
    def configure_singlepass(self, ref_time_mu: TInt64, dsp_freq_secular_ftw: TInt32):

        # set singlepass to values for DSP
        self.singlepass0.set_mu(ftw=self.singlepass0_default_freq_ftw - dsp_freq_secular_ftw,
                                pow_=self.dsp_phase_rsb_pow,
                                asf=self.dsp_ampl_rsb_asf,
                                phase_mode=PHASE_MODE_TRACKING,
                                ref_time_mu=ref_time_mu,
                                profile=self.profile_dsp
                                )

        self.singlepass1.set_mu(ftw=self.singlepass1_default_freq_ftw + dsp_freq_secular_ftw,
                                pow_=self.dsp_phase_bsb_pow,
                                asf=self.dsp_ampl_bsb_asf,
                                phase_mode=PHASE_MODE_TRACKING,
                                ref_time_mu=ref_time_mu,
                                profile=self.profile_dsp
                                )

        self.singlepass0.cpld.set_profile(self.profile_dsp)
        self.singlepass0.cpld.io_update.pulse_mu(8)

        self.singlepass0.set_att_mu(self.dsp_chromatic_att_mu)
        self.singlepass1.set_att_mu(self.dsp_chromatic_att_mu)

    @kernel(flags={"fast-math"})
    def set_default_singlepass_values(self):

        # reset singlepass to normal configuration
        self.singlepass0.set_mu(ftw=self.singlepass0_default_freq_ftw,
                                asf=self.singlepass0_default_amp_asf,
                                profile=self.profile_readout,
                                phase_mode=PHASE_MODE_CONTINUOUS
                                )

        self.singlepass1.set_mu(ftw=self.singlepass0_default_freq_ftw,
                                asf=self.singlepass1_default_amp_asf,
                                profile=self.profile_readout,
                                phase_mode=PHASE_MODE_CONTINUOUS)

        self.singlepass0.set_att_mu(self.singlepass0_default_att_mu)
        self.singlepass1.set_att_mu(self.singlepass1_default_att_mu)

        self.singlepass0.sw.on()
        self.singlepass1.sw.off()
        self.singlepass0.cpld.set_profile(self.profile_readout)
        self.singlepass0.cpld.io_update.pulse_mu(8)

    @kernel(flags={"fast-math"})
    def pulse_dsp(self) -> TNone:

        self.repump_qubit.set_profile(self.profile_dsp)
        self.repump_qubit.cpld.io_update.pulse_mu(8)

        self.qubit.set_mu(self.dsp_freq_carrier_ftw, asf=self.dsp_doublepass_asf,
                          profile=self.profile_dsp)
        self.qubit.set_att_mu(self.dsp_doublepass_att_mu)
        self.qubit.cpld.io_update.pulse_mu(8)

        # ensure doublepass is on
        with parallel:
            self.qubit.on()

            # turn on dissipation lasers (leave 397 off)
            self.repump_cooling.on()
            self.repump_qubit.on()

        for i in range(20):

            # turn on single pass
            self.singlepass0.sw.on()
            self.singlepass1.sw.on()

            delay_mu(self.cycle_time_mu)

            # turn off singlepass
            self.singlepass0.sw.off()
            self.singlepass1.sw.off()

            # # repump to |S1/2, mj=-1/2>
            # self.probe.on()
            # delay_mu(self.spinpol_time_mu)
            # # ensure spin pol laser is off to prevent AC stark shifts
            # self.probe.off()
            self.spin_polarization_re_subsequence.run_dma()

        # turn off the 729 and 397 sigma
        with parallel:
            self.qubit.off()
            self.probe.off()

        self.qubit.set_att_mu(self.default_doublepass_att_mu)
        self.repump_qubit.set_profile(self.profile_readout)
        self.repump_qubit.cpld.io_update.pulse_mu(8)
        self.repump_qubit.off()

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

        self.singlepass0_default_freq_ftw = self.singlepass0.frequency_to_ftw(self.singlepass0_default_freq_MHz * MHz)
        self.singlepass1_default_freq_ftw = self.singlepass1.frequency_to_ftw(self.singlepass1_default_freq_MHz * MHz)
