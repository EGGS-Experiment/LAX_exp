import numpy as np
from artiq.experiment import *
from artiq.coredevice.ad9910 import _AD9910_REG_CFR1, PHASE_MODE_CONTINUOUS

from LAX_exp.extensions import *
from LAX_exp.system.subsequences import SqueezeConfigurable
import LAX_exp.experiments.eggs_heating.IonSpectrumAnalyzer as IonSpectrumAnalyzer


class IonSpectrumAnalyzerDDS(IonSpectrumAnalyzer.IonSpectrumAnalyzer):
    """
    Experiment: Ion Spectrum Analyzer DDS

    ***todo: redocument***
    Cool the ions to the ground state of motion via sideband cooling,
    then apply bichromatic heating tones, and try to read out the fluorescence.
    """
    name = 'Ion Spectrum Analyzer DDS'


    def build_experiment(self):
        # run regular sideband cooling build
        super().build_experiment()

        # Ion Spectrum Analyzer
        self.setattr_argument("freq_ionSpecAnal_sideband_offset_khz_list",  Scannable(
                                                                                default=[
                                                                                    ExplicitScan([0]),
                                                                                    CenterScan(0, 5, 0.5, randomize=True)
                                                                                ],
                                                                                global_min=-8000, global_max=8000, global_step=10,
                                                                                unit="kHz", scale=1, ndecimals=3
                                                                            ), group=self.name)
        self.setattr_argument("time_readout_us_list",                       Scannable(
                                                                                default=[
                                                                                    ExplicitScan([95.]),
                                                                                    RangeScan(0, 1000, 200, randomize=True)
                                                                                ],
                                                                                global_min=1, global_max=100000, global_step=1,
                                                                                unit="us", scale=1, ndecimals=5
                                                                            ), group=self.name)

        # ISA - antisqueezing
        self.setattr_argument("enable_ISA_antisqueezing",                   BooleanValue(default=False), group='ISA.antisqueezing')
        self.setattr_argument("ampl_ISA_antisqueezing_rsb_pct",             NumberValue(default=40., ndecimals=2, step=10, min=0.0, max=99), group='ISA.antisqueezing')
        self.setattr_argument("ampl_ISA_antisqueezing_bsb_pct",             NumberValue(default=40., ndecimals=2, step=10, min=0.0, max=99), group='ISA.antisqueezing')
        self.setattr_argument("phase_ISA_antisqueezing_rsb_turns",          NumberValue(default=0.5, ndecimals=3, step=0.1, min=-1., max=1.), group='ISA.antisqueezing')
        self.setattr_argument("phase_ISA_antisqueezing_bsb_turns",          NumberValue(default=0.5, ndecimals=3, step=0.1, min=-1., max=1.), group='ISA.antisqueezing')

        self.setattr_argument("enable_ISA_antisqueezing_dipole",            BooleanValue(default=False), group='ISA.antisqueezing')
        self.setattr_argument("ampl_ISA_antisqueezing_dipole_pct",          NumberValue(default=20., ndecimals=2, step=10, min=0.0, max=99), group='ISA.antisqueezing')
        self.setattr_argument("phase_ISA_antisqueezing_dipole_turns",       NumberValue(default=0.5, ndecimals=3, step=0.1, min=-1., max=1.), group='ISA.antisqueezing')

        # extrinsic squeezing
        self.setattr_argument("freq_squeeze_khz",                           NumberValue(default=1542.2, ndecimals=3, step=10, min=1, max=400000), group='squeeze_configurable')
        self.setattr_argument("phase_antisqueeze_turns_list",               Scannable(
                                                                                    default=[
                                                                                        ExplicitScan([0.]),
                                                                                        RangeScan(0, 1.0, 6, randomize=True)
                                                                                    ],
                                                                                    global_min=0.0, global_max=1.0, global_step=1,
                                                                                    unit="turns", scale=1, ndecimals=3
                                                                                ), group='squeeze_configurable')
        self.setattr_argument("time_squeeze_us",                            NumberValue(default=50., ndecimals=3, step=100, min=1, max=1000000), group='squeeze_configurable')

        # tmp remove
        self.setattr_device("urukul0_cpld")
        self.setattr_device("urukul0_ch1")
        self.setattr_device("urukul0_ch2")
        self.setattr_device("urukul0_ch3")

        self.setattr_device("urukul1_cpld")
        self.setattr_device("urukul1_ch1")
        self.setattr_device("urukul1_ch2")
        self.setattr_device("urukul1_ch3")


    def prepare_experiment(self):
        super().prepare_experiment()

        # configure DDSs for ISA
        self._prepare_DDS_ISA()

    def _prepare_DDS_ISA(self):
        """
        todo: document
        """
        # alias urukuls for ease of use
        self.dds_cpld =                         self.urukul0_cpld
        self.dds =                              self.urukul0_ch1

        # attenuation
        self.att_eggs_heating_mu =              self.dds_cpld.att_to_mu(self.att_eggs_heating_db * dB)
        # preallocate register storage for urukul attenuation register
        self._reg_att_urukul0 =                 np.int32(0)
        self._reg_att_urukul1 =                 np.int32(0)

        # prepare an empty, 0 amplitude waveform
        self._freq_empty_ftw =                  self.dds.frequency_to_ftw(300 * MHz)

        # get latency values from dashboard
        # global urukul1 latency/phase compensation values
        self.phase_urukul1_adjust_turns =       self.get_dataset('dds.delay.phase_urukul1_adjust_turns')
        self.time_urukul1_system_latency_ns =   self.get_dataset('dds.delay.time_urukul1_system_latency_ns')
        # individual urukul1 channel latency compensation values
        self.time_urukul0_ch2_latency_ns =      self.get_dataset('dds.delay.time_urukul0_ch2_latency_ns')
        self.time_urukul0_ch3_latency_ns =      self.get_dataset('dds.delay.time_urukul0_ch3_latency_ns')

        self.time_urukul1_ch1_latency_ns =      self.get_dataset('dds.delay.time_urukul1_ch1_latency_ns') + self.time_urukul1_system_latency_ns
        self.time_urukul1_ch2_latency_ns =      self.get_dataset('dds.delay.time_urukul1_ch2_latency_ns') + self.time_urukul1_system_latency_ns
        self.time_urukul1_ch3_latency_ns =      self.get_dataset('dds.delay.time_urukul1_ch3_latency_ns') + self.time_urukul1_system_latency_ns

    @property
    def results_shape(self):
        return (self.repetitions * len(self.config_eggs_heating_list),
                8)

    @kernel(flags={"fast-math"})
    def initialize_experiment(self):
        # record general subsequences onto DMA
        self.initialize_subsequence.record_dma()
        self.sidebandcool_subsequence.record_dma()
        self.sidebandreadout_subsequence.record_dma()
        self.readout_subsequence.record_dma()

        # tmp remove
        self.ttl8.off()
        self.ttl9.off()
        # tmp remove


        '''PREPARE - DDS CPLD'''
        # close all internal switches
        delay_mu(100)
        with parallel:
            self.urukul0_ch1.sw.off()
            self.urukul0_ch2.sw.off()
            self.urukul0_ch3.sw.off()

            self.urukul1_ch1.sw.off()
            self.urukul1_ch2.sw.off()
            self.urukul1_ch3.sw.off()

        at_mu(now_mu() + 25000)
        self._reg_att_urukul0 = self.urukul0_cpld.get_att_mu()
        self._reg_att_urukul0 &= (0xFF << 0)
        self._reg_att_urukul0 |= ((self.att_eggs_heating_mu << 8) |
                                  (self.att_eggs_heating_mu << 16) |
                                  (self.att_eggs_heating_mu << 24))

        at_mu(now_mu() + 25000)
        self._reg_att_urukul1 = self.urukul1_cpld.get_att_mu()
        self._reg_att_urukul1 &= (0xFF << 0)
        self._reg_att_urukul1 |= ((self.att_eggs_heating_mu << 8) |
                                  (self.att_eggs_heating_mu << 16) |
                                  (self.att_eggs_heating_mu << 24))

        at_mu(now_mu() + 20000)
        self.urukul0_ch1.set_phase_mode(PHASE_MODE_CONTINUOUS)
        self.urukul0_ch2.set_phase_mode(PHASE_MODE_CONTINUOUS)
        self.urukul0_ch3.set_phase_mode(PHASE_MODE_CONTINUOUS)

        self.urukul1_ch1.set_phase_mode(PHASE_MODE_CONTINUOUS)
        self.urukul1_ch2.set_phase_mode(PHASE_MODE_CONTINUOUS)
        self.urukul1_ch3.set_phase_mode(PHASE_MODE_CONTINUOUS)


        '''PREPARE - DDS WAVEFORMS'''
        # set 0 amplitude waveforms to prevent leakage during DDS RF5's use of profile 0
        at_mu(now_mu() + 20000)
        self.urukul0_ch1.set_mu(self._freq_empty_ftw, asf=0x0, profile=0, pow_=0x0, phase_mode=PHASE_MODE_CONTINUOUS)
        self.urukul0_ch2.set_mu(self._freq_empty_ftw, asf=0x0, profile=0, pow_=0x0, phase_mode=PHASE_MODE_CONTINUOUS)
        self.urukul0_ch3.set_mu(self._freq_empty_ftw, asf=0x0, profile=0, pow_=0x0, phase_mode=PHASE_MODE_CONTINUOUS)

        at_mu(now_mu() + 20000)
        self.urukul1_ch1.set_mu(self._freq_empty_ftw, asf=0x0, profile=0, pow_=0x0, phase_mode=PHASE_MODE_CONTINUOUS)
        self.urukul1_ch2.set_mu(self._freq_empty_ftw, asf=0x0, profile=0, pow_=0x0, phase_mode=PHASE_MODE_CONTINUOUS)
        self.urukul1_ch3.set_mu(self._freq_empty_ftw, asf=0x0, profile=0, pow_=0x0, phase_mode=PHASE_MODE_CONTINUOUS)


        '''PREPARE - DDS REGISTERS'''
        at_mu(now_mu() + 10000)
        self.urukul0_ch1.write32(_AD9910_REG_CFR1, (1 << 16) | (1 << 13))
        self.urukul0_ch2.write32(_AD9910_REG_CFR1, (1 << 16) | (1 << 13))
        self.urukul0_ch3.write32(_AD9910_REG_CFR1, (1 << 16) | (1 << 13))

        at_mu(now_mu() + 10000)
        self.urukul1_ch1.write32(_AD9910_REG_CFR1, (1 << 16) | (1 << 13))
        self.urukul1_ch2.write32(_AD9910_REG_CFR1, (1 << 16) | (1 << 13))
        self.urukul1_ch3.write32(_AD9910_REG_CFR1, (1 << 16) | (1 << 13))

        at_mu(now_mu() + 10000)
        self.urukul0_ch1.set_cfr2(matched_latency_enable=1)
        self.urukul0_ch2.set_cfr2(matched_latency_enable=1)
        self.urukul0_ch3.set_cfr2(matched_latency_enable=1)

        at_mu(now_mu() + 10000)
        self.urukul1_ch1.set_cfr2(matched_latency_enable=1)
        self.urukul1_ch2.set_cfr2(matched_latency_enable=1)
        self.urukul1_ch3.set_cfr2(matched_latency_enable=1)


    # MAIN SEQUENCE
    @kernel(flags={"fast-math"})
    def run_main(self):
        self.core.reset()

        # MAIN LOOP
        for trial_num in range(self.repetitions):

            # sweep eggs rf configurations
            for config_vals in self.config_eggs_heating_list:

                '''CONFIGURE'''
                # extract values from config list
                freq_readout_ftw =          np.int32(config_vals[0])
                carrier_freq_hz =           config_vals[1]
                sideband_freq_hz =          config_vals[2]
                time_readout_mu =           np.int64(config_vals[3])
                ampl_rsb_frac =             config_vals[4]
                ampl_bsb_frac =             config_vals[5]
                ampl_dd_frac =              config_vals[6]
                offset_freq_hz =            config_vals[7]
                phase_antisqueeze_pow =     np.int32(config_vals[8])
                phase_rsb_turns =           config_vals[9]
                self.core.break_realtime()

                # configure EGGS tones and set readout
                self.dds_configure(carrier_freq_hz, sideband_freq_hz, offset_freq_hz, phase_rsb_turns, ampl_rsb_frac, ampl_bsb_frac, ampl_dd_frac)
                self.core.break_realtime()
                self.qubit.set_mu(freq_readout_ftw, asf=self.sidebandreadout_subsequence.ampl_sideband_readout_asf, profile=0, phase_mode=PHASE_MODE_CONTINUOUS)
                self.core.break_realtime()


                '''STATE PREPARATION'''
                # initialize ion in S-1/2 state
                self.initialize_subsequence.run_dma()
                # sideband cool
                self.sidebandcool_subsequence.run_dma()

                '''ION SPECTRUM ANALYZER'''
                self.dds_run()

                '''READOUT'''
                # set readout waveform for qubit
                self.qubit.set_profile(0)
                self.qubit.set_att_mu(self.sidebandreadout_subsequence.att_sideband_readout_mu)

                # population transfer pulse
                self.qubit.on()
                delay_mu(time_readout_mu)
                self.qubit.off()

                # read out fluorescence
                self.readout_subsequence.run_dma()

                # update dataset
                with parallel:
                    self.update_results(
                        freq_readout_ftw,
                        self.readout_subsequence.fetch_count(),
                        carrier_freq_hz,
                        sideband_freq_hz,
                        offset_freq_hz,
                        time_readout_mu,
                        phase_antisqueeze_pow,
                        phase_rsb_turns
                    )
                    self.core.break_realtime()

            # rescue ion as needed
            self.rescue_subsequence.run(trial_num)

            # support graceful termination
            with parallel:
                self.check_termination()
                self.core.break_realtime()

        # CLEANUP
        self.core.break_realtime()
        # reset all oscillator frequencies and amplitudes
        self.phaser_eggs.reset_oscillators()
        self.core.break_realtime()
        # set max attenuations for phaser outputs to reduce effect of internal noise
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.channel[0].set_att(31.5 * dB)
        at_mu(self.phaser_eggs.get_next_frame_mu())
        self.phaser_eggs.channel[1].set_att(31.5 * dB)
        self.core.break_realtime()


    # HELPER FUNCTIONS
    @kernel(flags={"fast-math"})
    def dds_configure(self, carrier_freq_hz: TFloat, sideband_freq_hz: TFloat, offset_freq_hz: TFloat,
                      phase_rsb_turns: TFloat,
                      ampl_rsb_frac: TFloat, ampl_bsb_frac: TFloat, ampl_carrier_frac: TFloat):
        """
        Configure the tones on Urukul for ISA.
        Create an RSB, BSB, and carrier signal

        Arguments:
            carrier_freq_hz         (float)     : the maximum waiting time (in machine units) for the trigger signal.
            sideband_freq_hz        (float)     : the holdoff time (in machine units)
            phase_rsb_turns         (float)     : the phase for the rsb tone (in turns)
        """
        # precalculate frequency values for each channel
        freq_rsb_hz =           carrier_freq_hz - sideband_freq_hz + offset_freq_hz
        freq_bsb_hz =           carrier_freq_hz + sideband_freq_hz + offset_freq_hz

        # calculate phase delays for each channel to account for inherent update latencies and system latencies
        # channel 0 (RSB)
        self.phase_ch0_osc0 =   phase_rsb_turns
        self.phase_ch1_osc0 =   (freq_rsb_hz * (self.time_urukul1_ch1_latency_ns * ns)
                                 + self.phase_urukul1_adjust_turns
                                 + phase_rsb_turns)

        # channel 1 (BSB)
        self.phase_ch0_osc1 =   (freq_bsb_hz * (self.time_urukul0_ch2_latency_ns * ns)
                                 + self.phase_eggs_heating_bsb_turns)
        self.phase_ch1_osc1 =   (freq_bsb_hz * (self.time_urukul1_ch2_latency_ns * ns)
                                 + self.phase_urukul1_adjust_turns
                                 + self.phase_eggs_heating_bsb_turns)

        # channel 2 (carrier) (note: ch1 has 0.5 turns to put carrier in dipole config)
        self.phase_ch0_osc2 =   carrier_freq_hz * (self.time_urukul0_ch3_latency_ns * ns)
        self.phase_ch1_osc2 =   (carrier_freq_hz * (self.time_urukul1_ch3_latency_ns * ns)
                                 + self.phase_urukul1_adjust_turns
                                 + 0.5)
        self.core.break_realtime()

        # set desired waveforms
        at_mu(now_mu() + 10000)
        with parallel:
            self.urukul0_ch1.set(carrier_freq_hz - sideband_freq_hz + offset_freq_hz, amplitude=ampl_rsb_frac, phase=self.phase_ch0_osc0,
                                    profile=3, phase_mode=PHASE_MODE_CONTINUOUS)
            self.urukul1_ch1.set(carrier_freq_hz - sideband_freq_hz + offset_freq_hz, amplitude=ampl_rsb_frac, phase=self.phase_ch1_osc0,
                                    profile=3, phase_mode=PHASE_MODE_CONTINUOUS)

        with parallel:
            self.urukul0_ch2.set(carrier_freq_hz + sideband_freq_hz + offset_freq_hz, amplitude=ampl_bsb_frac, phase=self.phase_ch0_osc1,
                                    profile=3, phase_mode=PHASE_MODE_CONTINUOUS)
            self.urukul1_ch2.set(carrier_freq_hz + sideband_freq_hz + offset_freq_hz, amplitude=ampl_bsb_frac, phase=self.phase_ch1_osc1,
                                    profile=3, phase_mode=PHASE_MODE_CONTINUOUS)

        with parallel:
            self.urukul0_ch3.set(carrier_freq_hz, amplitude=ampl_carrier_frac, phase=self.phase_ch0_osc2,
                                    profile=3, phase_mode=PHASE_MODE_CONTINUOUS)
            self.urukul1_ch3.set(carrier_freq_hz, amplitude=ampl_carrier_frac, phase=self.phase_ch1_osc2,
                                    profile=3, phase_mode=PHASE_MODE_CONTINUOUS)


        '''PREPARE - DDS REGISTERS'''
        at_mu(now_mu() + 10000)
        self.urukul0_ch1.write32(_AD9910_REG_CFR1, (1 << 16) | (1 << 13))
        self.urukul0_ch2.write32(_AD9910_REG_CFR1, (1 << 16) | (1 << 13))
        self.urukul0_ch3.write32(_AD9910_REG_CFR1, (1 << 16) | (1 << 13))

        at_mu(now_mu() + 10000)
        self.urukul1_ch1.write32(_AD9910_REG_CFR1, (1 << 16) | (1 << 13))
        self.urukul1_ch2.write32(_AD9910_REG_CFR1, (1 << 16) | (1 << 13))
        self.urukul1_ch3.write32(_AD9910_REG_CFR1, (1 << 16) | (1 << 13))


    @kernel(flags={"fast-math"})
    def dds_run(self):
        """
        Run the ISA via Urukuls.
        """
        # ensure attenuations are set correctly
        self.urukul0_cpld.set_all_att_mu(self._reg_att_urukul0)
        self.urukul1_cpld.set_all_att_mu(self._reg_att_urukul1)

        # note: we coarse align to previous SYNC_CLK period
        time_start_mu = now_mu() & ~7

        '''PULSE START'''
        # set start (active) profile
        at_mu(time_start_mu)
        with parallel:
            self.urukul0_cpld.set_profile(3)
            self.urukul1_cpld.set_profile(3)

        # open RF switches early since they have ~100 ns rise time
        at_mu(time_start_mu + ((416 + 63) - 210))
        with parallel:
            self.urukul0_ch1.sw.on()
            self.urukul0_ch2.sw.on()
            self.urukul0_ch3.sw.on()

            self.urukul1_ch1.sw.on()
            self.urukul1_ch2.sw.on()
            self.urukul1_ch3.sw.on()

        # send trigger when waveform begins
        at_mu(time_start_mu + (416 + 63))
        self.ttl8.on()
        delay_mu(self.time_eggs_heating_mu)


        '''PULSE STOP'''
        # close RF switches
        with parallel:
            self.urukul0_ch1.sw.off()
            self.urukul0_ch2.sw.off()
            self.urukul0_ch3.sw.off()

            self.urukul1_ch1.sw.off()
            self.urukul1_ch2.sw.off()
            self.urukul1_ch3.sw.off()

            self.ttl8.off()

        # switch to 0 amplitude profile
        delay_mu(8)
        with parallel:
            self.urukul0_cpld.set_profile(0)
            self.urukul1_cpld.set_profile(0)
