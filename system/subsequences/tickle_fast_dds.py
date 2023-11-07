import numpy as np
from artiq.experiment import *
from artiq.coredevice.urukul import CFG_PROFILE, CFG_IO_UPDATE
from artiq.coredevice.ad9910 import PHASE_MODE_ABSOLUTE, _AD9910_REG_CFR1

from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence

# TODO: INTEGRATE TTL SWITCHES FOR RF
# TODO: ENSURE WE'RE USING SET_CFR1 TO DO PHASE AUTOCLEAR


class TickleFastDDS(LAXSubsequence):
    """
    Subsequence: Tickle Fast DDS

    Heat the ion by applying an RF signal from a DDS at the secular frequency.
    Interferes two channels destructively to achieve fast switching times.
    """
    name = 'tickle_fast_dds'


    def build_subsequence(self):
        self.setattr_argument('att_ticklefast_dds_db', NumberValue(default=24, ndecimals=1, step=0.5, min=0, max=31.5), group='ticklefast')

        # get relevant devices
        self.dds_ch0 = self.get_device('urukul0_ch3')
        self.dds_ch1 = self.get_device('urukul1_ch3')

        # tmp remove
        self.setattr_device('ttl8')
        self.setattr_device('ttl9')
        # tmp remove

    def prepare_subsequence(self):
        # get DDS configuration parameters from dataset manager
        self.ampl_ticklefast_asf =          self.get_parameter('ampl_ticklefast_pct',
                                                               group='dds.ampl_pct', override=False, conversion_function=pct_to_asf)
        self.time_latency_ch1_system_ns =   self.get_parameter('time_urukul0_urukul1_ch3_latency_ns',
                                                               group='dds.delay', override=False)

        # prepare parameters for tickle pulse
        self.att_ticklefast_mu =            att_to_mu(self.att_ticklefast_dds_db * dB)

        # set empty holder variable for delay time
        self.time_delay_mu =                np.int64(0)

        # set empty variables for phase to compensate for ch0 & ch1 delays
        self.phase_ch1_inherent_turns =     0.
        self.phase_ch1_latency_turns =      np.float(0)
        self.phase_ch1_delay_turns =        np.float(0)
        self.phase_ch1_final_pow =          np.int32(0)

        # set parameter for DDS preparation latency before we can
        self.time_system_prepare_delay_mu = np.int64(420) & ~0x7
        self.time_system_cleanup_delay_mu = np.int64(420)

    @kernel(flags={"fast-math"})
    def initialize_subsequence(self):
        # set up DDSs for output
        with parallel:
            self.dds_ch0.set_att_mu(self.att_ticklefast_mu)
            self.dds_ch1.set_att_mu(self.att_ticklefast_mu)

            self.dds_ch0.set_phase_mode(PHASE_MODE_ABSOLUTE)
            self.dds_ch1.set_phase_mode(PHASE_MODE_ABSOLUTE)

            # tmp remove
            self.ttl8.off()
            self.ttl9.off()
            # tmp remove
        self.core.break_realtime()

        with parallel:
            self.dds_ch0.set_cfr2(matched_latency_enable=1)
            self.dds_ch1.set_cfr2(matched_latency_enable=1)


    @kernel(flags={"fast-math"})
    def run(self):
        time_start_mu = now_mu() & ~0x7

        # enable output switches and set initial profile
        at_mu(time_start_mu)
        self._set_cpld_profile_switch_on()

        # start output
        at_mu(time_start_mu + self.time_system_prepare_delay_mu)
        self.dds_ch0.cpld.set_profile(1)
        # cancel output
        # at_mu(time_start_mu + self.time_system_prepare_delay_mu + self.time_delay_mu)
        # self.dds_ch1.cpld.set_profile(1)

        # tmp remove
        at_mu(time_start_mu + self.time_system_prepare_delay_mu + 475)
        self.ttl8.on()
        # at_mu(time_start_mu + self.time_system_prepare_delay_mu + self.time_delay_mu + 475)
        # tmp remove
        # self.ttl9.on()

        # cleanup
        at_mu(time_start_mu + self.time_system_prepare_delay_mu + self.time_delay_mu + self.time_system_cleanup_delay_mu)
        self._set_cpld_profile_switch_off()

    @kernel(flags={"fast-math"})
    def configure(self, freq_ftw: TInt32, phase_pow: TInt32, time_delay_mu: TInt64):
        # store timings
        time_delay_tmp_mu =                 time_delay_mu & ~0x7
        self.time_delay_mu =                time_delay_mu

        # calculate phase values to compensate for ch0 & ch1 delays
        self.phase_ch1_latency_turns =      self.dds_ch0.ftw_to_frequency(freq_ftw) * (self.time_latency_ch1_system_ns * ns)
        self.phase_ch1_delay_turns =        self.dds_ch0.ftw_to_frequency(freq_ftw) * (time_delay_tmp_mu * ns)

        # combine compensation values into ch1 phase: ch0 vs ch1 inherent phase shift/relation,
        # ch0 vs ch1 inherent time delay, and 0.5 (for cancellation)
        self.phase_ch1_final_pow =          self.dds_ch0.turns_to_pow(self.phase_ch1_inherent_turns +
                                                                      self.phase_ch1_latency_turns +
                                                                      self.phase_ch1_delay_turns +
                                                                      0.5) + phase_pow

        # set waveforms for profiles
        at_mu(now_mu() + 50000)
        with parallel:
            with sequential:
                self.dds_ch0.set_mu(freq_ftw, asf=0x01, pow_=phase_pow, profile=0)
                self.dds_ch0.set_mu(freq_ftw, asf=self.ampl_ticklefast_asf, pow_=phase_pow, profile=1)
                self.dds_ch0.write32(_AD9910_REG_CFR1, (1 << 16) | (1 << 13))
            with sequential:
                self.dds_ch1.set_mu(freq_ftw, asf=0x01, pow_=self.phase_ch1_final_pow, profile=0)
                self.dds_ch1.set_mu(freq_ftw, asf=self.ampl_ticklefast_asf, pow_=self.phase_ch1_final_pow, profile=1)
                self.dds_ch1.write32(_AD9910_REG_CFR1, (1 << 16) | (1 << 13))
        self.core.break_realtime()

    @kernel
    def _set_cpld_profile_switch_on(self):
        with parallel:
            with sequential:
                cfg0 = ((self.dds_ch0.cpld.cfg_reg) & ~(7 << CFG_PROFILE)) | ((1 << 3) | (0 << CFG_PROFILE))
                self.dds_ch0.cpld.cfg_write(cfg0)

            with sequential:
                cfg1 = ((self.dds_ch1.cpld.cfg_reg) & ~(7 << CFG_PROFILE)) | ((1 << 3) | (0 << CFG_PROFILE))
                self.dds_ch1.cpld.cfg_write(cfg1)

    @kernel
    def _set_cpld_profile_switch_off(self):
        with parallel:
            with sequential:
                cfg0 = ((self.dds_ch0.cpld.cfg_reg) & ~((7 << CFG_PROFILE) | (1 << 3))) | (0 << CFG_PROFILE)
                self.dds_ch0.cpld.cfg_write(cfg0)

            with sequential:
                cfg1 = ((self.dds_ch1.cpld.cfg_reg) & ~((7 << CFG_PROFILE) | (1 << 3))) | (0 << CFG_PROFILE)
                self.dds_ch1.cpld.cfg_write(cfg1)

            # tmp remove
            with sequential:
                delay_mu(100)
                self.ttl8.off()
                self.ttl9.off()
            # tmp remove

    def analyze(self):
        print('\n\tresults:')
        print('\t\turukul0 profile 1 phase: {:.4f}'.format(self.dds_ch0.pow_to_turns(self.phase_ch1_final_pow)))
        print('\n')
