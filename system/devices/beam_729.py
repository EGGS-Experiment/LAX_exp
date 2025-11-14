from artiq.experiment import *
from artiq.coredevice import ad9910
from artiq.coredevice.urukul import DEFAULT_PROFILE

from LAX_exp.extensions import *
from LAX_exp.base import LAXDevice


class Beam729(LAXDevice):
    """
    Device: Qubit Beam (729nm, polarized)

    Uses the DDS channel to drive the 729nm AOM in double-pass configuration.
    """
    name = "qubit"
    core_device = ('beam', 'urukul0_ch0')
    devices = {
        'rf_switch':        'ttl14',
        'singlepass0':      'urukul0_ch1',
        'singlepass1':      'urukul0_ch2',
        'singlepass2':      'urukul0_ch3',
        'doublepass_inj':   'urukul1_ch3',
    }

    kernel_invariants = {
        # hardware devices
        "cpld", "sw",

        # DDS parameters - qubit/main/doublepass (chamber)
        "freq_qubit_ftw", "ampl_qubit_asf",

        # DDS parameters - singlepass0 (injection lock)
        "freq_singlepass0_default_ftw", "ampl_singlepass0_default_asf", "att_singlepass0_default_mu",

        # DDS parameters - singlepass1 (injection lock)
        "freq_singlepass1_default_ftw", "ampl_singlepass1_default_asf", "att_singlepass1_default_mu",

        # DDS parameters - doublepass (injection lock)
        "freq_doublepass_inj_default_ftw", "ampl_doublepass_inj_default_asf",
        "att_doublepass_inj_default_mu",
    }

    def prepare_device(self):
        """
        todo: document
        """
        # re-alias relevant base devices
        self.sw =   self.beam.sw
        self.cpld = self.beam.cpld

        # get main DDS (chamber doublepass) parameters
        self.freq_qubit_ftw = self.get_parameter('freq_qubit_mhz', group='beams.freq_mhz', override=False,
                                                 conversion_function=hz_to_ftw, units=MHz)
        self.ampl_qubit_asf = self.get_parameter('ampl_qubit_pct', group='beams.ampl_pct', override=False,
                                                 conversion_function=pct_to_asf)

        # get singlepass (injection lock) parameters
        self.freq_singlepass0_default_ftw = self.get_parameter('freq_729_singlepass0_mhz', group='beams.freq_mhz',
                                                               override=False,
                                                               conversion_function=hz_to_ftw, units=MHz)
        self.freq_singlepass1_default_ftw = self.get_parameter('freq_729_singlepass1_mhz', group='beams.freq_mhz',
                                                               override=False,
                                                               conversion_function=hz_to_ftw, units=MHz)
        self.freq_singlepass2_default_ftw = self.get_parameter('freq_729_singlepass2_mhz', group='beams.freq_mhz',
                                                               override=False,
                                                               conversion_function=hz_to_ftw, units=MHz)

        self.ampl_singlepass0_default_asf = self.get_parameter('ampl_729_singlepass0_pct', group='beams.ampl_pct',
                                                               override=False, conversion_function=pct_to_asf)
        self.ampl_singlepass1_default_asf = self.get_parameter('ampl_729_singlepass1_pct', group='beams.ampl_pct',
                                                               override=False, conversion_function=pct_to_asf)
        self.ampl_singlepass2_default_asf = self.get_parameter('ampl_729_singlepass2_pct', group='beams.ampl_pct',
                                                               override=False, conversion_function=pct_to_asf)

        self.att_singlepass0_default_mu = self.get_parameter('att_729_singlepass0_db', group='beams.att_db',
                                                             override=False, conversion_function=att_to_mu)
        self.att_singlepass1_default_mu = self.get_parameter('att_729_singlepass1_db', group='beams.att_db',
                                                             override=False, conversion_function=att_to_mu)
        self.att_singlepass2_default_mu = self.get_parameter('att_729_singlepass2_db', group='beams.att_db',
                                                             override=False, conversion_function=att_to_mu)

        # get doublepass (injection lock) parameters
        self.freq_doublepass_inj_default_ftw =  self.get_parameter('freq_729_doublepass_inj_mhz',
                                                                   group='beams.freq_mhz', override=False,
                                                                   conversion_function=hz_to_ftw, units=MHz)
        self.ampl_doublepass_inj_default_asf =  self.get_parameter('ampl_729_doublepass_inj_pct',
                                                                   group='beams.ampl_pct', override=False,
                                                                   conversion_function=pct_to_asf)
        self.att_doublepass_inj_default_mu =    self.get_parameter('att_729_doublepass_inj_db',
                                                                   group='beams.att_db', override=False,
                                                                   conversion_function=att_to_mu)

    def _prepare_device_values(self):
        """
        Check device parameters for validity.
        """
        # todo: check all ampls in [0, 50]
        # todo: check singlepass atts in [7, 31.5]
        # todo: check doublepass att in [6, 31.5]
        pass


    '''
    LAXDEVICE METHODS
    '''
    @kernel(flags={"fast-math"})
    def initialize_device(self) -> TNone:
        """
        Initialize all devices relevant to the 729nm.
        """
        # get CPLD attenuations so we don't override them
        self.cpld.get_att_mu()
        self.doublepass_inj.cpld.get_att_mu()
        self.core.break_realtime()

        # ensure phase_autoclear disabled on all beams to prevent phase accumulator reset
        self.set_cfr1()
        self.singlepass0.set_cfr1()
        self.singlepass1.set_cfr1()
        self.singlepass2.set_cfr1()
        self.doublepass_inj.set_cfr1()
        self.io_update()
        delay_mu(25000)

        # set matched_latency_enable on all relevant DDSs for consistency
        self.set_cfr2(matched_latency_enable=1)
        self.singlepass0.set_cfr2(matched_latency_enable=1)
        self.singlepass1.set_cfr2(matched_latency_enable=1)
        self.singlepass2.set_cfr2(matched_latency_enable=1)
        self.doublepass_inj.set_cfr2(matched_latency_enable=1)
        self.io_update()
        delay_mu(25000)

        # set up relevant AOMs to default values on ALL profiles
        # necessary b/c not all AOMs are configured/used for each experiment
        for i in range(8):
            self.singlepass0.set_mu(self.freq_singlepass0_default_ftw,
                                    asf=self.ampl_singlepass0_default_asf,
                                    profile=i, phase_mode=ad9910.PHASE_MODE_CONTINUOUS)
            delay_mu(25000) # 25 us
            self.singlepass1.set_mu(self.freq_singlepass1_default_ftw,
                                    asf=self.ampl_singlepass1_default_asf,
                                    profile=i, phase_mode=ad9910.PHASE_MODE_CONTINUOUS)
            delay_mu(25000) # 25 us
            self.singlepass2.set_mu(self.freq_singlepass2_default_ftw,
                                    asf=self.ampl_singlepass2_default_asf,
                                    profile=i, phase_mode=ad9910.PHASE_MODE_CONTINUOUS)
            delay_mu(25000) # 25 us
            self.doublepass_inj.set_mu(self.freq_doublepass_inj_default_ftw,
                                    asf=self.ampl_doublepass_inj_default_asf,
                                    profile=i, phase_mode=ad9910.PHASE_MODE_CONTINUOUS)
            delay_mu(25000) # 25 us

        # ensure events finish completion (since they're pretty heavy tbh)
        self.core.break_realtime()
        self.core.wait_until_mu(now_mu())
        self.core.break_realtime()
        delay_mu(500000) # 500 us

        # set AOMs for normal output/operation
        self.singlepass0.set_att_mu(self.att_singlepass0_default_mu)
        self.singlepass1.set_att_mu(self.att_singlepass1_default_mu)
        self.singlepass2.set_att_mu(self.att_singlepass2_default_mu)
        self.doublepass_inj.set_att_mu(self.att_doublepass_inj_default_mu)
        delay_mu(25000)
        self.singlepass0.sw.on()
        self.singlepass1.sw.off()
        self.singlepass2.sw.off()
        self.doublepass_inj.sw.on()
        delay_mu(25000)

    @kernel(flags={"fast-math"})
    def cleanup_device(self) -> TNone:
        """
        Restore all devices relevant to the 729nm back for normal operation.
        """
        # set up relevant AOMs to default values on ALL profiles
        # necessary b/c not all AOMs are configured/used for each experiment
        for i in range(8):
            self.singlepass0.set_mu(self.freq_singlepass0_default_ftw,
                                    asf=self.ampl_singlepass0_default_asf,
                                    profile=i, phase_mode=ad9910.PHASE_MODE_CONTINUOUS)
            delay_mu(50000) # 50 us
            self.singlepass1.set_mu(self.freq_singlepass1_default_ftw,
                                    asf=self.ampl_singlepass1_default_asf,
                                    profile=i, phase_mode=ad9910.PHASE_MODE_CONTINUOUS)
            delay_mu(50000) # 50 us
            self.doublepass_inj.set_mu(self.freq_doublepass_inj_default_ftw,
                                    asf=self.ampl_doublepass_inj_default_asf,
                                    profile=i, phase_mode=ad9910.PHASE_MODE_CONTINUOUS)
            delay_mu(50000) # 50 us
        self.core.break_realtime()

        # set AOMs for normal output/operation
        self.singlepass0.set_att_mu(self.att_singlepass0_default_mu)
        self.singlepass1.set_att_mu(self.att_singlepass1_default_mu)
        self.doublepass_inj.set_att_mu(self.att_doublepass_inj_default_mu)
        delay_mu(25000)

        self.singlepass0.sw.on()
        delay_mu(8)
        self.singlepass1.sw.off()
        delay_mu(8)
        self.doublepass_inj.sw.on()
        delay_mu(25000)

        # return to default profile on CPLD (this is the default profile used by user/GUIs)
        self.set_profile(DEFAULT_PROFILE)
        self.io_update()
        self.sw.off() # note: only turn off DDS int sw - leave ext sw OK for user
        delay_mu(10000)


    '''
    HARDWARE METHODS
    '''
    @kernel(flags={"fast-math"})
    def on(self) -> TNone:
        """
        todo: document
        """
        with parallel:
            # enable RF switch onboard Urukul
            self.sw.on()

            # enable external RF switch
            with sequential:
                self.rf_switch.off()
                delay_mu(TIME_ZASWA2_SWITCH_DELAY_MU)

    @kernel(flags={"fast-math"})
    def off(self) -> TNone:
        """
        todo: document
        """
        with parallel:
            # disable RF switch onboard Urukul
            self.sw.off()

            # disable external RF switch
            with sequential:
                self.rf_switch.on()
                delay_mu(TIME_ZASWA2_SWITCH_DELAY_MU)

    @kernel(flags={"fast-math"})
    def set_profile(self, profile_num: TInt32) -> TNone:
        """
        todo: document
        todo: do we really need this delay? shouldn't we hit io_update THEN delay?
        maybe: we don't really want set_profile to pulse io_update b/c it might
            fuck things up, e.g. clear phase accumulator
        :param profile_num: the AD9910 profile number to set. Must be an int in [0, 7].
        """
        self.cpld.set_profile(profile_num)
        delay_mu(TIME_AD9910_PROFILE_SWITCH_DELAY_MU)

    @kernel(flags={"fast-math"})
    def io_update(self) -> TNone:
        """
        Pulse the CPLDs IO_UPDATE pin.
        Can be used to clear the phase accumulator if the phase_autoclear
            flag is set in CFR1.
        """
        self.cpld.io_update.pulse_mu(8)
