from artiq.experiment import *
from artiq.coredevice import ad9910
from artiq.coredevice.urukul import DEFAULT_PROFILE

from LAX_exp.extensions import *
from LAX_exp.base import LAXDevice
from numpy import int64


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

        # DDS parameters - singlepass2 (injection lock)
        "freq_singlepass2_default_ftw", "ampl_singlepass2_default_asf", "att_singlepass2_default_mu",

        # DDS parameters - doublepass (injection lock)
        "freq_doublepass_inj_default_ftw", "ampl_doublepass_inj_default_asf",
        "att_doublepass_inj_default_mu",

        # switch delay time
        "switch_delay_time_mu",

        # note: ensure objects used for programmatic initialization (e.g. device_list, freq_ftw_list) are
        #   made kernel_invariant
    }

    def prepare_device(self):
        """
        todo: document
        """
        # re-alias relevant base devices
        self.sw = self.beam.sw
        self.cpld = self.beam.cpld

        self.dds_devices = ['singlepass0', 'singlepass1', 'singlepass2', 'doublepass_inj']

        # get main DDS (chamber doublepass) parameters
        self.freq_qubit_ftw = self.get_parameter('freq_qubit_mhz', group='beams.freq_mhz', override=False,
                                                 conversion_function=hz_to_ftw, units=MHz)
        self.ampl_qubit_asf = self.get_parameter('ampl_qubit_pct', group='beams.ampl_pct', override=False,
                                                 conversion_function=pct_to_asf)

        for device in self.dds_devices:
            setattr(self, f'freq_{device}_default_ftw',
                    self.get_parameter(f'freq_729_{device}_mhz', group='beams.freq_mhz',
                                       override=False,
                                       conversion_function=hz_to_ftw, units=MHz))

            setattr(self, f'ampl_{device}_default_asf',
                    self.get_parameter(f'ampl_729_{device}_pct', group='beams.ampl_pct',
                                       override=False, conversion_function=pct_to_asf))

            setattr(self, f'att_{device}_default_mu', self.get_parameter(f'att_729_{device}_db', group='beams.att_db',
                                                                         override=False, conversion_function=att_to_mu))

        self.device_list =  []
        self.freq_ftw_list = []
        self.ampl_asf_list =  []
        self.att_mu_list =  []
        for idx, device in enumerate(self.dds_devices):
            self.device_list.append(getattr(self, device))
            self.freq_ftw_list.append(getattr(self, f'freq_{device}_default_ftw'))
            self.ampl_asf_list.append(getattr(self, f'ampl_{device}_default_asf'))
            self.att_mu_list.append(getattr(self, f'att_{device}_default_mu'))

        # note: make all delays self.core.coarse_ref_period instead of 8ns
        self.switch_delay_time_mu = int64(8)

    def _check_device_values(self):
        """
        Check device parameters for validity.
        """
        # todo: check all ampls in [0, 50]
        # todo: check singlepass atts in [7, 31.5]
        # todo: check doublepass att in [6, 31.5]

        # clayton note: move "magic numbers" to fixed declarations in build_ for transparency
        #   and document clearly
        device_att_db = {'singlepass0': 7 * dB, 'singlepass1': 7 * dB, 'singlepass2': 7 * dB, 'doublepass_inj': 6 * dB}

        for device in self.dds_devices:
            device_attr = getattr(self, device)
            freq_mhz = device_attr.ftw_to_frequency(getattr(self, f'freq_{device}_default_ftw')) / MHz
            ampl_pct = device_attr.asf_to_amplitude(getattr(self, f'ampl_{device}_default_asf')) * 100.
            att_db = mu_to_att(getattr(self, f'att_{device}_default_mu'))

            if not 0 < ampl_pct < 50:
                raise ValueError(f"{device} must have amplitude between 0 and 50 percent")
            if not device_att_db[device] <= att_db <= 31.5:
                raise ValueError(f"{device} must have attenuation between {device_att_db[device]} and 31.5dB")

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

        self.set_cfr1()
        # ensure phase_autoclear disabled on all beams to prevent phase accumulator reset
        for idx in range(len(self.device_list)):
            device_attr = self.device_list[idx]
            device_attr.set_cfr1()
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
        for idx in range(len(self.device_list)):
            device_attr = self.device_list[idx]
            freq_attr = self.freq_ftw_list[idx]
            ampl_attr = self.ampl_asf_list[idx]
            att_attr = self.att_mu_list[idx]
            device_attr.set_att_mu(att_attr)
            delay_mu(25000)
            for i in range(8):
                device_attr.set_mu(freq_attr, asf=ampl_attr,
                                   profile=i, phase_mode=ad9910.PHASE_MODE_CONTINUOUS)
                delay_mu(25000)

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
        delay_mu(8)
        self.singlepass1.sw.off()
        delay_mu(8)
        self.singlepass2.sw.off()
        delay_mu(8)
        self.doublepass_inj.sw.on()
        delay_mu(25000)

    @kernel(flags={"fast-math"})
    def cleanup_device(self) -> TNone:
        """
        Restore all devices relevant to the 729nm back for normal operation.
        """
        # set up relevant AOMs to default values on ALL profiles
        # necessary b/c not all AOMs are configured/used for each experiment
        for idx in range(len(self.device_list)):
            device_attr = self.device_list[idx]
            freq_attr = self.freq_ftw_list[idx]
            ampl_attr = self.ampl_asf_list[idx]
            att_attr = self.att_mu_list[idx]
            device_attr.set_att_mu(att_attr)
            delay_mu(25000)
            for i in range(8):
                device_attr.set_mu(freq_attr, asf=ampl_attr, profile=i, phase_mode=ad9910.PHASE_MODE_CONTINUOUS)
                delay_mu(8000)

        self.singlepass0.sw.on()
        delay_mu(8)
        self.singlepass1.sw.off()
        delay_mu(8)
        self.singlepass2.sw.off()
        self.doublepass_inj.sw.on()
        delay_mu(25000)

        # return to default profile on CPLD (this is the default profile used by user/GUIs)
        self.set_profile(DEFAULT_PROFILE)
        self.io_update()
        self.sw.off()  # note: only turn off DDS int sw - leave ext sw OK for user
        delay_mu(10000)


    '''
    HARDWARE METHODS
    '''

    @kernel(flags={"fast-math"})
    def on(self) -> TNone:
        """
        todo: document
        """

        self.sw.on()
        delay_mu(8)
        # enable external RF switch
        self.rf_switch.off()
        delay_mu(TIME_ZASWA2_SWITCH_DELAY_MU)

    @kernel(flags={"fast-math"})
    def off(self) -> TNone:
        """
        todo: document
        """
        # disable RF switch onboard Urukul
        self.sw.off()
        delay_mu(8)
        # disable external RF switch
        self.rf_switch.on()
        delay_mu(TIME_ZASWA2_SWITCH_DELAY_MU)

    @kernel(flags={'fast-math'})
    def singlepass0_on(self) -> TNone:
        self.singlepass0.sw.on()
        delay_mu(self.switch_delay_time_mu)

    @kernel(flags={'fast-math'})
    def singlepass0_off(self) -> TNone:
        self.singlepass0.sw.off()
        delay_mu(self.switch_delay_time_mu)

    @kernel(flags={'fast-math'})
    def singlepass1_on(self) -> TNone:
        self.singlepass1.sw.on()
        delay_mu(self.switch_delay_time_mu)

    @kernel(flags={'fast-math'})
    def singlepass1_off(self) -> TNone:
        self.singlepass1.sw.off()
        delay_mu(self.switch_delay_time_mu)

    @kernel(flags={'fast-math'})
    def singlepass2_on(self) -> TNone:
        self.singlepass2.sw.on()
        delay_mu(self.switch_delay_time_mu)

    @kernel(flags={'fast-math'})
    def singlepass2_off(self) -> TNone:
        self.singlepass2.sw.off()
        delay_mu(self.switch_delay_time_mu)

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
