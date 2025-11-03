from artiq.experiment import *
from artiq.coredevice import ad9910
from artiq.coredevice import urukul

# todo: enable empty dds_profiles (make -1 and check)


class UrukulConfigure(EnvExperiment):
    """
    Tool: Urukul Configure

    Configure values for an Urukul channel (AD9910 ONLY).
    Urukul is set to the default profile (currently profile #7) afterwards.
    """
    name = 'Urukul Configure'
    kernel_invariants = {
        # hardware & parameters
        "dds"
    }

    def build(self):
        self.setattr_device("core")

        # DDS parameters
        dds_device_list = self._get_dds_devices()
        self.setattr_argument("dds_target", EnumerationValue(list(dds_device_list), default='urukul0_ch1'),
                              tooltip="DDS to configure.")

        # initialization
        self.setattr_argument("initialize_cpld",    BooleanValue(default=False), group="initialize",
                              tooltip="Initialize the parent urukul's CPLD.")
        self.setattr_argument("master_reset",  BooleanValue(default=False), group="initialize",
                              tooltip="Pulse the MASTER_RESET pin on the CPLD for ALL DDS outputs controlled by the CPLD.\n"
                                      "Clears all memory and resets registers to defaults for ALL AD9910s controlled by the Urukul CPLD.\n")
        self.setattr_argument("initialize_ad9910",  BooleanValue(default=False), group="initialize",
                              tooltip="Initialize the AD9910. Relocks the PLL, so relative DDS phases can be shuffled.")

        # CPLD parameters
        self.setattr_argument('dds_profiles',   PYONValue([0, 1, 2, 3, 4, 5, 6, 7]), group="config",
                              tooltip="DDS profiles to configure with the selected parameters.\n"
                                      "Can be an empty list.")
        self.setattr_argument("switch_on",      BooleanValue(default=True), group="config",
                              tooltip="Leave the channel switch ON/OFF after initialization.")

        # DDS (AD9910) parameters
        self.setattr_argument('att_db',     NumberValue(default=7., precision=1, step=0.5, min=0., max=31.5, unit="dB", scale=1.),
                              group="waveform",
                              tooltip="Attenuation (in dB) to set for the DDS channel.")
        self.setattr_argument('freq_mhz',   NumberValue(default=120.339, precision=7, step=5, min=0., max=400., unit="MHz", scale=1.),
                              group="waveform",
                              tooltip="Frequency to set for the DDS channel.")
        self.setattr_argument('ampl_pct',   NumberValue(default=50., precision=3, step=5., min=0.01, max=100., unit="%", scale=1.),
                              group="waveform",
                              tooltip="Amplitude (in %) to set for the DDS channel.")

    def prepare(self):
        """
        Precompute relevant values.
        """
        # get target DDS device
        try:
            self.dds = self.get_device(self.dds_target)
        except Exception as e:
            raise e

        # sanitize inputs
        if not isinstance(self.dds_profiles, list):
            raise ValueError("Invalid argument: dds_profiles must be a list.")
        # ensure dds_profiles are int in [0, 7]
        dds_profiles_valid = [isinstance(val, int) and (0 <= val <= 7) for val in self.dds_profiles]
        if not all(dds_profiles_valid):
            raise ValueError("Invalid DDS profile list. Must be an int in [0, 7].")
        # in event of empty dds_profiles list, set to -1 (so we can ignore it later on)
        if len(self.dds_profiles) == 0:
            self.dds_profiles = [-1]

    def _get_dds_devices(self):
        """
        Get all valid DDS (AD9910) devices from the device_db.
        :return: a set of all AD9910 devices.
        """
        is_local_dds_device = lambda v: (
                isinstance(v, dict) and (v.get('type') == 'local')
                and ('class' in v) and (v.get('class') == "AD9910")
        )

        # return sorted list of local DDS devices from device_db
        return sorted(set([
            k
            for k, v in self.get_device_db().items()
            if is_local_dds_device(v)
        ]))


    '''
    MAIN SEQUENCE
    '''
    @kernel(flags={"fast-math"})
    def run(self) -> TNone:
        # reset & prepare
        self.core.wait_until_mu(now_mu())
        self.core.break_realtime()
        self.core.reset()


        '''PREPARE'''
        # initialize relevant devices
        # note: allow events to complete before continuing
        if self.initialize_cpld:
            self.dds.cpld.init()
            self.core.break_realtime()
            self.core.wait_until_mu(now_mu())
            self.core.break_realtime()
            delay_mu(1000000) # 1ms

        if self.master_reset:
            self.dds.cpld.cfg_write(self.dds.cpld.cfg_reg | (1 << urukul.CFG_RST))
            delay_mu(100000) # 100 us
            self.dds.cpld.cfg_write(self.dds.cpld.cfg_reg & ~(1 << urukul.CFG_RST))

        if self.initialize_ad9910:
            self.dds.init()
            self.core.break_realtime()
            self.core.wait_until_mu(now_mu())
            self.core.break_realtime()
            delay_mu(1000000) # 1ms

        # get board attenuations to prevent unwanted overwrites
        self.dds.cpld.get_att_mu()
        self.core.break_realtime()

        # close DDS switch and set max att (prevent leakage while we update)
        self.dds.sw.off()
        self.dds.set_att(31.5 * dB)

        # clear CFRs to default state
        self.dds.set_cfr1()
        self.dds.cpld.io_update.pulse_mu(8)
        self.dds.set_cfr2()
        self.dds.cpld.io_update.pulse_mu(8)
        delay_mu(10000)


        '''SETUP'''
        # set parameters for DDS profiles
        for profile_num in self.dds_profiles:
            if profile_num != -1:
                self.dds.set(self.freq_mhz * MHz, amplitude=self.ampl_pct / 100.,
                             phase_mode=ad9910.PHASE_MODE_CONTINUOUS,
                             profile=profile_num)
                delay_mu(50000)


        '''CLEAN UP'''
        # set DDS to profile 7 (i.e. the default profile)
        self.dds.cpld.set_profile(urukul.DEFAULT_PROFILE)
        self.dds.cpld.io_update.pulse_mu(8)

        # set DDS output for use
        self.dds.set_att(self.att_db * dB)
        if self.switch_on:  self.dds.sw.on()
        else:   self.dds.sw.off()

        # ensure events complete submission
        self.core.wait_until_mu(now_mu())
        self.core.break_realtime()

