from artiq.experiment import *
from artiq.coredevice import urukul


class UrukulConfigure(EnvExperiment):
    """
    Tool: Urukul Configure

    Configure values for an Urukul channel (AD9910 ONLY).
    Urukul is set to the default profile (currently profile #7) afterwards.
    """
    name = 'Urukul Configure'
    kernel_invariants = {
        # hardware & parameters
        "dds", "freq_ftw", "ampl_asf"
    }

    def build(self):
        self.setattr_device("core")

        # DDS parameters
        dds_device_list = self._get_dds_devices()
        self.setattr_argument("dds_target", EnumerationValue(list(dds_device_list), default='urukul0_ch1'))

        # initialization
        self.setattr_argument("initialize_cpld",    BooleanValue(default=False), group="initialize",
                              tooltip="Initialize the parent urukul's CPLD.")
        self.setattr_argument("initialize_ad9910",  BooleanValue(default=False), group="initialize",
                              tooltip="Initialize the AD9910. Will clear all saved values on the AD9910.")

        # CPLD parameters
        self.setattr_argument('dds_profiles',   PYONValue([0, 1, 2, 3, 4, 5, 6, 7]), group="config",
                              tooltip="DDS profiles to configure with the selected parameters.")
        self.setattr_argument("switch_on",      BooleanValue(default=True), group="config",
                              tooltip="Leave the channel switch ON/OFF after initialization.")

        # DDS (AD9910) parameters
        self.setattr_argument('att_db',     NumberValue(default=7., precision=1, step=0.5, min=0., max=31.5, unit="dB", scale=1.), group="waveform")
        self.setattr_argument('freq_mhz',   NumberValue(default=120.339, precision=7, step=5, min=0., max=400., unit="MHz", scale=1.), group="waveform")
        self.setattr_argument('ampl_pct',   NumberValue(default=50., precision=3, step=5., min=0.01, max=100., unit="%", scale=1.), group="waveform")

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
        dds_profiles_valid = [isinstance(val, int) and (0 <= val <= 7) for val in self.dds_profiles]
        if not all(dds_profiles_valid):
            raise ValueError("Invalid DDS profile list. Must be an int in [0, 7].")

        # convert DDS parameters to machine units
        self.freq_ftw = self.dds.frequency_to_ftw(self.freq_mhz * MHz)
        self.ampl_asf = self.dds.amplitude_to_asf(self.ampl_pct / 100.)

    def _get_dds_devices(self):
        """
        Get all valid DDS (AD9910) devices from the device_db.
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

    @kernel(flags={"fast-math"})
    def run(self) -> TNone:
        # reset & prepare
        self.core.wait_until_mu(now_mu())
        self.core.break_realtime()
        self.core.reset()

        '''PREPARE'''
        # todo: run initialization etc.
        if self.initialize_cpld:
            self.dds.cpld.init()
        self.core.break_realtime()
        delay_mu(1000000) # 1ms

        if self.initialize_ad9910:
            self.dds.init()
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
            self.dds.set_mu(self.freq_ftw, asf=self.ampl_asf, profile=profile_num)
            delay_mu(10000)

        '''CLEAN UP'''
        # set DDS to profile 7 (i.e. the default profile)
        self.dds.cpld.set_profile(urukul.DEFAULT_PROFILE)
        self.dds.cpld.io_update.pulse_mu(8)

        # set DDS output for use
        self.dds.set_att(self.att_db * dB)
        if self.switch_on:
            self.dds.sw.on()
        else:
            self.dds.sw.off()

        # clean up
        self.core.wait_until_mu(now_mu())
        self.core.reset()

