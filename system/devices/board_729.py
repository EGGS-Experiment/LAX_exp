from artiq.experiment import *
from artiq.coredevice.urukul import DEFAULT_PROFILE

from LAX_exp.extensions import *
from LAX_exp.base import LAXDevice

class GenericBeam729(LAXDevice):

    """
        Device: Generic DDS for 729 Beam (729nm, polarized)
    """

    core_device = None
    devices = None
    core_device = None


    def prepare_device(self):
        # re-alias relevant base devices
        self.sw = self.beam.sw
        self.cpld = self.beam.cpld

        # get beam parameters
        self.freq_qubit_ftw = self.get_parameter('freq_qubit_mhz', group='beams.freq_mhz', override=False,
                                                 conversion_function=hz_to_ftw, units=MHz)
        self.ampl_qubit_asf = self.get_parameter('ampl_qubit_pct', group='beams.ampl_pct', override=False,
                                                 conversion_function=pct_to_asf)

    @kernel(flags={"fast-math"})
    def initialize_device(self) -> TNone:
        # get CPLD attenuations so we don't override them
        self.cpld.get_att_mu()
        self.core.break_realtime()

    @kernel(flags={"fast-math"})
    def cleanup_device(self) -> TNone:
        # set default profile on CPLD
        self.core.break_realtime()
        self.set_profile(DEFAULT_PROFILE)
        self.cpld.io_update.pulse_mu(8)
        self.off()

    @kernel(flags={"fast-math"})
    def on(self) -> TNone:
        with parallel:
            # enable RF switch onboard Urukul
            self.sw.on()

            # enable external RF switch
            with sequential:
                self.rf_switch.off()
                delay_mu(TIME_ZASWA2_SWITCH_DELAY_MU)

    @kernel(flags={"fast-math"})
    def off(self) -> TNone:
        with parallel:
            # disable RF switch onboard Urukul
            self.sw.off()

            # disable external RF switch
            with sequential:
                self.rf_switch.on()
                delay_mu(TIME_ZASWA2_SWITCH_DELAY_MU)

    @kernel(flags={"fast-math"})
    def set_profile(self, profile_num: TInt32) -> TNone:
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


class Beam729System(LAXDevice):
    name = "filter_cavity_qubit"
    core_device = None

    def build(self, **kwargs):
        super(Beam729System, self).build()
        self.setattr_device("doublepass")
        self.setattr_device("singlepass1")
        self.setattr_device("singlepass2")


class DoublePass729(GenericBeam729):
    name = "doublepass"
    # tmp!!!!!
    core_device = ('beam', 'urukul1_ch3')
    devices = {
            'rf_switch': 'ttl14'
    }
    # tmp!!!!
    kernel_invariants = {
        "cpld", "sw",
            "freq_qubit_ftw", "ampl_qubit_asf"
    }



class FirstSinglePass729(GenericBeam729):
    name = "singlepass1"
    # tmp!!!!!
    core_device = ('beam', 'urukul1_ch2')
    devices = {
      'rf_switch': 'ttl14'
    }
    kernel_invariants = {
        "cpld", "sw",
        "freq_qubit_ftw", "ampl_qubit_asf"
    }

class SecondSinglePass729(GenericBeam729):
    name = "singlepass1"
    # tmp!!!!!
    core_device = ('beam', 'urukul1_ch1')
      devices = {
        'rf_switch': 'ttl14'
    }
    kernel_invariants = {
        "cpld", "sw",
        "freq_qubit_ftw", "ampl_qubit_asf"
    }




