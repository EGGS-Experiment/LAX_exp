from artiq.experiment import *
from artiq.coredevice import ad9910
from numpy import int32, int64

from enum import Enum
# todo: make usable for multiple profiles a la urukul_configure


class RAM_MODE(Enum):
    """
    Enum class to convert the RAM mode bits to a string
    """
    RAM_MODE_DIRECTSWITCH =     ad9910.RAM_MODE_DIRECTSWITCH
    RAM_MODE_RAMPUP =           ad9910.RAM_MODE_RAMPUP
    RAM_MODE_BIDIR_RAMP =       ad9910.RAM_MODE_BIDIR_RAMP
    RAM_MODE_CONT_BIDIR_RAMP =  ad9910.RAM_MODE_CONT_BIDIR_RAMP
    RAM_MODE_CONT_RAMPUP =      ad9910.RAM_MODE_CONT_RAMPUP


class UrukulRead(EnvExperiment):
    """
    Tool: Urukul Read

    Read values from an Urukul channel.
    Warning: urukul profile will be set to target value and left there.
    """
    name = 'Urukul Read'

    def build(self):
        self.setattr_device("core")

        # DDS parameters
        dds_device_list = self._get_dds_devices()
        self.setattr_argument("dds_target",         EnumerationValue(list(dds_device_list), default='urukul2_ch1'))

        # profile parameters
        self.setattr_argument('dds_profile',        NumberValue(default=7, precision=0, step=1, min=0, max=7))
        self.setattr_argument("is_profile_ram",     BooleanValue(default=False))

        # RAM data
        self.setattr_argument("read_ram_data",      BooleanValue(default=False))

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

    def prepare(self):
        """
        Prepare values for speedy evaluation.
        """
        '''SANITIZE INPUTS'''
        # ensure profile number is valid, otherwise strange failure mode/dangerous overwrite
        if (type(self.dds_profile) is not int) or ((self.dds_profile < 0) | (self.dds_profile > 7)):
            raise Exception("Error: invalid profile number ({:d}). Must be integer in [0, 7].".format(self.dds_profile))

        # ensure we only read RAM data if profile is a RAM profile
        if (self.read_ram_data is True) & (self.is_profile_ram is False):
            raise Exception("Error: cannot grab RAM data from non-RAM profile.")


        '''GET DEVICES'''
        try:
            self.dds = self.get_device(self.dds_target)
            self.dds_cpld = self.dds.cpld

            # add DDS name to kernel invariants
            kernel_invariants = getattr(self, "kernel_invariants", set())
            self.kernel_invariants = kernel_invariants | {'dds', 'dds_cpld'}

        except Exception as e:
            print("Error: invalid DDS channel.")
            raise e


        '''PREPARE DATA STRUCTURES'''
        # set up data structure to get data for tone and ram profile
        self._profile_word =    int64(0)

        # set up dataset to store data for ram read
        self._ram_data_array =  [int32(0)] * 1024


    """
    MAIN FUNCTIONS
    """
    @kernel(flags={"fast-math"})
    def run_initialize(self) -> TNone:
        """
        Initialize devices prior to execution.
        """
        # reset core
        self.core.reset()
        self.core.break_realtime()

        # ensure all DDS outputs are OFF to prevent booboos
        self.dds.cpld.cfg_switches(0b0000)
        self.core.break_realtime()

        # store state of att register to prevent booboos
        self.dds.cpld.get_att_mu()
        self.core.break_realtime()

        # initialize urukul (idk why but everyone does it each time)
        self.dds.cpld.init()
        self.core.break_realtime()
        delay_mu(100000)

        # initialize DDS (idk why but everyone does it each time)
        # warning - this may ruin timing since it tunes sync_delay
        # self.dds.init()
        self.core.break_realtime()
        delay_mu(100000)

    @kernel(flags={"fast-math"})
    def run(self) -> TNone:
        # prepare devices for experiment
        self.run_initialize()

        # prepare DDS in desired profile (necessary since values are only taken from active profile)
        # todo: modularize - args are dds num (from list) and profile
        # todo: set up dds name list structure
        self.dds_cpld.set_profile(self.dds_profile)
        self.dds_cpld.io_update.pulse_mu(8)
        self.core.break_realtime()

        # todo: need to read directly from register (instead of using get_mu) since profile may be RAM
        self._profile_word = int64(self.dds.read64(ad9910._AD9910_REG_PROFILE0 + self.dds_profile))
        self.core.break_realtime()

        # # read data from RAM
        # if self.read_ram_data:
        #     self.read_ram(self._ram_data_array)
        #     delay_mu(10000000)
        #     self.core.break_realtime()

        # clean up by setting default profile
        self.dds_cpld.set_profile(ad9910.DEFAULT_PROFILE)
        self.dds_cpld.io_update.pulse_mu(8)
        self.core.wait_until_mu(now_mu())


    # @kernel(flags={"fast-math"})
    # def read_profile_data(self, dds_num: TInt32, profile_num: TInt32) -> TInt64:
    #     """
    #     todo: document
    #     """
    #     # get desired DDS and CPLD from device list
    #     dds_dev =   self.dds_device_list[dds_num]
    #     dds_cpld =  dds_dev.cpld
    #
    #     # prepare DDS in desired profile (necessary since values are only taken from active profile)
    #     dds_cpld.set_profile(profile_num)
    #     dds_cpld.io_update.pulse_mu(8)
    #     self.core.break_realtime()
    #
    #     # note: need to read directly from register (instead of using get_mu) since profile may be RAM
    #     _profile_word = int64(dds_dev.read64(ad9910._AD9910_REG_PROFILE0 + profile_num))
    #     self.core.break_realtime()
    #
    #     return _profile_word


    """
    DATA PROCESSING
    """
    def analyze(self):
        """
        Process results to extract DDS data.
        """
        # print & save results
        print("\tDDS:\t\t{}".format(self.dds_target))
        print("\tProfile:\t{}".format(self.dds_profile))

        # process profile data as single-tone waveform
        if not self.is_profile_ram:
            # extract values from profile word
            freq_mhz, ampl_pct, phas_turns = self.process_profile_word_singletone(self._profile_word)

            # print results
            print("\t\tFreq.: {:.4f} MHz".format(freq_mhz))
            print("\t\tAmpl.: {:.4f} %".format(ampl_pct))
            print("\t\tPhas.: {:.4f} turns\n".format(phas_turns))

            # store results in datasets
            self.set_dataset("freq_mhz", freq_mhz)
            self.set_dataset("ampl_pct", ampl_pct)
            self.set_dataset("phas_turns", phas_turns)

        # process profile data as RAM
        else:
            # extract values from profile word
            ram_mode_str, start_reg, stop_reg, step_interval_ns, nodwell_high, zero_crossing = self.process_profile_word_ram(self._profile_word)

            # print results
            print("\t\tRAM Mode:\t{:s}".format(ram_mode_str))
            print("\t\tStart reg.:\t{:d}".format(start_reg))
            print("\t\tStop reg.:\t{:d}".format(stop_reg))
            print("\t\tStep rate:\t{:.4f} ns".format(step_interval_ns))
            print("\t\tNo-dwell high:\t{:b}".format(nodwell_high))
            print("\t\tZero-crossing:\t{:b}".format(zero_crossing))

            # store results in datasets
            self.set_dataset("ram_mode", ram_mode_str)
            self.set_dataset("start_reg", start_reg)
            self.set_dataset("stop_reg", stop_reg)
            self.set_dataset("step_interval_ns", step_interval_ns)
            self.set_dataset("nodwell_high", nodwell_high)
            self.set_dataset("zero_crossing", zero_crossing)

            # todo: print and store RAM array data


    # @rpc
    def process_profile_word_singletone(self, profile_word: TInt64) -> TTuple([TFloat, TFloat, TFloat]):
        """
        Extract single-tone waveform from DDS profile word.
        """
        # extract values from profile word
        ftw = int32(profile_word & 0xFFFFFFFF)
        pow = int32((profile_word >> 32) & 0xFFFF)
        asf = int32((profile_word >> 48) & 0x3FFF)

        # convert values to human units
        freq_mhz =      self.dds.ftw_to_frequency(ftw) / MHz
        ampl_pct =      self.dds.asf_to_amplitude(asf) * 100.
        phas_turns =    self.dds.pow_to_turns(pow)

        # return single-tone waveform values
        return freq_mhz, ampl_pct, phas_turns

    @rpc
    def process_profile_word_ram(self, profile_word: TInt64) -> TTuple([TStr, TInt32, TInt32, TFloat, TBool, TBool]):
        """
        Extract RAM waveform from DDS profile word.
        """
        # extract values from profile word
        ram_mode_val =      int32(profile_word & (0xFFFF << 40))       # 16b is 0xFFFF
        start_reg =         int32(profile_word & (0b1111111111 << 14)) # 10b is 0x3FF
        stop_reg =          int32(profile_word & (0b1111111111 << 30)) # 10b is 0x3FF
        step_interval_mu =  int32(profile_word & 0xFFFFFFFF)
        nodwell_high =      bool(profile_word & (1 << 5))
        zero_crossing =     bool(profile_word & (1 << 3))

        # convert values to human units
        ram_mode_str =      str(RAM_MODE(ram_mode_val))
        step_interval_ns =  step_interval_mu * self.dds.sysclk_per_mu

        # return RAM waveform values
        return ram_mode_str, start_reg, stop_reg, step_interval_ns, nodwell_high, zero_crossing

