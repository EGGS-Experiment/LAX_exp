from artiq.experiment import *
from artiq.coredevice import ad9910
from numpy import int32, int64

from enum import Enum
# todo: also read urukul status register & config register


class RAM_MODE(Enum):
    """
    Enum class to convert the RAM mode bits to a string.
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
    kernel_invariants = {
        '_stored_profile_words', '_stored_ram_data_lists',
    }

    def build(self):
        self.setattr_device("core")

        # DDS parameters
        self.setattr_argument("dds_target",     EnumerationValue(list(self._get_dds_devices()), default='urukul0_ch0'),
                              tooltip="DDS channel to read from.")

        # profile parameters
        self.setattr_argument('dds_profiles',   PYONValue([0, 1, 2, 3, 4, 5, 6, 7]), group="config",
                              tooltip="DDS profiles to configure with the selected parameters.")

        # RAM profile configuration
        self.setattr_argument("is_profile_ram", BooleanValue(default=False),
                              tooltip="Whether the DDS profile should be interpreted as a RAM profile.")
        self.setattr_argument("read_ram_data",  BooleanValue(default=False),
                              tooltip="Whether the RAM data in the addresses specified by the RAM profile should be read.")

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
        ### SANITIZE INPUTS ###
        # ensure no duplicates in dds_profiles
        self.dds_profiles = list(set(self.dds_profiles))

        # ensure profile number is valid, otherwise strange failure mode/dangerous overwrite
        if not isinstance(self.dds_profiles, list):
            raise ValueError("Invalid argument: dds_profiles must be a list.")
        # ensure dds_profiles are int in [0, 7]
        dds_profiles_valid = [isinstance(val, int) and (0 <= val <= 7) for val in self.dds_profiles]
        if not all(dds_profiles_valid):
            raise ValueError("Invalid DDS profile list. Must be an int in [0, 7].")

        # ensure we only read RAM data if profile is a RAM profile
        if (self.read_ram_data is True) & (self.is_profile_ram is False):
            raise Exception("Error: cannot grab RAM data from non-RAM profile.")


        # get relevant devices
        try:
            self.dds = self.get_device(self.dds_target)
            self.dds_cpld = self.dds.cpld

            # add device names to kernel invariants
            kernel_invariants = getattr(self, "kernel_invariants", set())
            self.kernel_invariants = kernel_invariants | {'dds', 'dds_cpld'}
        except Exception as e:
            print("Error: invalid DDS channel.")
            raise e


        # set up data structure to get data for tone and ram profile
        self._stored_att_db = float(0.) # store RAM data indicated by the profile
        self._stored_profile_words = [] # store retrieved profile data
        self._stored_ram_data_lists = [] # store RAM data indicated by the profile


    """
    MAIN SEQUENCE
    """
    @kernel(flags={"fast-math"})
    def run_initialize(self) -> TNone:
        """
        Initialize devices prior to execution.
        """
        self.core.reset()

        # ensure all DDS outputs are OFF to prevent booboos
        self.dds.cpld.cfg_switches(0b0000)
        self.dds.sw.off()
        # todo: set switches of actual DDSs (i.e. dds.sw) as well
        self.core.break_realtime()

        # store state of att register to prevent booboos
        self.dds.cpld.get_att_mu()
        self.core.break_realtime()

        # # initialize urukul (idk why but everyone does it each time)
        # self.dds.cpld.init()
        # self.core.break_realtime()
        # delay_mu(100000)
        #
        # # initialize DDS (idk why but everyone does it each time)
        # # warning - this may ruin timing since it tunes sync_delay
        # # self.dds.init()
        # self.core.break_realtime()
        # delay_mu(100000)

    @kernel(flags={"fast-math"})
    def run(self) -> TNone:
        # prepare devices for experiment
        self.run_initialize()

        # get global/CPLD values
        self.core.break_realtime()
        self._stored_att_db = self.dds.get_att()
        # todo: read urukul status register (sw/smp/pll/ifc/proto_rev)
        # todo: read urukul config register (sw/led/profile/io_update/mask_nu/clk_sel/sync_sel/rst/clk_div)


        # prepare DDS in desired profile (necessary since values are only taken from active profile)
        for idx_profile in range(len(self.dds_profiles)):
            # get dds profile
            profile_num = self.dds_profiles[idx_profile]

            # set target DDS profile
            self.core.break_realtime()
            self.dds_cpld.set_profile(profile_num)
            self.dds_cpld.io_update.pulse_mu(8)

            # read directly from register (instead of using get_mu) since profile may be RAM
            self.core.break_realtime()
            profile_word_tmp = int64(self.dds.read64(ad9910._AD9910_REG_PROFILE0 + profile_num))
            self.store_profile_word(profile_word_tmp)

            # read data from RAM
            if self.read_ram_data:
                # parse profile word for RAM begin/end addresses and create commensurate list for storage
                _tmp0, addr_start, addr_stop, _tmp1, _tmp2, _tmp3 = self.process_profile_word_ram(profile_word_tmp)
                _ram_data_tmp = [int32(0)] * (addr_stop - addr_start + 1)

                # ensure RAM profile is valid before retrieving RAM data
                if (addr_start < addr_stop) and (0 <= addr_start <= 1023) and (0 <= addr_stop <= 1023):
                    # retrieve RAM data (thankfully, coredevice.ad9910 handles it for us)
                    self.core.break_realtime()
                    delay_mu(1000000)  # 10ms - heavy slack (full transfer is 1024 32b words @ 125/16 MHz => ~4.2ms)
                    self.dds.read_ram(_ram_data_tmp)
                    self.store_ram_data(_ram_data_tmp)  # store results in host via RPC

                # otherwise, store dummy value and move on
                else:
                    self.store_ram_data([-1])  # store results in host via RPC
                    continue

        # clean up by reverting to default profile
        self.core.break_realtime()
        self.dds_cpld.set_profile(ad9910.DEFAULT_PROFILE)
        self.dds_cpld.io_update.pulse_mu(8)
        self.core.wait_until_mu(now_mu())
        self.core.break_realtime()


    """
    DATA PROCESSING
    """
    def analyze(self):
        """
        Process results to extract DDS data.
        """
        print("DDS - {}".format(self.dds_target))

        # process global information
        self.set_dataset("dds_att_db", self._stored_att_db)
        print("\tAtt.:\t{} dB".format(self._stored_att_db))


        # process all profile data
        for idx_profile, profile_num in enumerate(self.dds_profiles):
            print("\tProfile {:d}:".format(profile_num))

            # process profile data as single-tone waveform
            if not self.is_profile_ram:
                # extract waveform values from 64b single-tone profile word
                freq_hz, ampl_frac, phas_turns = self.process_profile_word_singletone(
                    self._stored_profile_words[idx_profile])

                # store and print results
                self.set_dataset("profile_{:d}_freq_hz".format(profile_num), freq_hz)
                self.set_dataset("profile_{:d}_ampl_frac".format(profile_num), ampl_frac)
                self.set_dataset("profile_{:d}_phas_turns".format(profile_num), phas_turns)
                print("\t\tFreq.: {:.4f} MHz".format(freq_hz / MHz))
                print("\t\tAmpl.: {:.4f} %".format(ampl_frac * 100.))
                print("\t\tPhas.: {:.4f} turns".format(phas_turns))

            # process profile data as RAM
            else:
                # extract values from RAM profile word
                (ram_mode_str, start_addr, stop_addr, step_interval_ns,
                 nodwell_high, zero_crossing) = self.process_profile_word_ram(self._stored_profile_words[idx_profile])

                # store and print results
                self.set_dataset("profile_{:d}_ram_mode".format(profile_num), ram_mode_str)
                self.set_dataset("profile_{:d}_start_addr".format(profile_num), start_addr)
                self.set_dataset("profile_{:d}_stop_addr".format(profile_num), stop_addr)
                self.set_dataset("profile_{:d}_step_interval_ns".format(profile_num), step_interval_ns)
                self.set_dataset("profile_{:d}_nodwell_high".format(profile_num), nodwell_high)
                self.set_dataset("profile_{:d}_zero_crossing".format(profile_num), zero_crossing)
                print("\t\tRAM Mode:\t{:s}".format(ram_mode_str))
                print("\t\tStart reg.:\t{:d}".format(start_addr))
                print("\t\tStop reg.:\t{:d}".format(stop_addr))
                print("\t\tStep rate:\t{:.4f} ns".format(step_interval_ns))
                print("\t\tNo-dwell high:\t{:b}".format(nodwell_high))
                print("\t\tZero-crossing:\t{:b}".format(zero_crossing))

                # store retrieved RAM data (don't print out since it may be fat)
                if self.read_ram_data:
                    self.set_dataset("profile_{:d}_ram_data_raw".format(profile_num),
                                     self._stored_ram_data_lists[idx_profile])


    """
    HELPER FUNCTIONS
    """
    @portable
    def process_profile_word_singletone(self, profile_word: TInt64) -> TTuple([TFloat, TFloat, TFloat]):
        """
        Extract single-tone waveform from the DDS profile word.
        :param profile_word: 64b single-tone profile word from the AD9910.
        :return: a tuple of (freq_hz, ampl_frac, phas_turns)
        """
        # extract values from profile word
        asf = int32((profile_word >> 48) & 0x3FFF) # 14b (0x3FFF)
        pow = int32((profile_word >> 32) & 0xFFFF) # 16b (0xFFFF)
        ftw = int32(profile_word & ((1 << 32) - 1)) # 32b (0xFFFFFFFF)

        # convert values to human units and return
        freq_hz =       self.dds.ftw_to_frequency(ftw)
        ampl_frac =     self.dds.asf_to_amplitude(asf)
        phas_turns =    self.dds.pow_to_turns(pow)
        return freq_hz, ampl_frac, phas_turns

    @rpc
    def process_profile_word_ram(self, profile_word: TInt64) -> TTuple([TStr, TInt32, TInt32, TFloat, TBool, TBool]):
        """
        Extract RAM waveform from DDS profile word.
        :param profile_word: 64b RAM profile word from the AD9910.
        :return: a tuple of (ram_mode_str, start_addr, stop_addr, step_interval_ns, nodwell_high_flag, zero_crossing_flag)
        """
        # extract values from profile word
        step_interval_mu =  int32((profile_word >> 40) & 0xFFFF) # 16b
        stop_addr =         int32((profile_word >> 30) & 0x3FF) # 10b (0x3FF)
        start_addr =        int32((profile_word >> 14) & 0x3FF) # 10b (0x3FF)
        nodwell_high =      bool((profile_word >> 5) & 1)
        zero_crossing =     bool((profile_word >> 3) & 1)
        ram_mode_val =      int32(profile_word & 0x7) # 3b is 0x7

        # convert values to human units before returning
        ram_mode_str =      str(RAM_MODE(ram_mode_val))
        step_interval_ns =  float(step_interval_mu * self.dds.sysclk_per_mu)
        return ram_mode_str, start_addr, stop_addr, step_interval_ns, nodwell_high, zero_crossing

    @rpc(flags={"async"})
    def store_profile_word(self, profile_word: TInt64) -> TNone:
        """
        Store retrieved profile word on host-side.
        :param profile_word: the profile word to store.
        """
        self._stored_profile_words.append(profile_word)

    @rpc(flags={"async"})
    def store_ram_data(self, data_arr: TList(TInt32)) -> TNone:
        """
        Store retrieved RAM data on host-side.
        :param data_arr: the RAM data to store.
        """
        # simply append RAM data list to our holder
        self._stored_ram_data_lists.append(data_arr)

