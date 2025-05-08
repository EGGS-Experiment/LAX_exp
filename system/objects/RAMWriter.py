from artiq.experiment import *
from artiq.coredevice.ad9910 import AD9910, RAM_MODE_RAMPUP

from numpy import ceil


class RAMWriter(HasEnvironment):
    """
    Helper: RAM Writer

    Writes waveform data in blocks for AD9910 DDSs to unpredictable RTIOUnderflow errors.
    See here for details: https://github.com/m-labs/artiq/issues/1378
    """
    name = 'RAM Writer'
    kernel_invariants = {
        "dds", "block_size", "dds_profile", "time_block_write_slack_mu"
    }

    def build(self, dds_device: AD9910=None, dds_profile: TInt32=6, block_size: TInt32=88):
        """
        todo: document
        Arguments:
            dds_device: the DDS device object to write for. This MUST be an AD9910 object to prevent silly errors.
            dds_profile: the DDS profile (between [0, 7]) to use while writing RAM. Note that this profile does not
                have to be "reserved," since it is only used during the write_ram procedure. However, any
                preexisting data in this profile will be overwritten.
            block_size: the number of values to write at a time. Must be between [10, 150].
        """
        # get relevant devices
        self.setattr_device('core')

        # todo: check/sanitize input
        # store build arguments
        self.dds = dds_device
        self.block_size = block_size
        self.dds_profile = dds_profile

    def _prepare_argument_checks(self) -> TNone:
        """
        Check experiment arguments for validity.
        """
        if not isinstance(self.dds, AD9910):
            raise ValueError("Invalid device type for RAMWriter. dds_device must be an AD9910 object.")
        if not (10 <= self.block_size <= 150):
            raise ValueError("Invalid block_size for RAMWriter. Block size must be between [10, 150]")
        if not (0 <= self.dds_profile <= 7):
            raise ValueError("Invalid dds_profile for RAMWriter. Profile must be in [0, 7]")

    def prepare(self):
        """
        Prepare hardcoded values ahead of time.
        """
        self._prepare_argument_checks()
        self.time_block_write_slack_mu = self.core.seconds_to_mu(5 * us)

    @kernel(flags={"fast-math"})
    def write(self, ram_data: TList(TInt32), start_addr: TInt32) -> TNone:
        """
        todo: document
        """
        # stop DDS output and disable RAM mode before writing
        self.dds.sw.off()
        self.dds.set_cfr1(ram_enable=0)
        self.dds.cpld.io_update.pulse_mu(8)
        delay_mu(50000) # 50us

        # keep track of current index as we step through RAM data list
        index_current = 0
        num_vals = len(ram_data)

        # write RAM data to AD9910 in small blocks
        # num_write_operations = num_vals // self.block_size + min(num_vals % self.block_size, 0)
        for i in range(round(ceil(num_vals / self.block_size))):

            # update start/stop indices
            addr_current = start_addr + index_current
            update_end = min(self.block_size, num_vals - index_current)
            delay_mu(10000) # 10us

            # # tmp remove
            # delay_mu(1000000)
            # # print(index_current, addr_current, update_end)
            # print(addr_current, addr_current + update_end)
            # print(index_current, index_current + update_end, len(ram_data[index_current: index_current + update_end]))
            # # core_log(index_current, addr_current, update_end)
            # self.core.break_realtime()
            # delay_mu(1000000)
            # # tmp remove

            # step 1: set waveform address range in RAM profile
            # note: step and mode are unimportant b/c user SHOULD overwrite them later
            self.dds.set_profile_ram(
                start=addr_current, end=addr_current + update_end - 1,
                profile=self.dds_profile, step=0xFFF, mode=RAM_MODE_RAMPUP
            )
            self.dds.cpld.io_update.pulse_mu(8)
            delay_mu(512)

            # step 2: set target profile
            self.dds.cpld.set_profile(self.dds_profile)
            self.dds.cpld.io_update.pulse_mu(8)

            # step 3: write RAM
            delay_mu(self.time_block_write_slack_mu)
            self.dds.write_ram(ram_data[index_current: index_current + update_end])

            # update running start address
            index_current += self.block_size

            # allow submitted events to finish execution
            self.core.wait_until_mu(now_mu())
            self.core.break_realtime()

