import numpy as np
from artiq.experiment import *
from artiq.coredevice import ad9910


class RAMWriter(HasEnvironment):
    """
    Helper: RAM Writer

    Writes waveform data in blocks for AD9910 DDSs.
    """
    name = 'RAM Writer'
    kernel_invariants = {
        "dds", "block_size", "dds_profile", "time_block_write_slack_mu"
    }

    def build(self, dds_device=None, dds_profile=6, block_size=50):
        """
        todo: document
        """
        # get relevant devices
        self.setattr_device('core')

        # todo: check/sanitize input
        # store build arguments
        self.dds = dds_device
        self.block_size = block_size
        self.dds_profile = dds_profile

        # hardcoded values
        self.time_block_write_slack_mu = self.core.seconds_to_mu(250 * us)

    @kernel(flags={"fast-math"})
    def write(self, ram_data: TList(TInt32), start_addr: TInt32) -> TNone:
        """
        todo: document
        """
        # stop DDS output and disable RAM mode before writing
        self.dds.off()
        self.dds.set_cfr1(ram_enable=0)
        delay_mu(50000) # 50us

        # keep track of current index as we step through RAM data list
        index_current = 0
        num_vals = len(ram_data)

        # write RAM data to AD9910 in small blocks
        num_write_operations = num_vals // self.block_size + min(num_vals % self.block_size, 0)
        for i in range(num_write_operations):

            # update start/stop indices
            addr_current = start_addr + index_current
            update_end = min(self.block_size, num_vals - index_current) - 1
            delay_mu(10000) # 10us

            # # tmp remove
            # delay_mu(1000000)
            # print(index_current, addr_current, update_end)
            # # core_log(index_current, addr_current, update_end)
            # self.core.break_realtime()
            # delay_mu(1000000)
            # # tmp remove

            # set waveform address range in RAM profile
            # note: step and mode are unimportant b/c user SHOULD overwrite them later
            self.dds.set_profile_ram(
                start=addr_current, end=addr_current + update_end,
                profile=self.dds_profile, step=0x001, mode=ad9910.RAM_MODE_RAMPUP
            )
            self.dds.cpld.io_update.pulse_mu(8)
            delay_mu(512)
            self.dds.cpld.set_profile(self.dds_profile)
            self.dds.cpld.io_update.pulse_mu(8)

            # add sizable slack, then write RAM
            delay_mu(self.time_block_write_slack_mu) # 500us
            self.dds.write_ram(ram_data[index_current: index_current + update_end])

            # update running start address
            index_current += self.block_size

            # allow submitted events to finish execution
            self.core.wait_until_mu(now_mu())
            self.core.break_realtime()

