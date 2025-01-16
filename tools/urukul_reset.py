from artiq.experiment import *
from artiq.coredevice.urukul import CFG_RST


class UrukulReset(EnvExperiment):
    """
    Tool: Urukul Reset

    Reset an Urukul board via the MASTER_RESET
    """
    name = 'Urukul Reset'
    kernel_invariants = {
        'urukul_cpld', 'urukul_channels'
    }

    def build(self):
        self.setattr_device("core")

        # DDS parameters
        dds_device_list = self._get_urukul_devices()
        self.setattr_argument("urukul_target", EnumerationValue(list(dds_device_list), default='urukul0_cpld'))

    def _get_urukul_devices(self) -> TList(TStr):
        """
        Get all valid Urukul devices from the device_db.
        """
        def is_local_cpld_device(v):
            return isinstance(v, dict) and (v.get('type') == 'local') and (v.get('class') == 'CPLD')

        # get only local urukul CPLD devices from device_db
        return list([k for k, v in self.get_device_db().items() if is_local_cpld_device(v)])

    def _get_urukul_channels(self, cpld_name: TStr) -> TList(TStr):
        """
        Get all child DDS channels of a given urukul CPLD.
        Arguments:
            cpld_name: name of the parent urukul CPLD.
        """
        def is_local_dds_device(v):
            return (isinstance(v, dict) and (v.get('type') == 'local') and (v.get('class') == "AD9910") and
                    (v.get('arguments', {}).get('cpld_device') == cpld_name))

        # get only local urukul CPLD devices from device_db
        return list([k for k, v in self.get_device_db().items() if is_local_dds_device(v)])

    def prepare(self):
        """
        Prepare values for speedy evaluation.
        """
        try:
            # get urukul CPLD
            self.setattr_device(self.urukul_target)
            self.urukul_cpld = self.get_device(self.urukul_target)

            # get all child DDSs of the CPLD
            self.urukul_channels = list()
            for channel_name in self._get_urukul_channels(self.urukul_target):
                self.setattr_device(channel_name)
                self.urukul_channels.append(self.get_device(channel_name))

        except Exception as e:
            print("Error: invalid Urukul board.")
            raise e

    @kernel(flags={"fast-math"})
    def run(self) -> TNone:
        # reset core
        self.core.reset()

        # pulse MASTER_RESET signal on urukul's CPLD
        self.urukul_cpld.cfg_write(self.urukul_cpld.cfg_reg | (1 << CFG_RST))
        delay_mu(100000) # 100 us
        self.urukul_cpld.cfg_write(self.urukul_cpld.cfg_reg & ~(1 << CFG_RST))

        # re-initialize urukul CPLD
        delay_mu(5000000) # 5 ms
        self.urukul_cpld.init()

        # initialize channels on the urukul
        for channel in self.urukul_channels:
            channel.init()
            self.core.break_realtime()
            delay_mu(5000000) # 5 ms

        # synchronize timeline
        self.core.wait_until_mu(now_mu())
        self.core.break_realtime()

