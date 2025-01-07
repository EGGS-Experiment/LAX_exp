from artiq.experiment import *
from artiq.coredevice.urukul import CFG_RST


class UrukulReset(EnvExperiment):
    """
    Tool: Urukul Reset

    Reset an Urukul board via the MASTER_RESET
    """
    name = 'Urukul Reset'

    def build(self):
        self.setattr_device("core")

        # DDS parameters
        dds_device_list = self._get_urukul_devices()
        self.setattr_argument("urukul_target", EnumerationValue(list(dds_device_list), default='urukul0_cpld'))

    def _get_urukul_devices(self):
        """
        Get all valid Urukul devices from the device_db.
        """
        def is_local_dds_device(v):
            return isinstance(v, dict) and (v.get('type') == 'local') and ('class' in v) and (v.get('class') == "CPLD")

        # get only local urukul CPLD devices from device_db
        return set([k for k, v in self.get_device_db().items() if is_local_dds_device(v)])

    def prepare(self):
        """
        Prepare values for speedy evaluation.
        """
        '''GET DEVICES'''
        try:
            self.urukul_cpld = self.get_device(self.urukul_target)

            # add DDS name to kernel invariants
            kernel_invariants = getattr(self, "kernel_invariants", set())
            self.kernel_invariants = kernel_invariants | {'urukul_cpld'}

        except Exception as e:
            print("Error: invalid Urukul board.")
            raise e


    """
    MAIN FUNCTIONS
    """
    @kernel(flags={"fast-math"})
    def run(self) -> TNone:
        # reset core
        self.core.reset()
        self.core.break_realtime()

        # pulse MASTER_RESET signal on urukul's CPLD
        self.urukul_cpld.cfg_write(self.urukul_cpld.cfg_reg | (1 << CFG_RST))
        delay_mu(100000) # 100 us
        self.urukul_cpld.cfg_write(self.urukul_cpld.cfg_reg & ~(1 << CFG_RST))

        # re-initialize urukul CPLD
        delay_mu(5000000) # 5 ms
        self.urukul_cpld.init()

        # synchronize timeline
        self.core.wait_until_mu(now_mu())
        self.core.break_realtime()
        delay_mu(1000000) # 1 ms

