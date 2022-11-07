import logging

from dax.experiment import *

from dax.modules.led import LedModule
from dax.modules.cpld_init import CpldInitModule

from dax_example.modules.detection import DetectionModule
from dax_example.modules.trap import TrapModule

from dax_example.services.ion_load import IonLoadService
from dax_example.services.state import StateService
from dax_example.services.time_resolved import TimeResolvedService
from dax_example.services.operation import OperationService


class DaxExampleSystem(DaxSystem):
    SYS_ID = 'dax_example_system'
    SYS_VER = 1
    DAX_INFLUX_DB_KEY = None

    def build(self):
        # Adjust logging level
        self.logger.setLevel(min(self.logger.getEffectiveLevel(), logging.INFO))

        # Call super
        super(DaxExampleSystem, self).build()

        # Add standard modules
        self.led = LedModule(self, 'led', *('led0', 'led1'), init_kernel=False)
        self.cpld = CpldInitModule(self, 'cpld', init_kernel=False)
        self.update_kernel_invariants('led', 'cpld')

        # Add modules
        self.detection = DetectionModule(self, 'detection')
        self.trap = TrapModule(self, 'trap')
        self.update_kernel_invariants('detection', 'trap')

        # Add services
        self.ion_load = IonLoadService(self)
        self.state = StateService(self)
        self.time_resolved = TimeResolvedService(self)
        self.operation = OperationService(self)
        self.update_kernel_invariants('ion_load', 'state', 'time_resolved', 'operation')

    @kernel
    def init(self):
        # Call all initialization kernels at once, merging multiple smaller kernels into one
        self.led.init_kernel()
        self.cpld.init_kernel()
        self.detection.init_kernel()
        self.trap.init_kernel()
