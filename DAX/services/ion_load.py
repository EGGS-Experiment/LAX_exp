from dax.experiment import *

from dax_example.modules.detection import DetectionModule
from dax_example.modules.trap import TrapModule


class IonLoadService(DaxService):
    SERVICE_NAME = 'ion_load'

    def build(self) -> None:
        # Get modules
        self._detection = self.registry.find_module(DetectionModule)
        self._trap = self.registry.find_module(TrapModule)

    def init(self) -> None:
        pass

    def post_init(self) -> None:
        pass

    """Service functionality"""

    @host_only
    def load_ions(self, num_ions: int) -> None:
        """Load ions into the trap.

        :param num_ions: Number of ions to load
        """
        assert isinstance(num_ions, int) and 0 < num_ions <= 4

        # There is no example code for loading ions, just update the modules accordingly
        self._trap.set_num_ions(num_ions)
        self._detection.set_active_channels(list(range(num_ions)))
        self.logger.info(f'Loaded {num_ions} ion(s)')
