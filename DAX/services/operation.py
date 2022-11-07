import numpy as np

from dax.experiment import *
from dax.interfaces.operation import OperationInterface

from dax_example.modules.detection import DetectionModule
from dax_example.modules.trap import TrapModule
from dax_example.services.state import StateService


# noinspection PyAbstractClass
class OperationService(DaxService, OperationInterface):
    """A dummy operation service."""

    SERVICE_NAME = 'operation'

    def build(self):
        # Add kernel invariants for properties
        self.update_kernel_invariants('pi', 'num_qubits', '_channel_map')

        # Required modules and services
        self._detect = self.registry.find_module(DetectionModule)
        self._trap = self.registry.find_module(TrapModule)
        self._state = self.registry.get_service(StateService)
        self.update_kernel_invariants('_detect', '_trap', '_state')

    def init(self) -> None:
        pass

    def post_init(self) -> None:
        pass

    """Service functionality"""

    @rpc(flags={'async'})
    def _log_gate(self, name, *qubits, **kwargs):
        if any(q < 0 or q >= len(self._channel_map) for q in qubits):
            raise IndexError('Qubit index out of range')
        args = [f'{k}={v}' for k, v in kwargs.items()]
        args.extend((str(q) for q in qubits))
        self.logger.info(f'{name}({",".join(args)})')

    @kernel
    def i(self, qubit: TInt32):
        self._log_gate('I', qubit)

    @kernel
    def x(self, qubit: TInt32):
        self._log_gate('X', qubit)

    @kernel
    def y(self, qubit: TInt32):
        self._log_gate('Y', qubit)

    @kernel
    def z(self, qubit: TInt32):
        self._log_gate('Z', qubit)

    @kernel
    def h(self, qubit: TInt32):
        self._log_gate('H', qubit)

    @kernel
    def sqrt_x(self, qubit: TInt32):
        self._log_gate('sqrt_X', qubit)

    @kernel
    def sqrt_x_dag(self, qubit: TInt32):
        self._log_gate('sqrt_X_dag', qubit)

    @kernel
    def sqrt_y(self, qubit: TInt32):
        self._log_gate('sqrt_Y', qubit)

    @kernel
    def sqrt_y_dag(self, qubit: TInt32):
        self._log_gate('sqrt_Y_dag', qubit)

    @kernel
    def sqrt_z(self, qubit: TInt32):
        self._log_gate('sqrt_Z', qubit)

    @kernel
    def sqrt_z_dag(self, qubit: TInt32):
        self._log_gate('sqrt_Z_dag', qubit)

    @kernel
    def rx(self, theta: TFloat, qubit: TInt32):
        self._log_gate('rX', qubit, theta=theta)

    @kernel
    def ry(self, theta: TFloat, qubit: TInt32):
        self._log_gate('rY', qubit, theta=theta)

    @kernel
    def rz(self, theta: TFloat, qubit: TInt32):
        self._log_gate('rZ', qubit, theta=theta)

    @kernel
    def rphi(self, theta: TFloat, phi: TFloat, qubit: TInt32):
        self._log_gate('rphi', qubit, theta=theta, phi=phi)

    @kernel
    def xx(self, control: TInt32, target: TInt32):
        self._log_gate('xx', control, target)

    @kernel
    def xx_dag(self, control: TInt32, target: TInt32):
        self._log_gate('xx_dag', control, target)

    @kernel
    def rxx(self, theta: TFloat, control: TInt32, target: TInt32):
        self._log_gate('rxx', control, target, theta=theta)

    @kernel
    def cz(self, control: TInt32, target: TInt32):
        self._log_gate('cnot', control, target)

    @kernel
    def cnot(self, control: TInt32, target: TInt32):
        self._log_gate('cnot', control, target)

    @kernel
    def prep_0_all(self):
        self._trap.cool_pulse()

    @kernel
    def m_z_all(self):
        self._detect.detect_active()

    @kernel
    def get_measurement(self, qubit: TInt32) -> TInt32:
        return self._detect.measure(self._channel_map[qubit])

    @kernel
    def store_measurements(self, qubits: TList(TInt32)):
        self._state.measure_channels([self._channel_map[q] for q in qubits])

    @property
    def num_qubits(self) -> np.int32:
        return np.int32(self._trap.num_ions())

    @property
    def _channel_map(self):
        """A map to convert qubit index to channel."""
        return self._detect.active_channels()

    @host_only
    def set_realtime(self, realtime: bool) -> None:
        pass
