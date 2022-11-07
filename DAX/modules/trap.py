import numpy as np

import artiq.coredevice.ad9912
import artiq.coredevice.ttl

from dax.experiment import *


class TrapModule(DaxModule):
    # System dataset keys
    NUM_IONS_KEY = 'num_ions'
    COOL_TIME_KEY = 'cool_time'
    COOL_FREQ_KEY = 'cool_freq'
    COOL_ATT_KEY = 'cool_att'
    COOL_INIT_STATE_KEY = 'cool_init_state'

    def build(self):
        # Cooling DDS and switch
        self._cool_dds = self.get_device('urukul1_ch0', artiq.coredevice.ad9912.AD9912)
        self._cool_sw = self.get_device('ttl_urukul1_sw0', artiq.coredevice.ttl.TTLOut)
        self.update_kernel_invariants('_cool_dds', '_cool_sw')

    @host_only
    def init(self, *, force: bool = False) -> None:
        """Initialize this module.

        :param force: Force full initialization
        """
        # Obtain system datasets
        self._num_ions: np.int32 = self.get_dataset_sys(self.NUM_IONS_KEY, np.int32(2))
        self._cool_time: float = self.get_dataset_sys(self.COOL_TIME_KEY, 1 * us)
        self._cool_freq: float = self.get_dataset_sys(self.COOL_FREQ_KEY, 200 * MHz)
        self._cool_att: float = self.get_dataset_sys(self.COOL_ATT_KEY, 17 * dB)
        self._cool_init_state: bool = self.get_dataset_sys(self.COOL_INIT_STATE_KEY, True)
        self.update_kernel_invariants('_num_ions', '_cool_time', '_cool_freq', '_cool_att', '_cool_init_state')

        if force:
            # Initialize devices
            self.init_kernel()

    @kernel
    def init_kernel(self):
        # For initialization, always reset core first
        self.core.reset()

        # Configure cooling DDS and switch
        self._cool_dds.set_att(self._cool_att)
        self.cool_reset()
        self.cool_set(self._cool_init_state)

        # Wait until all events have been submitted, always required for initialization
        self.core.wait_until_mu(now_mu())

    @host_only
    def post_init(self):
        pass

    """Module functionality"""

    @kernel
    def cool_set(self, state: TBool):
        """Enable or disable cooling (asymmetric operation).

        :param state: State to set
        """
        self._cool_sw.set_o(state)

    @kernel
    def cool_on(self):
        """Enable cooling (asymmetric operation)."""
        self.cool_set(True)

    @kernel
    def cool_off(self):
        """Disable cooling (asymmetric operation)."""
        self.cool_set(False)

    @kernel
    def cool_pulse_mu(self, duration: TInt64 = 0):
        """Cool for a given period of time (symmetric operation).

        :param duration: The pulse duration in machine units
        """
        if duration <= 0:
            # Use default cool time
            duration = self.core.seconds_to_mu(self._cool_time)

        try:
            # Switch devices
            self.cool_on()
            delay_mu(duration)
            self.cool_off()
        except RTIOUnderflow:
            self.core.break_realtime()
            self.cool_off()
            self.core.wait_until_mu(now_mu())
            raise

    @kernel
    def cool_pulse(self, duration: TFloat = 0.0):
        """Cool for a given period of time (symmetric operation).

        :param duration: The pulse duration
        """
        self.cool_pulse_mu(self.core.seconds_to_mu(duration))

    @kernel
    def cool_config(self, freq: TFloat, phase: TFloat = 0.0):
        """Change frequency and phase of microwave DDS.

        :param freq: Frequency
        :param phase: Phase, default = 0.0
        """
        delay(200 * us)
        self._cool_dds.set(freq, phase=phase)

    @kernel
    def cool_config_mu(self, ftw: TInt64, pow_: TInt32 = 0):
        """Change frequency and phase of microwave DDS.

        :param ftw: Frequency tuning word
        :param pow_: Phase offset word, default = 0
        """
        delay(200 * us)
        self._cool_dds.set_mu(ftw, pow=pow_)

    @kernel
    def cool_reset(self):
        """Reset microwave DDS to default frequency and zero phase."""
        delay(200 * us)
        self._cool_dds.set(self._cool_freq)

    """Helper functions"""

    @portable
    def num_ions(self) -> TInt32:
        """Return the number of ions.

        :return: Number of ions
        """
        return self._num_ions

    @host_only
    def set_num_ions(self, num_ions: np.int32):
        """Update the number of ions after loading.

        :param num_ions: The number of ions to set
        """
        # Store value locally and in the dataset
        self._num_ions = np.int32(num_ions)
        self.set_dataset_sys(self.NUM_IONS_KEY, self._num_ions)
