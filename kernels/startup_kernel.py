"""
Copied over from DAX via the startup kernel page in their wiki: https://gitlab.com/duke-artiq/dax/-/snippets/1966946/raw/master/startup_kernel.py
Added minor formatting changes.
"""
from artiq.experiment import *


class StartupKernel(EnvExperiment):
    """
    Generic auto-generated startup kernel class.

    Devices can be excluded from initialization by adding their device keys to :attr:`EXCLUDE`.

    To add a new device class in need of initialization at startup:
    1. Add the device class name to :attr:`DEVICE_CLASSES`.
    2. Add a kernel function for device initialization and add appropriate delays, see existing functions below.

    Kernel functions for device initialization of device type ``DeviceType`` should be named ``devicetype``
    and the list of devices of that type can be accessed through the attribute ``self.devicetype_devices``.
    """
    
    """
    USER CONFIGURATION
    """
    
    # Key of the core device (only useful for DRTIO systems)
    CORE_DEVICE = 'core'

    # Set of devices to exclude (aliases will be resolved)
    EXCLUDE = set()

    # List of devices classes to be initialized
    DEVICE_CLASSES = [
        'CPLD',         # note: Urukul CPLDs have to be initialized before the DDS channels themselves
        'AD9910',
        'Phaser',
        'Fastino',
        'Sampler'
    ]
    
    # The amount of time reserved for an initialization to finish.
    DELAY_INIT_MU = 5000000 # 5 ms


    """
    INITIALIZATION FUNCTIONS
    """
    @kernel
    def cpld(self):
        for d in self.cpld_devices:
            delay_mu(self.DELAY_INIT_MU)
            d.init()

    @kernel
    def ad53xx(self):
        for d in self.ad53xx_devices:
            delay_mu(self.DELAY_INIT_MU)
            d.init()

    @kernel
    def ad9910(self):
        for d in self.ad9910_devices:
            delay_mu(self.DELAY_INIT_MU)
            d.init()

            # set waveform
            delay_mu(self.DELAY_INIT_MU)
            d.set_mu(0x1C28F5C2, asf=0x1FFF)    # 110MHz, 50% amplitude

    @kernel
    def ad9912(self):
        for d in self.ad9912_devices:
            delay_mu(self.DELAY_INIT_MU)
            d.init()

    @kernel
    def ad9914(self):
        for d in self.ad9914_devices:
            delay_mu(self.DELAY_INIT_MU)
            d.init()

    @kernel
    def adf5356(self):
        for d in self.adf5356_devices:
            delay_mu(self.DELAY_INIT_MU)
            d.init()

    @kernel
    def fastino(self):
        for d in self.fastino_devices:
            delay_mu(self.DELAY_INIT_MU)
            d.init()

    @kernel
    def mirny(self):
        for d in self.mirny_devices:
            delay_mu(self.DELAY_INIT_MU)
            d.init()

    @kernel
    def phaser(self):
        for d in self.phaser_devices:
            delay_mu(self.DELAY_INIT_MU)
            d.init()

    @kernel
    def sampler(self):
        for d in self.sampler_devices:
            delay_mu(self.DELAY_INIT_MU)
            d.init()

    @kernel
    def suservo(self):
        for d in self.suservo_devices:
            delay_mu(self.DELAY_INIT_MU)
            d.init()

    @kernel
    def zotino(self):
        for d in self.zotino_devices:
            delay_mu(self.DELAY_INIT_MU)
            d.init()


    """
    NO USER CHANGES REQUIRED BELOW THIS COMMENT
    """
    
    kernel_invariants = {'core', '_init_functions'}

    def build(self):
        assert isinstance(self.EXCLUDE, set), 'Exclude keys must be a set'
        assert all(isinstance(k, str) for k in self.EXCLUDE), 'All excluded keys must be of type str'
        assert isinstance(self.CORE_DEVICE, str), 'Core device key must be of type str'
        assert isinstance(self.DEVICE_CLASSES, list), 'Device classes must be a list'
        assert all(isinstance(c, str) for c in self.DEVICE_CLASSES), 'All device classes must be of type str'
        assert len(set(self.DEVICE_CLASSES)) == len(self.DEVICE_CLASSES), 'Device classes list contains duplicate'
        assert all(hasattr(self, c.lower()) for c in self.DEVICE_CLASSES), 'Not all device classes have an init fn'

        try:
            from dax.sim.ddb import DAX_SIM_CONFIG_KEY
        except ImportError:
            pass  # DAX not available
        else:
            if DAX_SIM_CONFIG_KEY in self.get_device_db():
                # DAX.sim is enabled
                print('WARNING: DAX.sim is enabled')

        # Message
        print('Generating startup kernel...')

        # Core device
        self.core = self.get_device(self.CORE_DEVICE)
        # Resolve all excluded keys
        exclude = self.resolve_keys(self.EXCLUDE)
        # Get all local devices
        all_devices = self.get_local_devices()
        # List of initialization functions
        self._init_functions = []

        # For every device class with init functions, see if we have such devices and set attributes accordingly
        for device_class in self.DEVICE_CLASSES:
            # Get the device list for this class
            all_device_keys = [k for k, v in all_devices.items() if v['class'] == device_class]
            excluded_device_keys = [k for k in all_device_keys if k in exclude]
            devices = {k: self.get_device(k) for k in all_device_keys if k not in exclude}

            if devices:
                # Message user a device type was found
                print('Initializing {} device(s) of type {}: {}'.format(len(devices), device_class, ', '.join(devices)))

                # Expected name of the initialization function for this device
                init_func_name = device_class.lower()
                # Add initialization function to list
                self._init_functions.append(getattr(self, init_func_name))

                # Add supplementary attribute to self
                init_attr_name = '{}_devices'.format(init_func_name)
                setattr(self, init_attr_name, list(devices.values()))
                self.kernel_invariants.add(init_attr_name)

            if excluded_device_keys:
                # Message user devices were excluded
                print('Excluded {} device(s) of type {}: {}'.format(
                    len(excluded_device_keys), device_class, ', '.join(excluded_device_keys)))

        if not self._init_functions:
            # Handle case where there are no init commands
            self._init_commands = self._no_init_commands
            # Message
            print('No devices were initialized')

        # Report about devices that were not initialized
        not_initialized = {k: v for k, v in all_devices.items() if v['class'] not in self.DEVICE_CLASSES}
        print('Device types not initialized: {}'.format(', '.join({v['class'] for v in not_initialized.values()})))
        print('Devices not initialized: {}'.format(', '.join(not_initialized.keys())))

    def resolve_keys(self, keys):
        """
        Resolve all aliases of a set of keys and return the union.
        """

        resolved_keys = set(keys)
        ddb = self.get_device_db()

        def resolve(key):
            item = ddb[key]
            if isinstance(item, str):
                resolved_keys.add(item)
                resolve(item)

        try:
            for k in keys:
                resolve(k)
        except KeyError as e:
            raise KeyError(f'Key "{e.args[0]}" could not be resolved in the device DB') from None

        return resolved_keys

    def get_local_devices(self):
        """
        Get local devices from the device DB.
        """

        def is_local_device(v):
            return isinstance(v, dict) and v.get('type') == 'local' and 'class' in v

        # Filter local devices from device DB
        return {k: v for k, v in self.get_device_db().items() if is_local_device(v)}

    @kernel
    def run(self):
        """
        Kernel function that runs all the required initialization functions.
        """

        # RTIO reset
        self.core.reset()

        # Process init commands
        self._init_commands()

        # RTIO sync
        self.core.wait_until_mu(now_mu())

    @kernel
    def _init_commands(self):
        """
        Call all initialization functions.
        """
        for f in self._init_functions:
            f()
            self.core.break_realtime()  # Guarantee positive slack

    @kernel
    def _no_init_commands(self):
        """
        Dummy function used in case there are no devices to initialize.
        """
        pass