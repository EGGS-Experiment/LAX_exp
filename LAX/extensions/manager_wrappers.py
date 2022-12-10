# todo: create parameter group that stores parameters separately for each class object using attributes
# todo: specify save location
# todo: create experiment group that saves arguments, save_location, etc. there
# todo: allow writing by enevp to diff location

from importlib import import_module
from artiq.master.worker_db import DatasetManager, DeviceManager, _write

logger = logging.getLogger("artiq.master.experiments")
from LAX_exp.LAX.extensions.device_db_ext import device_db


class LAXDatasetManager:
    """
    Wraps around the DatasetManager object used by the master to instantiate HasEnvironment objects.
    """

    # WRAPPERS
    def __init__(self, base, parent):
        # set base (wrapped) object
        self.__dict__['base'] =                 base

        # get parent (i.e. object which owns the HasEnvironment object, typically an EnvExperiment)
        self.__dict__['parent'] =               parent

        # add parameter management
        self.__dict__['parameters'] =           dict()

    def __getattr__(self, attr):
        return getattr(self.base, attr)

    def __setattr__(self, attr, value):
        setattr(self.base, attr, value)


    # MODIFIED
    def get(self, key, archive=False, parameter=False):
        """
        Add handling for parameter arguments.
        """
        # regular handling
        if key in self.local:
            return self.local[key]

        data = self.ddb.get(key)
        if archive:
            if key in self.archive:
                logger.warning("Dataset '%s' is already in archive, "
                               "overwriting", key, stack_info=True)
            self.archive[key] = data

        # handle parameters as well
        if parameter:
            if key in self.parameter:
                logger.warning("Dataset '%s' is already in archive, "
                               "overwriting", key, stack_info=True)
            self.parameter[key] = archive

        return data

    def write_hdf5(self, f):
        """
        Add handling for stored parameters.
        """
        # handle datasets normally
        datasets_group = f.create_group("datasets")
        for k, v in self.local.items():
            _write(datasets_group, k, v)

        # don't handle archive group

        # store parameters in a separate group
        parameters_group = f.create_group("parameters")
        for k, v in self.parameters.items():
            _write(parameters_group.attrs, k, v)


class LAXDeviceManager:
    """
    Wraps around the DeviceManager object used by the master to instantiate HasEnvironment objects.
    """

    # WRAPPERS
    def __init__(self, base, parent):
        # set base (wrapped) object
        self.__dict__['base'] =                 base

        # get parent (i.e. object which owns the HasEnvironment object, typically an EnvExperiment)
        self.__dict__['parent'] =               parent

        # add device holders for LAX
        self.__dict__['ddb_lax'] =              device_db
        self.__dict__['active_lax_devices'] =   dict()

    def __getattr__(self, attr):
        return getattr(self.base, attr)

    def __setattr__(self, attr, value):
        setattr(self.base, attr, value)


    # MODIFIED
    def get(self, name):
        """
        Redirect requests for LAXDevices
        """
        # redirect requests for LAXDevices
        if name in self.ddb_lax:

            # get device description
            dev_desc = self.ddb_lax[name]

            # return object if it already exists
            if dev_desc in self.active_lax_devices:
                return self.active_lax_devices[dev_desc]

            # instantiate object and add to holding dictionary
            else:
                dev_obj = self.__create_lax_device(dev_desc)
                self.active_lax_devices[name] = dev_obj

        # get device normally
        else:
            return self.base.get(name)

    def __create_lax_device(self, desc):
        if desc["type"] == "local":
            module = import_module(desc["module"])
            device_class = getattr(module, desc["class"])
            return device_class(self.parent, **desc.get("arguments", {}))
