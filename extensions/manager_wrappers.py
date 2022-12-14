"""
Base classes for building PyQt5 GUI clients for LabRAD.
"""
import logging
from importlib import import_module
from artiq.master.worker_db import _write

logger = logging.getLogger("artiq.master.experiments")
from LAX_exp.extensions.device_db_ext import device_db_ext


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

        # add parameter and argument management
        self.__dict__['parameters'] =           dict()
        self.__dict__['arguments'] =            dict()

    def __getattr__(self, attr):
        return getattr(self.base, attr)

    def __setattr__(self, attr, value):
        setattr(self.base, attr, value)


    # MODIFIED
    def set(self, key, value, broadcast=False, persist=False, archive=True, parameter=False, argument=False):
        """
        Added handling for arguments.
        """
        if argument:
            self.arguments[key] = value
        elif parameter:
            self.parameters[key] = value
        else:
            self.base.set(key, value, broadcast, persist, archive)

    def get(self, key, archive=False, parameter=False):
        """
        Add handling for parameters.
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

        # handle parameters
        if parameter:
            if key in self.parameters:
                logger.warning("Parameter '%s' is already in parameter storage, "
                               "overwriting", key, stack_info=True)
            self.parameters[key] = data

        return data

    def write_hdf5(self, f):
        """
        Add handling for stored parameters.
        """
        self.base.write_hdf5(f)
        # todo: remove handling for archive groups

        # store parameters in a separate group as attributes
        parameters_group = f.create_group("parameters")
        for k, v in self.parameters.items():
            _write(parameters_group.attrs, k, v)

        # store arguments in a separate group as attributes
        arguments_group = f.create_group("arguments")
        for k, v in self.arguments.items():
            _write(arguments_group.attrs, k, v)


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
        self.__dict__['ddb_lax'] =              device_db_ext
        self.__dict__['active_lax_devices'] =   list()

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
            for desc, obj in self.active_lax_devices:
                if dev_desc == desc:
                    return obj

            # instantiate object and add to holding dictionary
            try:
                dev_obj = self.__create_lax_device(dev_desc)
                self.active_lax_devices.append((dev_desc, dev_obj))
                return dev_obj
            except Exception as e:
                print('Unable to create LAX Device: {}'.format(name))
                print('\tError: {}'.format(e))

        # get device normally
        else:
            return self.base.get(name)

    def __create_lax_device(self, desc):
        """
        Create and instantiate a given LAXDevice.
        """
        if desc["type"] == "local":
            module = import_module(desc["module"])
            device_class = getattr(module, desc["class"])
            return device_class(self.parent, **desc.get("arguments", {}))