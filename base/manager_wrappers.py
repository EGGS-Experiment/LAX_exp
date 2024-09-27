"""
Base classes for building PyQt5 GUI clients for LabRAD.
"""
import logging
from importlib import import_module
from artiq.master.worker_db import _write

from LAX_exp.system.device_db_ext import device_db_ext

logger = logging.getLogger("artiq.master.experiments")


def _write_to_group(group, k, v):
    """
    Wrap ARTIQs base "_write" function with slightly more error handling
    and add support for dictionaries.
    """
    try:
        if type(v) is not dict:
            _write(group, k, v)
        else:
            # convert any dicts to text
            _write(group, k, str(v))
    except Exception as e:
        logger.warning("Error: unable to write key {} to group {}: {}".format(k, group, e))


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
                # logger.warning("Parameter '%s' is already in parameter storage, "
                #                "overwriting", key, stack_info=True)
                pass
            self.parameters[key] = data

        return data

    def write_hdf5(self, f):
        """
        Add handling for stored parameters.
        """
        # store datasets in a separate group
        datasets_group = f.create_group("datasets")
        for k, v in self.local.items():
            _write_to_group(datasets_group, k, v)

        # store archived datasets (e.g. calibrations) in the archive group
        archive_group = f.create_group("archive")
        for k, v in self.archive.items():
            _write_to_group(archive_group, k, v)

        # store parameters in a separate group as attributes
        parameters_group = f.create_group("parameters")
        for k, v in self.parameters.items():
            _write_to_group(parameters_group.attrs, k, v)

        # store arguments in a separate group as attributes
        arguments_group = f.create_group("arguments")
        for k, v in self.arguments.items():
            _write_to_group(arguments_group.attrs, k, v)


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
                print('\tError: {}'.format(repr(e)))

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
