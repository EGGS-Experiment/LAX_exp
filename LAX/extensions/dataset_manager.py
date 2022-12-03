# todo: device groups w/attributes
# todo: sequence groups w/attributes
# todo: redo expid
# todo: update parameters
# todo: specify save location
# todo: specify
from artiq.master.worker_db import DatasetManager

class LAXDatasetManager(DatasetManager):
    """
    todo: document
    """

    def set(self, key, value, broadcast=False, persist=False, archive=True, parameter=False):
        if persist:
            broadcast = True

        if broadcast:
            self._broadcaster[key] = persist, value
        elif key in self._broadcaster.raw_view:
            del self._broadcaster[key]

        if archive:
            self.local[key] = value
        elif key in self.local:
            del self.local[key]

        # todo: notate
        if parameter:
            self.parameters[key] = value

    def write_hdf5(self, f):
        datasets_group = f.create_group("datasets")
        for k, v in self.local.items():
            _write(datasets_group, k, v)

        archive_group = f.create_group("archive")
        for k, v in self.archive.items():
            _write(archive_group, k, v)

        # todo: notate
        parameters_group = datasets_group.attrs
        for k, v in self.parameters.items():
            _write(parameters_group, k, v)

    def _write(self, group, k, v):
        # Add context to exception message when the user writes a dataset that is
        # not representable in HDF5.
        try:
            group[k] = v
        except TypeError as e:
            raise TypeError("Error writing dataset '{}' of type '{}': {}".format(
                k, type(v), e))
