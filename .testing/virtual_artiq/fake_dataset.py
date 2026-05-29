class FakeDatasetTarget:
    """
    Small replacement for ARTIQ's set_dataset/get_dataset behavior.
    """

    def __init__(self):
        self.data = {}

    def set_dataset(
        self,
        key,
        value,
        persist=False,
        archive=True,
    ):
        self.data[key] = value

    def get_dataset(self, key, default=None):
        return self.data.get(key, default)

    def __getitem__(self, key):
        return self.data[key]

    def __contains__(self, key):
        return key in self.data
