class FakePayload:
    """
    Fake payload experiment.

    fail_on_calls={2} means the second payload run fails.
    """

    def __init__(self, fail_on_calls=None):
        self.calls = 0
        self.fail_on_calls = set(fail_on_calls or [])

    def __call__(self):
        self.calls += 1

        if self.calls in self.fail_on_calls:
            raise RuntimeError(f"Fake payload failed on call {self.calls}")

        return {
            "payload_call": self.calls,
            "success": True,
        }
