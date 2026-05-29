class FakeSidebandCalibration:
    """
    Fake sideband calibration.

    Each call returns slightly different values so tests can check
    whether a new calibration actually happened.
    """

    def __init__(self):
        self.calls = 0

    def run_calibration(self, name: str) -> dict:
        self.calls += 1

        return {
            "red_sideband_freq_hz": 80.0e6 + self.calls * 1.0e3,
            "blue_sideband_freq_hz": 82.0e6 + self.calls * 1.0e3,
            "pi_time_s": 25.0e-6,
            "call_index": self.calls,
        }
