class FakeClock:
    def __init__(self, start: float = 0.0):
        self._now = float(start)

    def time(self) -> float:
        return self._now

    def advance(self, seconds: float) -> None:
        self._now += float(seconds)

    def set(self, value: float) -> None:
        self._now = float(value)
