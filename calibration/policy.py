from dataclasses import dataclass
from typing import Optional


@dataclass
class CalibrationState:
    last_calibration_time: Optional[float] = None
    runs_since_calibration: int = 0
    previous_run_failed: bool = False


@dataclass(frozen=True)
class CalibrationDecision:
    should_calibrate: bool
    reason: str = "not_needed"


@dataclass(frozen=True)
class CalibrationPolicy:
    every_n_experiments: Optional[int] = None
    max_age_s: Optional[float] = None
    calibrate_after_error: bool = True
    calibrate_if_missing: bool = True

    def should_calibrate(
        self,
        state: CalibrationState,
        now: float,
    ) -> CalibrationDecision:
        if self.calibrate_if_missing and state.last_calibration_time is None:
            return CalibrationDecision(True, "missing")

        if self.calibrate_after_error and state.previous_run_failed:
            return CalibrationDecision(True, "after_error")

        if (
            self.every_n_experiments is not None
            and self.every_n_experiments > 0
            and state.runs_since_calibration >= self.every_n_experiments
        ):
            return CalibrationDecision(True, "run_count")

        if (
            self.max_age_s is not None
            and state.last_calibration_time is not None
            and now - state.last_calibration_time >= self.max_age_s
        ):
            return CalibrationDecision(True, "time_interval")

        return CalibrationDecision(False, "not_needed")
