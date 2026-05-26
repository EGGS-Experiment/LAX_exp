import time
from typing import Any, Callable, Optional

from .policy import CalibrationPolicy, CalibrationState
from .records import CalibrationRecord


class CalibrationManager:
    """
    Pure Python calibration manager.

    This class does not know about DDS, PMT, Urukul, RTIO, or kernels.
    That is intentional.

    It only decides:
        - when to calibrate
        - how to save the result
        - how to respond to payload experiment errors
    """

    def __init__(
        self,
        name: str,
        policy: CalibrationPolicy,
        calibration_backend: Any,
        store: Any,
        clock: Optional[Callable[[], float]] = None,
    ):
        self.name = name
        self.policy = policy
        self.calibration_backend = calibration_backend
        self.store = store
        self.clock = clock or time.time
        self.state = CalibrationState()
        self.latest_record: Optional[CalibrationRecord] = None

    def before_payload(self) -> Optional[CalibrationRecord]:
        now = self.clock()
        decision = self.policy.should_calibrate(self.state, now)

        if not decision.should_calibrate:
            return None

        return self._run_calibration(decision.reason)

    def run_payload(self, payload_callable: Callable[[], Any]) -> Any:
        self.before_payload()

        try:
            result = payload_callable()
        except Exception:
            self.state.runs_since_calibration += 1
            self.state.previous_run_failed = True

            if self.policy.calibrate_after_error:
                self._run_calibration("after_error")

            raise

        self.state.runs_since_calibration += 1
        self.state.previous_run_failed = False
        return result

    def force_calibration(self, reason: str = "manual") -> CalibrationRecord:
        return self._run_calibration(reason)

    def _run_calibration(self, reason: str) -> CalibrationRecord:
        now = self.clock()

        raw_result = self.calibration_backend.run_calibration(self.name)
        record = self._normalize_record(raw_result, now, reason)

        self.store.save(record)

        self.latest_record = record
        self.state.last_calibration_time = record.created_at
        self.state.runs_since_calibration = 0
        self.state.previous_run_failed = False

        return record

    def _normalize_record(
        self,
        raw_result: Any,
        now: float,
        reason: str,
    ) -> CalibrationRecord:
        if isinstance(raw_result, CalibrationRecord):
            return raw_result

        if isinstance(raw_result, dict):
            return CalibrationRecord(
                name=self.name,
                values=raw_result,
                created_at=now,
                valid=True,
                reason=reason,
            )

        raise TypeError(
            "Calibration backend must return either a CalibrationRecord "
            "or a dict of calibration values."
        )
