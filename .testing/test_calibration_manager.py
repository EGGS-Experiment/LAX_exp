import json

import pytest

from LAX_exp.calibration import (
    CalibrationManager,
    CalibrationPolicy,
    DatasetJsonCalibrationStore,
)
from virtual_artiq import (
    FakeClock,
    FakeDatasetTarget,
    FakePayload,
    FakeSidebandCalibration,
)


def make_manager(tmp_path, policy):
    clock = FakeClock(start=0.0)
    dataset = FakeDatasetTarget()
    backend = FakeSidebandCalibration()

    store = DatasetJsonCalibrationStore(
        dataset_target=dataset,
        json_path=tmp_path / "calibrations.json",
    )

    manager = CalibrationManager(
        name="sideband",
        policy=policy,
        calibration_backend=backend,
        store=store,
        clock=clock.time,
    )

    return manager, clock, dataset, backend, store


def test_first_payload_calibrates_if_missing(tmp_path):
    policy = CalibrationPolicy()
    manager, clock, dataset, backend, store = make_manager(tmp_path, policy)

    payload = FakePayload()
    result = manager.run_payload(payload)

    assert result["success"] is True
    assert backend.calls == 1
    assert dataset["calibration.sideband.call_index"] == 1
    assert dataset["calibration.sideband.valid"] is True


def test_calibration_is_saved_to_json(tmp_path):
    policy = CalibrationPolicy()
    manager, clock, dataset, backend, store = make_manager(tmp_path, policy)

    manager.before_payload()

    json_path = tmp_path / "calibrations.json"

    with json_path.open("r", encoding="utf-8") as f:
        data = json.load(f)

    latest = data["latest"]["sideband"]

    assert latest["name"] == "sideband"
    assert latest["values"]["call_index"] == 1
    assert len(data["history"]) == 1


def test_calibrates_every_n_experiments(tmp_path):
    policy = CalibrationPolicy(every_n_experiments=2)
    manager, clock, dataset, backend, store = make_manager(tmp_path, policy)

    payload = FakePayload()

    manager.run_payload(payload)
    assert backend.calls == 1

    manager.run_payload(payload)
    assert backend.calls == 1

    manager.run_payload(payload)
    assert backend.calls == 2
    assert dataset["calibration.sideband.call_index"] == 2


def test_calibrates_after_time_interval(tmp_path):
    policy = CalibrationPolicy(max_age_s=10.0)
    manager, clock, dataset, backend, store = make_manager(tmp_path, policy)

    manager.before_payload()
    assert backend.calls == 1

    clock.advance(5.0)
    manager.before_payload()
    assert backend.calls == 1

    clock.advance(5.0)
    manager.before_payload()
    assert backend.calls == 2


def test_calibrates_immediately_after_payload_error(tmp_path):
    policy = CalibrationPolicy(calibrate_after_error=True)
    manager, clock, dataset, backend, store = make_manager(tmp_path, policy)

    payload = FakePayload(fail_on_calls={1})

    with pytest.raises(RuntimeError):
        manager.run_payload(payload)

    assert backend.calls == 2
    assert dataset["calibration.sideband.call_index"] == 2
    assert manager.state.runs_since_calibration == 0
    assert manager.state.previous_run_failed is False
