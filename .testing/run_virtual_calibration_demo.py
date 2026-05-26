from pathlib import Path
import sys


TEST_ROOT = Path(__file__).resolve().parent
REPO_ROOT = TEST_ROOT.parent          # /Users/yin/vscode/LAX_exp
PACKAGE_PARENT = REPO_ROOT.parent     # /Users/yin/vscode

for path in (PACKAGE_PARENT, TEST_ROOT):
    path_str = str(path)
    if path_str not in sys.path:
        sys.path.insert(0, path_str)

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


def main():
    clock = FakeClock(start=0.0)
    dataset = FakeDatasetTarget()
    backend = FakeSidebandCalibration()
    payload = FakePayload(fail_on_calls={4})

    store = DatasetJsonCalibrationStore(
        dataset_target=dataset,
        json_path=Path(".testing/output/virtual_calibrations.json"),
    )

    policy = CalibrationPolicy(
        every_n_experiments=3,
        max_age_s=60.0,
        calibrate_after_error=True,
        calibrate_if_missing=True,
    )

    manager = CalibrationManager(
        name="sideband",
        policy=policy,
        calibration_backend=backend,
        store=store,
        clock=clock.time,
    )

    for i in range(6):
        print(f"\nPayload attempt {i + 1}")

        try:
            result = manager.run_payload(payload)
            print(f"Payload result: {result}")
        except RuntimeError as exc:
            print(f"Payload failed: {exc}")

        clock.advance(20.0)

        print(f"Calibration calls: {backend.calls}")
        print(
            f"Runs since calibration: {manager.state.runs_since_calibration}")
        print(
            f"Latest red sideband: {dataset.get_dataset('calibration.sideband.red_sideband_freq_hz')}")

    print("\nFinal fake dataset:")
    for key, value in sorted(dataset.data.items()):
        print(f"{key} = {value}")


if __name__ == "__main__":
    main()
