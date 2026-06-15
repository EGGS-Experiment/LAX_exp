from pathlib import Path
import sys

from artiq.experiment import EnvExperiment, NumberValue, BooleanValue
THIS_FILE = Path(__file__).resolve()
REPO_ROOT = THIS_FILE.parents[2]
PACKAGE_PARENT = REPO_ROOT.parent
from LAX_exp.calibration import (
    CalibrationManager,
    CalibrationPolicy,
    DatasetJsonCalibrationStore,
)

class FakeSidebandCalibration:
    def __init__(self):
        self.calls = 0

    def run_calibration(self, name):
        self.calls += 1

        return {
            "red_sideband_freq_hz": 80.0e6 + self.calls * 1.0e3,
            "blue_sideband_freq_hz": 82.0e6 + self.calls * 1.0e3,
            "pi_time_s": 25.0e-6,
            "call_index": self.calls,
        }
        
class HostOnlyCalibrationTest(EnvExperiment):
    def build(self):
        self.num_payloads = self.get_argument(
            "num_payloads",
            NumberValue(default=6, ndecimals=0, step=1),
        )

        self.fail_on_fourth_payload = self.get_argument(
            "fail_on_fourth_payload",
            BooleanValue(default=True),
        )
    def prepare(self):
        self.calibration_backend = FakeSidebandCalibration()

        store = DatasetJsonCalibrationStore(
            dataset_target=self,
            json_path=REPO_ROOT / ".testing" / "output" / "artiq_host_calibrations.json",
            dataset_prefix="calibration",
            persist=True,
            archive=True,
        )

        policy = CalibrationPolicy(
            every_n_experiments=3,
            max_age_s=60.0,
            calibrate_after_error=True,
            calibrate_if_missing=True,
        )

        self.manager = CalibrationManager(
            name="sideband",
            policy=policy,
            calibration_backend=self.calibration_backend,
            store=store,
        )

    def run(self):
        for i in range(int(self.num_payloads)):
            print("")
            print("Payload attempt", i + 1)

            try:
                self.manager.run_payload(lambda: self.fake_payload(i + 1))
            except RuntimeError as exc:
                print("Payload failed:", exc)

            red_freq = self.get_dataset("calibration.sideband.red_sideband_freq_hz")
            call_index = self.get_dataset("calibration.sideband.call_index")
            reason = self.get_dataset("calibration.sideband.reason")

            print("Calibration calls:", self.calibration_backend.calls)
            print("Runs since calibration:", self.manager.state.runs_since_calibration)
            print("Latest red sideband:", red_freq)
            print("Latest call index:", call_index)
            print("Latest reason:", reason)

    def fake_payload(self, payload_index):
        if self.fail_on_fourth_payload and payload_index == 4:
            raise RuntimeError("controlled fake host-only payload failure")

        return {
            "payload_index": payload_index,
            "success": True,
        }