import json
from pathlib import Path
from typing import Any, MutableMapping, Optional

from .records import CalibrationRecord


class DatasetJsonCalibrationStore:
    """
    Saves calibration results in two places.

    1. ARTIQ-style datasets
    2. A JSON file

    The dataset target can be either:
        - an ARTIQ EnvExperiment object with set_dataset()
        - a normal Python dict for virtual testing
    - dataset_target: Any object that supports set_dataset(key, value) or dict-like item assignment
    - json_path: Path to the JSON file where calibration records will be saved
    - dataset_prefix: Prefix for dataset keys (default: "calibration")
    - persist=True
        Store the dataset in the ARTIQ master's on-disk dataset database.
        This makes the value survive across experiments and master restarts.
        ARTIQ says persist also implies broadcast.
        Use `persist=True` when you want later experiments to reuse the value:

    - archive=True
        Save this dataset into the current experiment run's local result file,
        usually the HDF5 result file for that run.
        Use `archive=True` when you want the calibration value copied into the result file of the current run.
    """

    def __init__(
        self,
        dataset_target: Any,
        json_path: str | Path,
        dataset_prefix: str = "calibration",
        persist: bool = True,
        archive: bool = True,
    ):
        self.dataset_target = dataset_target
        self.json_path = Path(json_path)
        self.dataset_prefix = dataset_prefix
        self.persist = persist
        self.archive = archive

    def save(self, record: CalibrationRecord) -> None:
        self._save_to_datasets(record)
        self._save_to_json(record)

    def load_latest(self, name: str) -> Optional[CalibrationRecord]:
        if not self.json_path.exists():
            return None

        with self.json_path.open("r", encoding="utf-8") as f:
            data = json.load(f)

        latest = data.get("latest", {}).get(name)
        if latest is None:
            return None

        return CalibrationRecord.from_dict(latest)

    def _save_to_datasets(self, record: CalibrationRecord) -> None:
        base_key = f"{self.dataset_prefix}.{record.name}"

        for key, value in record.values.items():
            self._set_dataset(f"{base_key}.{key}", value)

        self._set_dataset(f"{base_key}.created_at", record.created_at)
        self._set_dataset(f"{base_key}.valid", record.valid)
        self._set_dataset(f"{base_key}.reason", record.reason)

    def _set_dataset(self, key: str, value: Any) -> None:
        if hasattr(self.dataset_target, "set_dataset"):
            try:
                self.dataset_target.set_dataset(
                    key,
                    value,
                    persist=self.persist,
                    archive=self.archive,
                )
            except TypeError:
                self.dataset_target.set_dataset(key, value)
            return

        if isinstance(self.dataset_target, MutableMapping):
            self.dataset_target[key] = value
            return

        raise TypeError(
            "dataset_target must be an ARTIQ-like object with set_dataset() "
            "or a mutable mapping like dict."
        )

    def _save_to_json(self, record: CalibrationRecord) -> None:
        self.json_path.parent.mkdir(parents=True, exist_ok=True)

        if self.json_path.exists():
            with self.json_path.open("r", encoding="utf-8") as f:
                data = json.load(f)
        else:
            data = {
                "latest": {},
                "history": [],
            }

        record_dict = record.to_dict()

        data["latest"][record.name] = record_dict
        data["history"].append(record_dict)

        with self.json_path.open("w", encoding="utf-8") as f:
            json.dump(data, f, indent=2, sort_keys=True, default=str)
