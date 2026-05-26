from dataclasses import dataclass, field
from typing import Any, Mapping


@dataclass(frozen=True)
class CalibrationRecord:
    name: str
    values: Mapping[str, Any]
    created_at: float
    valid: bool = True
    reason: str = "manual"
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "values": dict(self.values),
            "created_at": self.created_at,
            "valid": self.valid,
            "reason": self.reason,
            "metadata": dict(self.metadata),
        }

    @classmethod
    def from_dict(cls, data: dict) -> "CalibrationRecord":
        return cls(
            name=data["name"],
            values=data.get("values", {}),
            created_at=float(data["created_at"]),
            valid=bool(data.get("valid", True)),
            reason=data.get("reason", "manual"),
            metadata=data.get("metadata", {}),
        )
