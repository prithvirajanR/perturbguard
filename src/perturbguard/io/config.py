from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import yaml


@dataclass
class PerturbGuardConfig:
    schema: dict[str, str] = field(default_factory=dict)
    controls: dict = field(default_factory=dict)
    raw: dict = field(default_factory=dict)


def load_config(path: str | Path | None) -> PerturbGuardConfig:
    if path is None:
        return PerturbGuardConfig()
    path = Path(path)
    data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if not isinstance(data, dict):
        raise ValueError(f"Config {path} must be a YAML top-level mapping.")
    schema = data.get("schema") or {}
    controls = data.get("controls") or {}
    if not isinstance(schema, dict):
        raise ValueError(f"Config {path} field 'schema' must be a mapping.")
    if not isinstance(controls, dict):
        raise ValueError(f"Config {path} field 'controls' must be a mapping.")
    return PerturbGuardConfig(
        schema=dict(schema),
        controls=dict(controls),
        raw=data,
    )
