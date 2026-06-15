"""High-level CRYSTAL23 input generation workflows.

The functions in this module create a complete, case-centered run tree
from compact Python objects or from a YAML/JSON workflow file.  The
intention is to make the computational protocol reproducible without
requiring users to manually rename CRYSTAL input decks or rebuild the
same property files by hand.
"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

import yaml

from ga2o3dichroism.crystal import CphfSettings, insert_cphf_block, property_input, runpcry_command, runpprop_command
from ga2o3dichroism.workflow import Case, DEFAULT_PROPERTY_KINDS, default_cases


@dataclass(frozen=True)
class CphfJob:
    """Definition of one CRYSTAL23 CPHF/CPKS input deck.

    Parameters
    ----------
    name:
        Name of the generated CRYSTAL job, without the ``.d12`` suffix.
    settings:
        Static or dynamic response settings written into the ``CPHF`` block.
    label:
        Human-readable description stored in the manifest.
    """

    name: str
    settings: CphfSettings = field(default_factory=CphfSettings.static)
    label: str = ""

    @classmethod
    def static(cls, case_name: str, *, name: str | None = None, label: str = "") -> "CphfJob":
        """Create a static CPHF/CPKS job for ``case_name``."""
        return cls(name=name or f"CPHF_{case_name}_STATIC", settings=CphfSettings.static(), label=label or "static CPHF/CPKS response")

    @classmethod
    def dynamic(
        cls,
        case_name: str,
        *,
        start_nm: float,
        stop_nm: float | None = None,
        steps: int = 1,
        damping: float | None = None,
        name: str | None = None,
        label: str = "",
    ) -> "CphfJob":
        """Create a frequency-dependent CPHF/CPKS job for ``case_name``."""
        suffix = _dynamic_suffix(start_nm=start_nm, stop_nm=stop_nm)
        return cls(
            name=name or f"CPHF_{case_name}_{suffix}",
            settings=CphfSettings.dynamic(start_nm=start_nm, stop_nm=stop_nm, steps=steps, damping=damping),
            label=label or "dynamic CPHF/CPKS response",
        )


@dataclass(frozen=True)
class InputCase:
    """Complete input-generation definition for one CRYSTAL23 case.

    Parameters
    ----------
    name:
        Canonical CRYSTAL job name used for generated files and run commands.
    source_file:
        Existing SCF/geometry ``.d12`` template, relative to ``input_dir``.
    label:
        Human-readable case description.
    properties:
        Property inputs to generate as ``.d3`` files.  Accepted values are the
        same names handled by :func:`ga2o3dichroism.crystal.property_input`.
    cphf_jobs:
        Full CPHF/CPKS ``.d12`` inputs generated from ``source_file``.
    """

    name: str
    source_file: str
    label: str = ""
    properties: tuple[str, ...] = DEFAULT_PROPERTY_KINDS
    cphf_jobs: tuple[CphfJob, ...] = ()

    @classmethod
    def from_case(
        cls,
        case: Case,
        *,
        properties: Iterable[str] = DEFAULT_PROPERTY_KINDS,
        cphf_jobs: Iterable[CphfJob] = (),
    ) -> "InputCase":
        """Create an input-generation case from a lightweight workflow case."""
        return cls(
            name=case.name,
            source_file=case.source_file,
            label=case.label,
            properties=tuple(properties),
            cphf_jobs=tuple(cphf_jobs),
        )


@dataclass(frozen=True)
class GeneratedFile:
    """Metadata for one generated input file."""

    kind: str
    case_name: str
    job_name: str
    path: Path
    command: str = ""
    label: str = ""


@dataclass(frozen=True)
class GenerationResult:
    """Summary returned after generating a CRYSTAL23 run tree."""

    root: Path
    files: tuple[GeneratedFile, ...]
    manifest: Path
    scripts: tuple[Path, ...]


def publication_cases() -> tuple[InputCase, ...]:
    """Return the default manuscript-oriented Ga-vacancy workflow.

    The two main neutral 2×2×2 Ga-vacancy cases receive static and dynamic
    CPHF/CPKS inputs in addition to the standard post-SCF property inputs.
    Remaining bundled example cases receive only the post-SCF property files.
    """
    selected: list[InputCase] = []
    for case in default_cases():
        cphf_jobs: tuple[CphfJob, ...] = ()
        if case.name in {"GaO_TETRA_VAC", "GaO_OCTA_VAC"}:
            cphf_jobs = (
                CphfJob.static(case.name),
                CphfJob.dynamic(case.name, start_nm=7500, stop_nm=12500, steps=10, damping=0.002),
            )
        selected.append(InputCase.from_case(case, cphf_jobs=cphf_jobs))
    return tuple(selected)


def load_workflow(path: str | Path) -> tuple[InputCase, ...]:
    """Load a case-generation workflow from a YAML or JSON file.

    The accepted schema is intentionally small: a top-level ``defaults`` mapping
    and a ``cases`` list.  Each case can override ``properties`` and define a
    ``cphf`` list with ``mode: static`` or ``mode: dynamic``.
    """
    workflow_path = Path(path)
    data = _read_mapping(workflow_path)
    defaults = data.get("defaults", {}) or {}
    default_properties = tuple(defaults.get("properties", DEFAULT_PROPERTY_KINDS))
    raw_cases = data.get("cases")
    if not isinstance(raw_cases, list) or not raw_cases:
        raise ValueError(f"Workflow file must define a non-empty 'cases' list: {workflow_path}")

    parsed_cases: list[InputCase] = []
    for raw_case in raw_cases:
        if not isinstance(raw_case, dict):
            raise ValueError("Each workflow case must be a mapping.")
        name = str(raw_case["name"])
        raw_source = raw_case.get("source") or raw_case.get("source_file")
        if raw_source is None or not str(raw_source).strip():
            raise ValueError(f"Case {name!r} must define 'source'.")
        source = str(raw_source)
        properties = tuple(raw_case.get("properties", default_properties))
        cphf_jobs = tuple(_parse_cphf_job(case_name=name, raw_job=raw_job) for raw_job in raw_case.get("cphf", []) or [])
        parsed_cases.append(
            InputCase(
                name=name,
                source_file=source,
                label=str(raw_case.get("label", "")),
                properties=properties,
                cphf_jobs=cphf_jobs,
            )
        )
    return tuple(parsed_cases)


def generate_crystal_inputs(
    *,
    input_dir: str | Path,
    out_dir: str | Path,
    cores: int = 128,
    cases: Iterable[InputCase] | None = None,
) -> GenerationResult:
    """Generate all CRYSTAL23 input files and run scripts for ``cases``.

    The generated directory is organized by case so that each folder contains
    the SCF ``.d12``, property ``.d3`` files, optional CPHF/CPKS ``.d12`` files,
    and a local ``run_case.sh`` script.  Global scripts are also written under
    ``scripts/`` for batch submission or manual inspection.
    """
    input_path = Path(input_dir)
    root = Path(out_dir)
    case_root = root / "crystal23"
    scripts_root = root / "scripts"
    case_root.mkdir(parents=True, exist_ok=True)
    scripts_root.mkdir(parents=True, exist_ok=True)

    selected_cases = tuple(cases or publication_cases())
    generated: list[GeneratedFile] = []
    scf_commands = _global_script_header()
    cphf_commands = _global_script_header()
    property_commands = _global_script_header()
    all_commands = _global_script_header()

    for case in selected_cases:
        source_path = input_path / case.source_file
        if not source_path.exists():
            raise FileNotFoundError(f"Missing source input for {case.name}: {source_path}")
        case_dir = case_root / case.name
        case_dir.mkdir(parents=True, exist_ok=True)
        case_commands = ["#!/usr/bin/env bash", "set -euo pipefail", ""]

        scf_path = case_dir / f"{case.name}.d12"
        scf_text = source_path.read_text(encoding="utf-8")
        scf_path.write_text(scf_text, encoding="utf-8")
        scf_command = runpcry_command(cores, case.name)
        generated.append(GeneratedFile("scf", case.name, case.name, scf_path, scf_command, case.label))
        _append_cd_command(scf_commands, case.name, scf_command)
        _append_cd_command(all_commands, case.name, scf_command)
        case_commands.append(scf_command)

        for cphf_job in case.cphf_jobs:
            cphf_path = case_dir / f"{cphf_job.name}.d12"
            cphf_path.write_text(insert_cphf_block(scf_text, cphf_job.settings), encoding="utf-8")
            cphf_command = runpcry_command(cores, cphf_job.name)
            generated.append(GeneratedFile("cphf", case.name, cphf_job.name, cphf_path, cphf_command, cphf_job.label))
            _append_cd_command(cphf_commands, case.name, cphf_command)
            _append_cd_command(all_commands, case.name, cphf_command)
            case_commands.append(cphf_command)

        for property_kind in case.properties:
            normalized_kind = _normalize_property_kind(property_kind)
            property_job = f"{normalized_kind.upper()}_{case.name}"
            property_path = case_dir / f"{property_job}.d3"
            property_path.write_text(property_input(normalized_kind, case.name), encoding="utf-8")
            property_command = runpprop_command(cores, property_job, case.name)
            generated.append(GeneratedFile("property", case.name, property_job, property_path, property_command, normalized_kind))
            _append_cd_command(property_commands, case.name, property_command)
            _append_cd_command(all_commands, case.name, property_command)
            case_commands.append(property_command)

        _write_script(case_dir / "run_case.sh", case_commands)

    manifest = root / "manifest.csv"
    _write_manifest(manifest, root, generated)
    scripts = (
        scripts_root / "run_scf.sh",
        scripts_root / "run_cphf.sh",
        scripts_root / "run_properties.sh",
        scripts_root / "run_all.sh",
    )
    _write_script(scripts[0], scf_commands)
    _write_script(scripts[1], cphf_commands)
    _write_script(scripts[2], property_commands)
    _write_script(scripts[3], all_commands)
    return GenerationResult(root=root, files=tuple(generated), manifest=manifest, scripts=scripts)


def _read_mapping(path: Path) -> dict[str, Any]:
    if path.suffix.lower() in {".json"}:
        data = json.loads(path.read_text(encoding="utf-8"))
    else:
        data = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"Workflow file must contain a mapping: {path}")
    return data


def _parse_cphf_job(*, case_name: str, raw_job: Any) -> CphfJob:
    if not isinstance(raw_job, dict):
        raise ValueError(f"CPHF job for {case_name!r} must be a mapping.")
    mode = str(raw_job.get("mode", "static")).lower()
    name = raw_job.get("name")
    label = str(raw_job.get("label", ""))
    if mode == "static":
        return CphfJob.static(case_name, name=str(name) if name else None, label=label)
    if mode == "dynamic":
        if "start_nm" not in raw_job:
            raise ValueError(f"Dynamic CPHF job for {case_name!r} must define 'start_nm'.")
        return CphfJob.dynamic(
            case_name,
            start_nm=float(raw_job["start_nm"]),
            stop_nm=float(raw_job["stop_nm"]) if raw_job.get("stop_nm") is not None else None,
            steps=int(raw_job.get("steps", 1)),
            damping=float(raw_job["damping"]) if raw_job.get("damping") is not None else None,
            name=str(name) if name else None,
            label=label,
        )
    raise ValueError(f"Unknown CPHF mode for {case_name!r}: {mode!r}")


def _normalize_property_kind(kind: str) -> str:
    key = kind.lower().strip()
    aliases = {"band": "bands", "dos": "doss"}
    return aliases.get(key, key)


def _dynamic_suffix(*, start_nm: float, stop_nm: float | None) -> str:
    def fmt(value: float) -> str:
        text = f"{value:g}".replace(".", "p")
        return f"{text}NM"

    if stop_nm is None:
        return fmt(start_nm)
    return f"{fmt(start_nm)}_{fmt(stop_nm)}"


def _global_script_header() -> list[str]:
    return [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        r'SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"',
        r'ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"',
        "",
    ]


def _append_cd_command(lines: list[str], case_name: str, command: str) -> None:
    lines.extend([f'cd "${{ROOT_DIR}}/crystal23/{case_name}"', command, ""])


def _write_script(path: Path, lines: list[str]) -> None:
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")
    path.chmod(0o755)


def _write_manifest(path: Path, root: Path, files: Iterable[GeneratedFile]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=("kind", "case", "job", "path", "command", "label"))
        writer.writeheader()
        for item in files:
            writer.writerow(
                {
                    "kind": item.kind,
                    "case": item.case_name,
                    "job": item.job_name,
                    "path": item.path.relative_to(root).as_posix(),
                    "command": item.command,
                    "label": item.label,
                }
            )
