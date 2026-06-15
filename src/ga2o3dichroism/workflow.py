"""Small workflow helpers for staging CRYSTAL23 jobs."""

from __future__ import annotations

import shutil
from dataclasses import dataclass
from pathlib import Path

from ga2o3dichroism.crystal import property_input, runpcry_command, runpprop_command

DEFAULT_PROPERTY_KINDS: tuple[str, ...] = ("bands", "doss", "anbd", "echg", "ech3", "ppan")


@dataclass(frozen=True)
class Case:
    """One reproducible CRYSTAL calculation case."""

    name: str
    source_file: str
    label: str


def default_cases() -> tuple[Case, ...]:
    """Return the default β-Ga₂O₃ input cases bundled with the examples."""
    return (
        Case("GaO_TETRA_VAC", "Ga2O3_VAC_TETRA_2x2x2.d12", "2x2x2 neutral tetrahedral Ga vacancy"),
        Case("GaO_OCTA_VAC", "Ga2O3_VAC_OCTA_2x2x2.d12", "2x2x2 neutral octahedral Ga vacancy"),
        Case("GaO_1x2x2_TETRA_VAC", "Ga2O3_1x2x2_TETRA_VAC_SCF.d12", "1x2x2 tetrahedral Ga vacancy"),
        Case("GaO_1x2x2_OCTA_VAC", "Ga2O3_1x2x2_OCTA_VAC_SCF.d12", "1x2x2 octahedral Ga vacancy"),
        Case("GaO_1x2x2_TETRA_VACS", "Ga2O3_1x2x2_TETRA_VACS_SCF.d12", "1x2x2 paired tetrahedral Ga vacancies"),
        Case("GaO_1x2x2_OCTA_OCTA_VACS", "Ga2O3_1x2x2_OCTA_OCTA_VACS_SCF.d12", "1x2x2 paired octahedral Ga vacancies"),
        Case("GaO_1x2x2_OCTA_TETRA_VACS", "Ga2O3_1x2x2_OCTA_TETRA_VACS_SCF.d12", "1x2x2 mixed octahedral/tetrahedral Ga vacancies"),
        Case("GaO_OVAC_1", "Ga2O3_Ovac_1_SCF.d12", "oxygen vacancy site 1"),
        Case("GaO_OVAC_2", "Ga2O3_Ovac_2_SCF_SPIN.d12", "oxygen vacancy site 2"),
        Case("GaO_OVAC_3", "Ga2O3_Ovac_3_SCF.d12", "oxygen vacancy site 3"),
    )


def stage_crystal_workflow(
    *,
    input_dir: str | Path,
    out_dir: str | Path,
    cores: int = 128,
    cases: tuple[Case, ...] | None = None,
    property_kinds: tuple[str, ...] = DEFAULT_PROPERTY_KINDS,
) -> Path:
    """Create a clean run tree for CRYSTAL23 jobs."""
    input_path = Path(input_dir)
    out_path = Path(out_dir)
    scf_path = out_path / "scf"
    prop_path = out_path / "properties"
    scripts_path = out_path / "scripts"
    scf_path.mkdir(parents=True, exist_ok=True)
    prop_path.mkdir(parents=True, exist_ok=True)
    scripts_path.mkdir(parents=True, exist_ok=True)

    selected_cases = cases or default_cases()
    scf_commands = ["#!/usr/bin/env bash", "set -euo pipefail", ""]
    prop_commands = ["#!/usr/bin/env bash", "set -euo pipefail", ""]

    for case in selected_cases:
        source = input_path / case.source_file
        if not source.exists():
            raise FileNotFoundError(f"Missing source input for {case.name}: {source}")
        shutil.copy2(source, scf_path / f"{case.name}.d12")
        scf_commands.append(runpcry_command(cores, case.name))
        for kind in property_kinds:
            prop_name = f"{kind.upper()}_{case.name}"
            (prop_path / f"{prop_name}.d3").write_text(property_input(kind, case.name), encoding="utf-8")
            prop_commands.append(runpprop_command(cores, prop_name, case.name))

    _write_script(scripts_path / "run_scf.sh", scf_commands)
    _write_script(scripts_path / "run_properties.sh", prop_commands)
    return out_path


def _write_script(path: Path, lines: list[str]) -> None:
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    path.chmod(0o755)
