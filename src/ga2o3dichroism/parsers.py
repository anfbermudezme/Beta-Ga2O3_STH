"""Lightweight parsers for CRYSTAL text outputs used in plotting."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from ga2o3dichroism.constants import EV_NM, HARTREE_TO_EV

_FLOAT_PATTERN = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eEdD][-+]?\d+)?"


@dataclass(frozen=True)
class NumericTable:
    """Numeric table read from a whitespace CRYSTAL output file."""

    path: Path
    values: np.ndarray
    columns: tuple[str, ...]

    def to_dataframe(self) -> pd.DataFrame:
        """Return the table as a pandas DataFrame."""
        return pd.DataFrame(self.values, columns=self.columns)


def read_numeric_table(path: str | Path, *, min_columns: int = 2) -> NumericTable:
    """Read numeric rows from a CRYSTAL text file."""
    file_path = Path(path)
    rows = [row for row in _numeric_rows(file_path.read_text(errors="ignore")) if len(row) >= min_columns]
    if not rows:
        raise ValueError(f"No numeric rows with at least {min_columns} columns found in {file_path}.")
    width = max(len(row) for row in rows)
    data = np.full((len(rows), width), np.nan, dtype=float)
    for index, row in enumerate(rows):
        data[index, : len(row)] = row
    columns = tuple(["x", *[f"y{i}" for i in range(1, width)]])
    return NumericTable(path=file_path, values=data, columns=columns)


def read_doss(path: str | Path, *, energy_unit: str = "hartree") -> pd.DataFrame:
    """Read a CRYSTAL DOSS-like file as a DataFrame."""
    table = read_numeric_table(path)
    df = table.to_dataframe().rename(columns={"x": "energy"})
    df["energy_ev"] = df["energy"] * HARTREE_TO_EV if energy_unit.lower().startswith("hart") else df["energy"]
    return df


def read_band(path: str | Path, *, energy_unit: str = "ev") -> pd.DataFrame:
    """Read a simple BAND-like text file as long-form data."""
    table = read_numeric_table(path)
    raw = table.to_dataframe().rename(columns={"x": "k"})
    factor = HARTREE_TO_EV if energy_unit.lower().startswith("hart") else 1.0
    records: list[dict[str, float | int]] = []
    for band_index, column in enumerate(raw.columns[1:], start=1):
        for k_value, energy in zip(raw["k"], raw[column]):
            if not np.isnan(energy):
                records.append({"k": float(k_value), "band": band_index, "energy_ev": float(energy) * factor})
    return pd.DataFrame.from_records(records)


def parse_cphf_output(path: str | Path) -> dict[str, float | str | None]:
    """Extract wavelength, energy, and tensor-like values from a CPHF output."""
    file_path = Path(path)
    text = file_path.read_text(errors="ignore")
    wavelength_nm = _find_wavelength_nm(file_path.name, text)
    energy_ev = EV_NM / wavelength_nm if wavelength_nm else _find_energy_ev(text)
    absorption = _find_last_triplet(text, keywords=("ABSORPTION", "ALPHA"))
    dielectric = _find_last_triplet(text, keywords=("DIELECTRIC", "EPSILON", "PERMITTIVITY"))
    return {
        "file": file_path.name,
        "path": str(file_path),
        "wavelength_nm": wavelength_nm,
        "energy_ev": energy_ev,
        "abs_x": absorption[0] if absorption else None,
        "abs_y": absorption[1] if absorption else None,
        "abs_z": absorption[2] if absorption else None,
        "eps_x": dielectric[0] if dielectric else None,
        "eps_y": dielectric[1] if dielectric else None,
        "eps_z": dielectric[2] if dielectric else None,
    }


def parse_cphf_folder(folder: str | Path) -> pd.DataFrame:
    """Parse all ``*.out`` files in a folder and return one DataFrame."""
    folder_path = Path(folder)
    rows = [parse_cphf_output(path) for path in sorted(folder_path.glob("*.out"))]
    if not rows:
        raise FileNotFoundError(f"No .out files found in {folder_path}.")
    return pd.DataFrame(rows).sort_values(["wavelength_nm", "energy_ev"], na_position="last")


def _numeric_rows(text: str) -> Iterable[list[float]]:
    for line in text.splitlines():
        numbers = []
        for token in line.replace("D", "E").replace("d", "E").replace(",", " ").split():
            try:
                numbers.append(float(token))
            except ValueError:
                continue
        if numbers:
            yield numbers


def _find_wavelength_nm(filename: str, text: str) -> float | None:
    for source in (filename, text[:5000]):
        match = re.search(rf"({_FLOAT_PATTERN})\s*(?:nm|nanometer|nanometre)", source, re.IGNORECASE)
        if match:
            return float(match.group(1).replace("D", "E").replace("d", "E"))
    return None


def _find_energy_ev(text: str) -> float | None:
    match = re.search(rf"({_FLOAT_PATTERN})\s*eV", text, re.IGNORECASE)
    if match:
        return float(match.group(1).replace("D", "E").replace("d", "E"))
    return None


def _find_last_triplet(text: str, *, keywords: tuple[str, ...]) -> tuple[float, float, float] | None:
    result = None
    for line in text.splitlines():
        if any(keyword in line.upper() for keyword in keywords):
            values = [float(value.replace("D", "E").replace("d", "E")) for value in re.findall(_FLOAT_PATTERN, line)]
            if len(values) >= 3:
                result = (values[-3], values[-2], values[-1])
    return result
