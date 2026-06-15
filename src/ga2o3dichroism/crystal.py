"""CRYSTAL23 input builders used in the Ga₂O₃ defect-dichroism workflow."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

PropertyKind = Literal["bands", "doss", "anbd", "echg", "ech3", "ppan"]


@dataclass(frozen=True)
class BandPath:
    """Integer reciprocal-space path for a CRYSTAL ``BAND`` input.

    Parameters
    ----------
    denominator:
        Integer denominator used by CRYSTAL for the k-point coordinates.
    subdivisions:
        Number of subdivisions along the full path.
    first_band, last_band:
        Band-index window written to the property input.
    points:
        Ordered list of labels and integer coordinates. Consecutive points are
        converted into CRYSTAL line segments.
    """

    denominator: int = 120
    subdivisions: int = 240
    first_band: int = 1
    last_band: int = 180
    iplo: int = 1
    lpr66: int = 0
    points: tuple[tuple[str, tuple[int, int, int]], ...] = (
        ("X1", (-38, 82, 0)),
        ("Y", (0, 60, 0)),
        ("G", (0, 0, 0)),
        ("N", (-60, 60, 0)),
        ("X", (-60, 60, 60)),
        ("G", (0, 0, 0)),
        ("M", (0, 60, 60)),
        ("I", (-38, 82, 60)),
        ("L", (0, 60, 60)),
        ("F", (-60, 60, 60)),
        ("Y", (0, 60, 0)),
        ("G", (0, 0, 0)),
        ("Z", (0, 0, 60)),
        ("F1", (-60, 60, 60)),
        ("Z", (0, 0, 60)),
        ("I1", (-38, 82, 60)),
    )

    @property
    def n_segments(self) -> int:
        """Number of line segments in the path."""
        return len(self.points) - 1


@dataclass(frozen=True)
class CphfSettings:
    """Settings for a CRYSTAL ``CPHF`` / ``CPKS`` block.

    ``wavelengths_nm`` controls dynamic calculations. If it is ``None``, a
    static CPHF block is generated. If it contains one value, a single dynamic
    wavelength is written. If it contains two values, CRYSTAL receives the start
    and stop wavelengths with ``steps`` points.
    """

    wavelengths_nm: tuple[float, ...] | None = None
    steps: int = 1
    damping: float | None = None
    maxcycle: int = 200
    tolalpha: int = 6

    @classmethod
    def static(cls, *, maxcycle: int = 200, tolalpha: int = 6) -> "CphfSettings":
        """Create settings for a static dielectric-response calculation."""
        return cls(wavelengths_nm=None, maxcycle=maxcycle, tolalpha=tolalpha)

    @classmethod
    def dynamic(
        cls,
        *,
        start_nm: float,
        stop_nm: float | None = None,
        steps: int = 1,
        damping: float | None = None,
        maxcycle: int = 200,
        tolalpha: int = 6,
    ) -> "CphfSettings":
        """Create settings for a frequency-dependent CPHF/CPKS calculation."""
        wavelengths = (float(start_nm),) if stop_nm is None else (float(start_nm), float(stop_nm))
        return cls(
            wavelengths_nm=wavelengths,
            steps=int(steps),
            damping=damping,
            maxcycle=maxcycle,
            tolalpha=tolalpha,
        )


def runpcry_command(cores: int, job_name: str) -> str:
    """Return the MPI command used to run a CRYSTAL23 wavefunction job."""
    return f"runPcry23 {cores} {job_name}"


def runpprop_command(cores: int, property_job: str, wavefunction_job: str) -> str:
    """Return the MPI command used to run a CRYSTAL23 properties job."""
    return f"runPprop23 {cores} {property_job} {wavefunction_job}"


def band_input(case_name: str, path: BandPath | None = None) -> str:
    """Build a CRYSTAL ``BAND`` property input for one calculation case."""
    path = path or BandPath()
    lines = [
        "BAND",
        f"BANDS_{case_name}"[:72],
        f"{path.n_segments} {path.denominator} {path.subdivisions} "
        f"{path.first_band} {path.last_band} {path.iplo} {path.lpr66}",
    ]
    for (_, start), (_, stop) in zip(path.points[:-1], path.points[1:]):
        lines.append("{} {} {}   {} {} {}".format(*start, *stop))
    lines.append("END")
    return "\n".join(lines) + "\n"


def doss_input(
    case_name: str,
    *,
    shrink: tuple[int, int, int] = (6, 12, 1),
    first_band: int = 1,
    last_band: int = 180,
    n_points: int = 2000,
) -> str:
    """Build a CRYSTAL ``DOSS`` property input."""
    return text_block(
        "NEWK",
        f"{shrink[0]} {shrink[1]} {shrink[2]}",
        "DOSS",
        f"0 {n_points} {first_band} {last_band} 2 14 0",
        "END",
    )


def anbd_input(case_name: str, *, shrink: tuple[int, int, int] = (6, 12, 1)) -> str:
    """Build a compact ``ANBD`` property input."""
    return text_block("NEWK", f"{shrink[0]} {shrink[1]} {shrink[2]}", "ANBD", "END")


def echg_input(case_name: str, *, npx: int = 200, npy: int = 200) -> str:
    """Build a 2D charge/spin-density map input using ``ECHG``."""
    return text_block("ECHG", "0", f"{npx} {npy}", "0.0 0.0 0.0", "1.0 0.0 0.0", "0.0 1.0 0.0", "END")


def ech3_input(case_name: str, *, n_points: int = 160) -> str:
    """Build a 3D charge/spin-density cube input using ``ECH3``."""
    return text_block("ECH3", str(n_points), "END")


def ppan_input(case_name: str) -> str:
    """Build a Mulliken population-analysis property input using ``PPAN``."""
    return text_block("PPAN", "END")


def property_input(kind: PropertyKind | str, case_name: str) -> str:
    """Build one property-input text block for a named CRYSTAL case."""
    key = kind.lower()
    if key in {"band", "bands"}:
        return band_input(case_name)
    if key in {"dos", "doss"}:
        return doss_input(case_name)
    if key == "anbd":
        return anbd_input(case_name)
    if key == "echg":
        return echg_input(case_name)
    if key == "ech3":
        return ech3_input(case_name)
    if key == "ppan":
        return ppan_input(case_name)
    raise ValueError(f"Unknown CRYSTAL property kind: {kind!r}")


def cphf_block(settings: CphfSettings | None = None) -> str:
    """Return a CRYSTAL ``CPHF`` block."""
    settings = settings or CphfSettings.static()
    lines = ["CPHF"]
    if settings.wavelengths_nm is not None:
        lines.append("DYNAMIC")
        lines.append(str(settings.steps))
        lines.extend(f"{value:.10g}" for value in settings.wavelengths_nm)
        if settings.damping is not None:
            lines.append("DAMPING")
            lines.append(f"{settings.damping:.10g}")
    lines.extend(["MAXCYCLE", str(settings.maxcycle), "TOLALPHA", str(settings.tolalpha), "END"])
    return "\n".join(lines) + "\n"


def insert_cphf_block(input_text: str, settings: CphfSettings | None = None) -> str:
    """Insert a CPHF block before the first geometry ``END`` in a full input."""
    lines = input_text.rstrip().splitlines()
    geometry_start = _geometry_start_index(lines)
    insert_at = _first_geometry_end_index(lines, geometry_start)
    block_lines = cphf_block(settings).rstrip().splitlines()
    return "\n".join([*lines[:insert_at], *block_lines, *lines[insert_at:]]) + "\n"


def text_block(*lines: str) -> str:
    """Join CRYSTAL input lines and terminate with a newline."""
    return "\n".join(lines) + "\n"


def _geometry_start_index(lines: list[str]) -> int:
    for index, line in enumerate(lines):
        if line.strip().upper() in {"CRYSTAL", "SLAB", "POLYMER", "MOLECULE", "EXTERNAL"}:
            return index
    raise ValueError("No CRYSTAL geometry declaration found in input text.")


def _first_geometry_end_index(lines: list[str], start: int) -> int:
    for index in range(start + 1, len(lines)):
        if lines[index].strip().upper() == "END":
            return index
    raise ValueError("No geometry-closing END record found in input text.")
