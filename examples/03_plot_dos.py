"""Plot a CRYSTAL DOSS-like file after it has been copied into data/dos."""

from pathlib import Path
from ga2o3dichroism.parsers import read_doss
from ga2o3dichroism.plotting import plot_dos

dos = read_doss("data/dos/GaO_TETRA_VAC.DOSS", energy_unit="hartree")
fig, _ = plot_dos(dos, energy_window=(-1, 5))
Path("figures").mkdir(parents=True, exist_ok=True)
fig.savefig("figures/dos_tetra.png", dpi=600, bbox_inches="tight")
