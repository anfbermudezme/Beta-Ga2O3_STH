"""Generate a dynamic CPHF/CPKS input from an existing CRYSTAL .d12 file."""

from pathlib import Path
from ga2o3dichroism.crystal import CphfSettings, insert_cphf_block

source = Path("examples/crystal23/scf/Ga2O3_VAC_TETRA_2x2x2.d12")
settings = CphfSettings.dynamic(start_nm=7500, stop_nm=12500, steps=10, damping=0.002)
text = insert_cphf_block(source.read_text(), settings)
Path("runs/scf").mkdir(parents=True, exist_ok=True)
Path("runs/scf/CPHF_GaO_TETRA_VAC_7500_12500.d12").write_text(text)
