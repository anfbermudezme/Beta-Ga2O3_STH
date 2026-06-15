"""Create a clean CRYSTAL23 run tree from the bundled example inputs."""

from pathlib import Path
from ga2o3dichroism import stage_crystal_workflow

stage_crystal_workflow(input_dir=Path("examples/crystal23/scf"), out_dir=Path("runs"), cores=128)
