# Reproducibility design

This repository is organized so that the CRYSTAL23 workflow can be regenerated from a small, declarative workflow file instead of manually copying and renaming inputs.

## Core principle

The source `.d12` templates remain in `examples/crystal23/scf/`. The command

```bash
ga2o3d generate --workflow examples/workflows/publication.yml --input-dir examples/crystal23/scf --out runs --cores 128
```

creates all reproducibility inputs required to rerun the computational protocol:

- SCF / geometry `.d12` inputs,
- CPHF/CPKS full `.d12` inputs,
- post-SCF property `.d3` inputs for `BANDS`, `DOSS`, `ANBD`, `ECHG`, `ECH3`, and `PPAN`,
- executable shell scripts with the CRYSTAL23 MPI commands,
- `manifest.csv` with file paths, job names, commands, and descriptions.

## Why this is extensible

New defects, charge states, functionals, optical-response windows, or property calculations can be added by editing `examples/workflows/publication.yml`. The Python package does not need to be modified for ordinary workflow extensions.

This makes the workflow auditable: the manuscript calculations are represented as structured data, while the code only implements reusable CRYSTAL23 input builders.

## What is intentionally not stored

Large raw CRYSTAL outputs, wavefunction/restart files, and generated figures should remain outside the Python package. They can be regenerated locally under `runs/`, `data/`, and `figures/`.
