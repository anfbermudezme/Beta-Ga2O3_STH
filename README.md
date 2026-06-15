# <img src="https://github.com/user-attachments/assets/09b208df-02c0-48e6-8013-882ceb73d628" alt="beta-Ga2O3 structure" width="72"/> Ga₂O₃ defect dichroism

<p>
  <a href="LICENSE"><img alt="License: MIT" src="https://img.shields.io/badge/license-MIT-blue.svg"></a>
  <img alt="Version: 0.1.1" src="https://img.shields.io/badge/version-0.1.1-success.svg">
  <a href="https://www.researchsquare.com/article/rs-9008105/v1"><img alt="Manuscript status: under review" src="https://img.shields.io/badge/manuscript-under%20review-orange.svg"></a>
  <a href="https://github.com/anfbermudezme/ga2o3-defect-dichroism"><img alt="Reproducible CRYSTAL23 workflow" src="https://img.shields.io/badge/CRYSTAL23-input%20generator-6f42c1.svg"></a>
</p>

Readable tools for preparing **CRYSTAL23** inputs and reproducing the Python analysis used for the β-Ga₂O₃ vacancy-bound self-trapped-hole / dichroism study.

The associated manuscript is currently **under review**. A public Research Square version is available here: <https://www.researchsquare.com/article/rs-9008105/v1>.

## License

This repository is released under the **MIT License**. See [`LICENSE`](LICENSE).

The MIT License applies to the software source code and documentation in this repository, unless a file states otherwise. Scientific conclusions, manuscript text, figures, raw simulation outputs, third-party datasets, and external software such as CRYSTAL23 remain under their own terms. See [`NOTICE`](NOTICE) for the project notice and citation request.

## How to use it

### Install

```bash
python -m pip install -U pip
python -m pip install -e ".[dev]"
pytest -q

# Optional release tooling check
bump-my-version show current_version
```

The command-line entry point is:

```bash
ga2o3d --help
```

### 1. Generate the complete CRYSTAL23 input tree

```bash
ga2o3d generate \
  --workflow examples/workflows/publication.yml \
  --input-dir examples/crystal23/scf \
  --out runs \
  --cores 128
```

This creates a case-centered run tree:

```text
runs/
├── crystal23/
│   ├── GaO_TETRA_VAC/
│   │   ├── GaO_TETRA_VAC.d12
│   │   ├── CPHF_GaO_TETRA_VAC_STATIC.d12
│   │   ├── CPHF_GaO_TETRA_VAC_7500_12500.d12
│   │   ├── BANDS_GaO_TETRA_VAC.d3
│   │   ├── DOSS_GaO_TETRA_VAC.d3
│   │   ├── ANBD_GaO_TETRA_VAC.d3
│   │   ├── ECHG_GaO_TETRA_VAC.d3
│   │   ├── ECH3_GaO_TETRA_VAC.d3
│   │   ├── PPAN_GaO_TETRA_VAC.d3
│   │   └── run_case.sh
│   └── GaO_OCTA_VAC/
├── scripts/
│   ├── run_scf.sh
│   ├── run_cphf.sh
│   ├── run_properties.sh
│   └── run_all.sh
└── manifest.csv
```

The generated `manifest.csv` records every generated input, the CRYSTAL job name, and the command used to run it. See [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md) for the workflow-design rationale.

### 2. Run the SCF / geometry jobs

Run a single case:

```bash
cd runs/crystal23/GaO_TETRA_VAC
runPcry23 128 GaO_TETRA_VAC
```

Or run all generated SCF jobs:

```bash
bash runs/scripts/run_scf.sh
```

### 3. Run post-SCF properties

```bash
cd runs/crystal23/GaO_TETRA_VAC
runPprop23 128 BANDS_GaO_TETRA_VAC GaO_TETRA_VAC
runPprop23 128 DOSS_GaO_TETRA_VAC  GaO_TETRA_VAC
runPprop23 128 ANBD_GaO_TETRA_VAC  GaO_TETRA_VAC
runPprop23 128 ECHG_GaO_TETRA_VAC  GaO_TETRA_VAC
runPprop23 128 ECH3_GaO_TETRA_VAC  GaO_TETRA_VAC
runPprop23 128 PPAN_GaO_TETRA_VAC  GaO_TETRA_VAC
```

Or run all generated property jobs:

```bash
bash runs/scripts/run_properties.sh
```

### 4. Run CPHF/CPKS inputs

Static and dynamic CPHF/CPKS `.d12` files are generated automatically for the main Ga-vacancy cases defined in `examples/workflows/publication.yml`.

```bash
cd runs/crystal23/GaO_TETRA_VAC
runPcry23 128 CPHF_GaO_TETRA_VAC_STATIC
runPcry23 128 CPHF_GaO_TETRA_VAC_7500_12500
```

Or run all generated CPHF/CPKS jobs:

```bash
bash runs/scripts/run_cphf.sh
```

### 5. Plot or parse results

Bands:

```bash
ga2o3d plot-bands --input data/bands/GaO_TETRA_VAC.BAND --out figures/bands_tetra.png
```

DOS:

```bash
ga2o3d plot-dos --input data/dos/GaO_TETRA_VAC.DOSS --out figures/dos_tetra.png
```

CPHF/CPKS output folder to CSV:

```bash
ga2o3d parse-cphf --folder data/cphf/tetra --out data/processed/cphf_tetra.csv
```

Mulliken oxygen-sublattice spin summary:

```bash
ga2o3d plot-mulliken --out figures/mulliken_spin.png
```

## Extending the workflow

Add new calculations by editing `examples/workflows/publication.yml`; the Python code does not need to be modified.

A minimal new case looks like this:

```yaml
- name: GaO_NEW_DEFECT
  source: Ga2O3_NEW_DEFECT.d12
  label: short description of the defect model
  properties:
    - bands
    - doss
    - anbd
    - echg
    - ech3
    - ppan
  cphf:
    - mode: dynamic
      name: CPHF_GaO_NEW_DEFECT_7500_12500
      start_nm: 7500
      stop_nm: 12500
      steps: 10
      damping: 0.002
```

The generator will automatically create the SCF `.d12`, CPHF/CPKS `.d12`, post-SCF `.d3` files, run scripts, and manifest entries.

## Versioning

This project uses [`bump-my-version`](https://callowayproject.github.io/bump-my-version/) for release version updates. The current package version is **0.1.1**.

Show the active version and possible bump targets:

```bash
bump-my-version show current_version
bump-my-version show-bump
```

Create a release bump from a clean Git working tree:

```bash
bump-my-version bump patch   # 0.1.1 -> next patch release
bump-my-version bump minor   # 0.1.1 -> next minor release
bump-my-version bump major   # 0.1.1 -> next major release
```

The configured bump updates the synchronized version strings in `pyproject.toml`, `src/ga2o3dichroism/__init__.py`, `README.md`, `CITATION.cff`, `NOTICE`, `VERSIONING.md`, and `package-summary.json`, then creates a release commit and a `v<version>` Git tag.

## Python examples

Create a CRYSTAL property input:

```python
from ga2o3dichroism.crystal import property_input

text = property_input("bands", case_name="GaO_TETRA_VAC")
print(text)
```

Generate the full input tree from Python:

```python
from ga2o3dichroism.generate import generate_crystal_inputs, load_workflow

cases = load_workflow("examples/workflows/publication.yml")
result = generate_crystal_inputs(
    input_dir="examples/crystal23/scf",
    out_dir="runs",
    cores=128,
    cases=cases,
)
print(result.manifest)
```

Read and plot a DOSS file:

```python
from ga2o3dichroism.parsers import read_doss
from ga2o3dichroism.plotting import plot_dos

dos = read_doss("data/dos/GaO_TETRA_VAC.DOSS")
fig, ax = plot_dos(dos, energy_window=(-1, 5))
fig.savefig("figures/dos_tetra.png", dpi=600, bbox_inches="tight")
```

Keep large raw CRYSTAL outputs, wavefunction files, and exported figures outside the Python package. Use local `data/`, `runs/`, and `figures/` folders.
