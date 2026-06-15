# <img src="https://github.com/user-attachments/assets/09b208df-02c0-48e6-8013-882ceb73d628" alt="beta-Ga2O3 structure" width="72"/> Ga₂O₃ defect dichroism

<p>
  <a href="LICENSE"><img alt="License: MIT" src="https://img.shields.io/badge/license-MIT-blue.svg"></a>
  <a href="https://www.researchsquare.com/article/rs-9008105/v1"><img alt="Manuscript status: under review" src="https://img.shields.io/badge/manuscript-under%20review-orange.svg"></a>
  <a href="https://github.com/anfbermudezme/ga2o3-defect-dichroism"><img alt="Reproducible CRYSTAL23 workflow" src="https://img.shields.io/badge/CRYSTAL23-workflow-6f42c1.svg"></a>
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
```

The command-line entry point is:

```bash
ga2o3d --help
```

### 1. Create a clean CRYSTAL23 run folder

```bash
ga2o3d stage --input-dir examples/crystal23/scf --out runs --cores 128
```

This copies the supplied `.d12` inputs into `runs/scf`, creates property inputs in `runs/properties`, and writes simple shell scripts in `runs/scripts`.

### 2. Run the SCF / geometry jobs

```bash
cd runs/scf
runPcry23 128 GaO_TETRA_VAC
runPcry23 128 GaO_OCTA_VAC
```

### 3. Run post-SCF properties

```bash
cd runs/properties
runPprop23 128 BANDS_GaO_TETRA_VAC GaO_TETRA_VAC
runPprop23 128 DOSS_GaO_TETRA_VAC  GaO_TETRA_VAC
runPprop23 128 ANBD_GaO_TETRA_VAC  GaO_TETRA_VAC
runPprop23 128 ECHG_GaO_TETRA_VAC  GaO_TETRA_VAC
runPprop23 128 ECH3_GaO_TETRA_VAC  GaO_TETRA_VAC
```

### 4. Generate a CPHF/CPKS full input

CPHF/CPKS is handled as a full `.d12` input by inserting the `CPHF ... END` block in the geometry section.

Static response:

```bash
ga2o3d make-cphf examples/crystal23/scf/Ga2O3_VAC_TETRA_2x2x2.d12 \
  --out runs/scf/CPHF_GaO_TETRA_VAC.d12
```

Dynamic response from 7500 to 12500 nm:

```bash
ga2o3d make-cphf examples/crystal23/scf/Ga2O3_VAC_TETRA_2x2x2.d12 \
  --out runs/scf/CPHF_GaO_TETRA_VAC_7500_12500.d12 \
  --dynamic --steps 10 --start-nm 7500 --stop-nm 12500 --damping 0.002
```

Then run it as a normal CRYSTAL calculation:

```bash
cd runs/scf
runPcry23 128 CPHF_GaO_TETRA_VAC_7500_12500
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

## Python examples

Create a CRYSTAL property input:

```python
from ga2o3dichroism.crystal import property_input

text = property_input("bands", case_name="GaO_TETRA_VAC")
print(text)
```

Generate a full dynamic CPHF input:

```python
from pathlib import Path
from ga2o3dichroism.crystal import CphfSettings, insert_cphf_block

scf_text = Path("examples/crystal23/scf/Ga2O3_VAC_TETRA_2x2x2.d12").read_text()
settings = CphfSettings.dynamic(start_nm=7500, stop_nm=12500, steps=10, damping=0.002)
Path("CPHF_GaO_TETRA_VAC.d12").write_text(insert_cphf_block(scf_text, settings))
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
