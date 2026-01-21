# <img src="https://github.com/user-attachments/assets/09b208df-02c0-48e6-8013-882ceb73d628" alt="β-Ga₂O₃ structure" width="64"/> **β‑Ga₂O₃ defect dichroism post‑processing toolkit** 

This repository bundles two things in one place:

1) **Research data** (CRYSTAL23 inputs/outputs + curated post‑processing files) for pristine and defective **β‑Ga₂O₃**, including **Ga vacancies** (tetrahedral / octahedral), **vacancy‑bound self-trapped holes**, **oxygen vacancies** (O₁/O₂/O₃), and a finite‑size **1×2×2** check.

2) A lightweight, testable **Python library + CLI** named **`ga2o3`** to **parse** and **plot** the key files produced in this workflow:
- CRYSTAL `.BAND` (spin bands)
- CRYSTAL `.DOSS` (spin/projection DOS)
- CRYSTAL CPHF/CPKS `.out` (dielectric tensor → absorption, polarization‑resolved)
- CASTEP `.f25` (2D maps for charge/spin density visualization)

## ✨ What `ga2o3` gives you

- **One command** (`ga2o3`) to generate consistent, publication‑style plots.
- **Dataset-aware path resolution**: it understands your on‑disk layout (cases like `pristine`, `tetra`, `octa`, `O_1`, `tetra_1x2x2`, …).
- **Robust parsers** for common CRYSTAL/CASTEP text formats (permissive readers; avoids brittle regex-only parsing).
- **Headless-safe plotting** (matplotlib `Agg` backend by default → works on clusters/CI).
- **Tests included** (`pytest`) with fixture files to prevent regressions when you refactor.

---

## 🗺 Repository map

High-level layout (key parts):

```text
.
├── data/                  # normalized/curated files for plotting (recommended GA2O3_DATA root)
│   ├── BANDS/             # per-case folders, each contains ALPHA.BAND / BETA.BAND
│   ├── DOS/               # flat folder with .DOSS files (case encoded in filename)
│   ├── CPKS_CPHF/          # per-case folders, each contains many discrete CPHF .out files
│   └── formation_energy/  # energies spreadsheet + outputs (parsing not implemented yet)
├── crystal23/             # raw CRYSTAL23 inputs/outputs (full provenance)
├── 3DPlots/               # cube/cmdx visualization assets (VESTA/VMD/etc.)
├── src/ga2o3/             # the Python package
└── tests/                 # pytest suite + small sample data files
```

**Important:** the `ga2o3.dataset.Ga2O3Dataset` class expects a root containing:
`BANDS/`, `DOS/`, `CPKS_CPHF/`, `ANBD/`, `formation_energy/`.  
In this repo, that root is **`./data`**, so you almost always want:

```bash
export GA2O3_DATA="$(pwd)/data"
```

---

## 🚀 Installation

### Requirements
- Python **≥ 3.9**
- numpy, pandas, matplotlib (installed automatically)

### Install (editable, recommended during development)
From the repository root (folder containing `pyproject.toml`):

```bash
python -m pip install -U pip
python -m pip install -e ".[dev]"
```

### Run tests
```bash
pytest -q
```

---

## ⚡ Quickstart (CLI)

### 1) Point the toolkit to the dataset
```bash
export GA2O3_DATA="$(pwd)/data"
```

### 2) Bands (spin-polarized)
```bash
ga2o3 bands --case pristine --out figs/bands_pristine.png --emin -1 --emax 5
ga2o3 bands --case tetra    --out figs/bands_tetra.png    --emin -1 --emax 5
ga2o3 bands --case octa     --out figs/bands_octa.png     --emin -1 --emax 5
```

Optional:
- `--shift <eV>` subtracts a rigid energy offset (useful if you precomputed VBM shift).
- `--legend` shows an alpha/beta legend.

### 3) DOS (spin-resolved; auto-finds ALPHA/BETA by filename)
```bash
ga2o3 dos --case tetra --out figs/dos_tetra.png --emin -0.5 --emax 0 --xmin -1000 --xmax 1000
ga2o3 dos --case octa  --out figs/dos_octa.png  --emin -0.5 --emax 0 --xmin -1000 --xmax 1000
```

If you want explicit input files:
```bash
ga2o3 dos \
  --case tetra \
  --alpha data/DOS/GaO_TETRA_18HF_DOS_dat_ALPHA.DOSS \
  --beta  data/DOS/GaO_TETRA_18HF_DOS_dat_BETA.DOSS \
  --out   figs/dos_tetra.png
```

### 4) Polarization-resolved absorption from CPHF outputs
This reads **all** `*.out` files in `data/CPKS_CPHF/<case>/` and overlays curves.

```bash
ga2o3 cphf --outdir figs/cphf --yscale log --correction ratio --wl-min 465
```

Defaults are tuned to your workflow:
- `--cases Pristine TETRA OCTA`
- energy correction modes:
  - `none`   (no correction)
  - `scissor` (E' = E + Δ)
  - `ratio`  (E' = E * (Eg_target/Eg_pbe) + extra_shift)
- default gaps: `Eg_pbe=2.62 eV`, `Eg_target=4.43 eV`

Outputs written to `--outdir`:
- `absorption_iso.png`
- `absorption_xx.png`
- `absorption_yy.png`
- `absorption_zz.png`

### 5) Charge density map from CASTEP `.f25`
This command takes a direct file path (it does not require `GA2O3_DATA`):

```bash
ga2o3 density crystal23/ECHG/output/Ga2O3_PRISTINE_ECHG.f25 --out figs/charge_pristine.png
```

Percentile clipping controls contrast:
```bash
ga2o3 density crystal23/ECHG/output/Ga2O3_2x2x2_VAC_TETRA_ECHG.f25 --out figs/charge_tetra.png --pmin 0.01 --pmax 99.2
```

### CLI help
```bash
ga2o3 --help
ga2o3 bands  --help
ga2o3 dos    --help
ga2o3 cphf   --help
ga2o3 density --help
```

---

## 🧠 Using the library in Python

### Set the dataset root
```python
from ga2o3.dataset import Ga2O3Dataset

ds = Ga2O3Dataset.from_env()   # reads $GA2O3_DATA
# ds = Ga2O3Dataset(root="data")  # equivalent (explicit)
```

### Bands
```python
from ga2o3.plot.bands import plot_bands_only

alpha, beta = ds.band_files("pristine")
fig, ax = plot_bands_only(alpha, beta, e_min=-1, e_max=5, band_legend=False)
fig.savefig("figs/bands_pristine.png", dpi=600, bbox_inches="tight")
```

### DOS
```python
from ga2o3.plot.dos import plot_spin_resolved_dos

alpha, beta = ds.dos_spin_files("tetra")
fig, ax = plot_spin_resolved_dos(alpha, beta, e_range=(-0.5, 0.0), dos_range=(-1000, 1000))
fig.savefig("figs/dos_tetra.png", dpi=600, bbox_inches="tight")
```

### CPHF absorption overlays
```python
from ga2o3.plot.cphf import plot_absorption_overlays, save_absorption_figs

figs = plot_absorption_overlays(
    ds,
    cases=["Pristine", "TETRA", "OCTA"],
    wl_min=465.0,
    y_scale="log",
    correction="ratio",   # "none" | "scissor" | "ratio"
)
save_absorption_figs(figs, "figs/cphf", fmt="png", dpi=600)
```

### Density maps from `.f25`
```python
from ga2o3.plot.density import plot_charge_density_f25

plot_charge_density_f25(
    "crystal23/ECHG/output/Ga2O3_PRISTINE_ECHG.f25",
    out_path="figs/charge_pristine.png",
    charge_percentiles=(0.01, 99.2),
    dpi=900,
)
```

---

## 📁 Expected data layout (what the parsers assume)

`ga2o3` is designed around this layout.

```text
<GA2O3_DATA>/
  BANDS/
    pristine/
      ALPHA.band   (or ALPHA.BAND)
      BETA.band    (or BETA.BAND)
      ...
    tetra/
      ALPHA.BAND
      BETA.BAND
      ...
    octa/
      ALPHA.BAND
      BETA.BAND
      ...
    O_1/
      ALPHA.BAND
      BETA.BAND
      ...
    tetra_1x2x2/
    octa_1x2x2/
  DOS/
    *ALPHA*.DOSS
    *BETA*.DOSS
    ...
  CPKS_CPHF/
    pristine/
      *.out
    tetra/
      *.out
    octa/
      *.out
```

### Case name matching
Case names are matched **case-insensitively** and accept a few aliases:

- `bulk` → `pristine`
- `tetrahedral` → `tetra`
- `octahedral` → `octa`
- `o1`/`o_1` → `O_1` (and similarly `O_2`, `O_3`)
- `tetra1x2x2` → `tetra_1x2x2` (and similarly octa)

---

## 🔬 Scientific background (what the data represent)

### System and supercells
- Host: monoclinic **β‑Ga₂O₃** (space group **C2/m**).
- Lattice parameters fixed to experiment (as used in the production calculations):
  - a = 12.215 Å  
  - b = 3.050 Å  
  - c = 5.820 Å  
  - β = 103.978°
- Production supercell: **2×2×2** (160 atoms).
- Finite-size check: **1×2×2** (80 atoms).

### Ground-state electronic structure
- **Spin‑polarized hybrid DFT** (HSE‑type screened hybrid) for:
  - pristine band gap benchmark
  - Ga vacancy states
  - vacancy‑bound O‑centered self-trapped holes
- Typical SCF convergence: **1×10⁻⁷ Hartree** (as set in the CRYSTAL inputs).

### Band-structure path
Band structures follow the monoclinic high‑symmetry path:

```text
X1–Y–Γ–N–X–Γ–M–I–L–F–Y–Γ–Z–F1–Z–I1
```

Energies are usually plotted with **VBM set to 0 eV** (or a consistent reference shift).

### Self-trapped holes fingerprint criteria
A vacancy-bound O-centered self-trapped hole is identified by both:

1) **Real-space localization:** magnetization density
   \[
   m(\mathbf{r})=\rho_\alpha(\mathbf{r})-\rho_\beta(\mathbf{r})
   \]
   localized on an O neighbor (O‑2p character).

2) **Reciprocal-space signature:** a nearly dispersionless, spin‑polarized **in‑gap (flat) band** near the VBM.

### Optical response workflow (CPHF/CPKS)
- Polarization-resolved clamped-ion dielectric response is computed using periodic **CPHF/CPKS**.
- Practical constraint: hybrid functionals are not used for this step; optics are computed at **PBE** level.
- The photon-energy axis is aligned to the hybrid gap using a rigid **scissor** or multiplicative **ratio** correction (implemented in `ga2o3` via `--correction`).

### Supplementary datasets included
- Neutral oxygen vacancies at inequivalent sites: **V_OI / V_OII / V_OIII** (cases `O_1`, `O_2`, `O_3`).
- Finite-size robustness check: `tetra_1x2x2`, `octa_1x2x2`.

---

## 📚 References (methods/software)
- **CRYSTAL23**: A. Erba et al., *J. Chem. Theory Comput.* **19**, 6891 (2023).
- Screened hybrid functionals (HSE family): Heyd–Scuseria–Ernzerhof.
- Supercell defect formalism: Freysoldt, Neugebauer, Van de Walle, *Rev. Mod. Phys.* **86**, 253 (2014).
- Periodic CPHF/CPKS dielectric response in CRYSTAL: Ferrero / Lacivita / co-workers (see CRYSTAL documentation and related papers).

---

## 🧩 Scope / limitations (so expectations are clear)

- `ga2o3` is **not** a full CRYSTAL output parser; it focuses on the file types above.
- `formation_energy/` is included for provenance; **formation-energy parsing is not implemented yet** in this package.
- CASTEP `.f25` parsing assumes the presence of `MAPN` blocks and reads up to two maps by default.

---

## 📜 License
MIT License. See `LICENSE`.

---

## 🙏 Acknowledgments
CRYSTAL23 developers and maintainers for the electronic-structure engine and periodic CPHF/CPKS implementation.
Universidad Complutense de Madrid
