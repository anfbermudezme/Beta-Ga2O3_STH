# CRYSTAL23 calculation package — β-Ga₂O₃ Ga vacancies, vacancy-bound hole polarons, and polarization fingerprints

This repository contains **CRYSTAL23** input files, output logs, and post-processing code used to generate the *first-principles calculations* reported in:

- **Main paper:** “Site-Selective Small Polarons induced by gallium vacancies and their polarization fingerprints in β-Ga₂O₃”.
- **Supplementary material:** oxygen-vacancy comparison and finite-size (1×2×2) robustness checks.

The workflows implemented here cover:
1) Pristine β-Ga₂O₃ benchmark (hybrid-DFT band structure / gap).
2) Gallium vacancies on inequivalent sublattices: **V_GaI** (tetrahedral) and **V_GaII** (octahedral).
3) Vacancy-bound **oxygen-centered small-hole polarons**: real-space spin localization + in-gap flat bands.
4) Polarization-resolved optical response via **PBE-CPKS** (clamped-ion dielectric tensor), aligned to the hybrid gap.
5) Supplementary checks:
   - Neutral oxygen vacancies **V_OI / V_OII / V_OIII** (non-magnetic, deep-donor mid-gap states).
   - Ga-vacancy polaron states in a smaller **1×2×2** supercell (finite-size / higher concentration test).

---

## 1. Code, data, and reproducibility

### 1.1 Software
- **CRYSTAL23** (periodic Kohn–Sham with localized Gaussian basis sets; main engine used throughout).
- Optional visualization/post-processing tools (depending on your workflow):
  - A plotting environment (Python/matplotlib, gnuplot, etc.)
  - A structure/spin-density viewer (e.g., VESTA) if you export grids / cube-like files.

### 1.2 What is included
- CRYSTAL23 inputs for bulk and defect supercells.
- CRYSTAL23 outputs used to extract:
  - total energies, band structures, DOS/PDOS
  - spin density / magnetization signatures of polarons
  - dielectric response from CPKS
- Scripts used to parse energies and generate plots (if provided in this repository).

### 1.3 What you may need to edit to run on your system
- CRYSTAL23 executable path / MPI launcher.
- Filenames and scratch directories.
- Any cluster-specific batch scripts (SLURM/PBS).

---

## 2. Directory structure

- 1x2x2/                 # Finite-size test supercell calculations (Supplementary Information)
- 3D/                    # 3D spin-density and real-space grid data
- ANBD/                  # Repository-specific ANBD analysis (interpret as needed)
- BandsAndDOS_code/      # Scripts for extracting and plotting band structures and DOS
- ECHG_code/             # Scripts for ECHG, charge-density, and spin-density processing
- FormationEnergy/       # Scripts to parse outputs and compute defect formation energies
- HSE-Type_SCF/          # Hybrid-DFT (HSE-type) SCF input and output files
- HSE06_full_relaxed/    # Fully relaxed HSE06 structures and related outputs
- figures/               # Source files for figures used in the manuscript
- CPKS/                  # CPKS optical results
- README.md              # Project README

---

## 3. Physical model and key numerical settings

### 3.1 Structure and supercells
- Host: monoclinic **β-Ga₂O₃**, space group **C2/m**.
- Lattice parameters fixed to experiment:
  - a = 12.215 Å
  - b = 3.050 Å
  - c = 5.820 Å
  - β = 103.978°
- Main production supercell: **2×2×2** (160 atoms).
  - One vacancy per 160 atoms ≈ 0.6 at.%.
  - One Ga vacancy per 64 Ga sites ≈ 1.6% on the Ga sublattice (ordered-array artifact; used as a dilute-limit compromise).
- Supplementary finite-size test: **1×2×2** (80 atoms).
  - One Ga vacancy per 32 Ga sites ≈ 3.1% on the Ga sublattice (stronger periodic-image interactions).

### 3.2 Electronic structure method (ground state)
- **Spin-polarized hybrid DFT** using an **HSE-type screened hybrid functional** (Heyd–Scuseria–Ernzerhof family).
- The screening parameter and exact-exchange fraction are set in the CRYSTAL23 input so that:
  - the pristine β-Ga₂O₃ band gap matches optical benchmarks, and
  - O-centered small polarons are stabilized with realistic trapping behavior.
- SCF electronic convergence threshold: **1×10⁻⁷ Hartree**.
- Brillouin-zone sampling: **Monkhorst–Pack** meshes with consistent k-point density across different supercells.

### 3.3 Band-structure path
Band structures were computed along the monoclinic high-symmetry path:

    X1–Y–Γ–N–X–Γ–M–I–L–F–Y–Γ–Z–F1–Z–I1

Energies are typically plotted with **VBM set to 0 eV**.

### 3.4 Polaron identification (how we decide a “hole polaron exists”)
A vacancy-bound small-hole polaron is identified by both:
1) **Real-space signature:** localized magnetization density
       m(r) = ρ_α(r) − ρ_β(r)
   concentrated on a specific oxygen neighbor (O-2p-like localization).
2) **Reciprocal-space signature:** a nearly dispersionless (flat), spin-polarized **mid-gap band** just above the VBM.

Energy quantities used in the analysis:
- Polaron level position: ΔE = (VBM → center of polaron flat band)
- Polaron-to-CBM excitation: ΔE = (polaron flat band → CBM)

## 4. Optical response workflow (main paper)

### 4.1 Method: CPKS dielectric tensor (clamped-ion)
Goal:
- Compute polarization-resolved near-edge optical response trends for pristine, V_GaI, and V_GaII.

Method:
- Use CRYSTAL23 periodic **CPHF/CPKS** implementation to compute the purely electronic (clamped-ion) dielectric tensor ε(ω).

Important constraint:
- The current CPKS implementation does not support hybrid functionals.
- Therefore, dielectric response is computed at **PBE** level.

Polarizations:
- E ∥ a
- E ∥ b
- E ∥ c
- Orientational average:
    ε(ω) = (1/3) [ε_xx(ω) + ε_yy(ω) + ε_zz(ω)]

Energy alignment:
- Apply an a posteriori rigid shift to the photon-energy axis so the PBE absorption edge aligns to the hybrid-DFT band gap.

Qualitative validation trends:
- For E ∥ b: V_GaII absorbance > V_GaI across the shown UV range.
- For E ∥ c: V_GaI dominates and can exceed V_GaII by ~5× at long wavelengths.

---

## 5. Supplementary calculations included in this repository

### 5.1 Oxygen vacancies: V_OI, V_OII, V_OIII (neutral)
What was computed:
- Spin-polarized calculations for neutral oxygen vacancies at the three inequivalent oxygen sites:
  - V_OI, V_OII, V_OIII

Key result (robust check):
- All oxygen-vacancy configurations converge to a **spin-degenerate (non-magnetic)** ground state:
  - magnetization m(r) ≈ 0 (no exchange splitting)
  - oxygen vacancies alone do **not** stabilize oxygen-centered small-hole polarons in the conditions considered

Electronic-structure fingerprint:
- Oxygen vacancies introduce occupied, weakly dispersive **mid-gap donor bands**.
- Donor-level energies depend on which oxygen sublattice hosts the vacancy, implying a distribution of donor-to-CBM photoionization thresholds and sub-gap optical responses.

### 5.2 Finite-size / higher-concentration check: 1×2×2 Ga-vacancy supercell
What was computed:
- Repeat V_GaI and V_GaII calculations in a smaller **1×2×2** supercell (80 atoms).

Why this matters:
- Smaller supercells increase periodic image interactions; polaronic in-gap levels can broaden and shift.

---

## 6. References (methods)
- A. Erba et al., “CRYSTAL23: A program for computational solid state physics and chemistry”, J. Chem. Theory Comput. 19, 6891 (2023).
- HSE screened hybrid functional family (Heyd–Scuseria–Ernzerhof).
- Supercell defect formalism reference (Freysoldt, Neugebauer, Van de Walle, Rev. Mod. Phys. 86, 253 (2014)).
- CPKS/CPHF periodic dielectric response in CRYSTAL (Ferrero et al., Lacivita et al.).
