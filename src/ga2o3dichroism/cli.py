"""Command-line interface for the Ga₂O₃ defect-dichroism helpers."""

from __future__ import annotations

import argparse
from pathlib import Path

from ga2o3dichroism.crystal import CphfSettings, insert_cphf_block, property_input
from ga2o3dichroism.parsers import parse_cphf_folder, read_band, read_doss
from ga2o3dichroism.plotting import plot_bands, plot_dos, plot_mulliken_summary
from ga2o3dichroism.workflow import default_cases, stage_crystal_workflow


def main(argv: list[str] | None = None) -> int:
    """Run the ``ga2o3d`` command-line interface."""
    parser = argparse.ArgumentParser(prog="ga2o3d", description="Ga2O3 defect dichroism helper CLI")
    sub = parser.add_subparsers(dest="cmd", required=True)
    sub.add_parser("cases")

    stage = sub.add_parser("stage")
    stage.add_argument("--input-dir", type=Path, default=Path("examples/crystal23/scf"))
    stage.add_argument("--out", type=Path, default=Path("runs"))
    stage.add_argument("--cores", type=int, default=128)

    prop = sub.add_parser("make-property")
    prop.add_argument("kind", choices=("bands", "doss", "anbd", "echg", "ech3", "ppan"))
    prop.add_argument("case_name")
    prop.add_argument("--out", type=Path, required=True)

    cphf = sub.add_parser("make-cphf")
    cphf.add_argument("input", type=Path)
    cphf.add_argument("--out", type=Path, required=True)
    cphf.add_argument("--dynamic", action="store_true")
    cphf.add_argument("--steps", type=int, default=1)
    cphf.add_argument("--start-nm", type=float, default=7500.0)
    cphf.add_argument("--stop-nm", type=float)
    cphf.add_argument("--damping", type=float)

    bands = sub.add_parser("plot-bands")
    bands.add_argument("--input", type=Path, required=True)
    bands.add_argument("--out", type=Path, required=True)
    bands.add_argument("--energy-unit", choices=("ev", "hartree"), default="ev")
    bands.add_argument("--emin", type=float, default=-1.0)
    bands.add_argument("--emax", type=float, default=5.0)

    dos = sub.add_parser("plot-dos")
    dos.add_argument("--input", type=Path, required=True)
    dos.add_argument("--out", type=Path, required=True)
    dos.add_argument("--energy-unit", choices=("ev", "hartree"), default="hartree")
    dos.add_argument("--emin", type=float, default=-1.0)
    dos.add_argument("--emax", type=float, default=5.0)

    parse = sub.add_parser("parse-cphf")
    parse.add_argument("--folder", type=Path, required=True)
    parse.add_argument("--out", type=Path, required=True)

    mulliken = sub.add_parser("plot-mulliken")
    mulliken.add_argument("--out", type=Path, required=True)

    args = parser.parse_args(argv)
    if args.cmd == "cases":
        for case in default_cases():
            print(f"{case.name:32s} {case.source_file:42s} {case.label}")
        return 0
    if args.cmd == "stage":
        print(stage_crystal_workflow(input_dir=args.input_dir, out_dir=args.out, cores=args.cores))
        return 0
    if args.cmd == "make-property":
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(property_input(args.kind, args.case_name), encoding="utf-8")
        print(args.out)
        return 0
    if args.cmd == "make-cphf":
        settings = CphfSettings.dynamic(start_nm=args.start_nm, stop_nm=args.stop_nm, steps=args.steps, damping=args.damping) if args.dynamic else CphfSettings.static()
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(insert_cphf_block(args.input.read_text(), settings), encoding="utf-8")
        print(args.out)
        return 0
    if args.cmd == "plot-bands":
        fig, _ = plot_bands(read_band(args.input, energy_unit=args.energy_unit), energy_window=(args.emin, args.emax))
        args.out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(args.out, dpi=600, bbox_inches="tight")
        return 0
    if args.cmd == "plot-dos":
        fig, _ = plot_dos(read_doss(args.input, energy_unit=args.energy_unit), energy_window=(args.emin, args.emax))
        args.out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(args.out, dpi=600, bbox_inches="tight")
        return 0
    if args.cmd == "parse-cphf":
        df = parse_cphf_folder(args.folder)
        args.out.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(args.out, index=False)
        print(args.out)
        return 0
    if args.cmd == "plot-mulliken":
        plot_mulliken_summary(out=args.out)
        print(args.out)
        return 0
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
