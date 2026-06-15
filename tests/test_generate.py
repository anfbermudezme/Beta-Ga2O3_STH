from pathlib import Path

from ga2o3dichroism.generate import CphfJob, InputCase, generate_crystal_inputs, load_workflow


def test_load_publication_workflow_contains_cphf_jobs() -> None:
    cases = load_workflow("examples/workflows/publication.yml")
    tetra = next(case for case in cases if case.name == "GaO_TETRA_VAC")
    assert tetra.source_file == "Ga2O3_VAC_TETRA_2x2x2.d12"
    assert "bands" in tetra.properties
    assert len(tetra.cphf_jobs) == 2


def test_generate_crystal_inputs_writes_d12_d3_scripts_and_manifest(tmp_path: Path) -> None:
    case = InputCase(
        name="GaO_TETRA_VAC",
        source_file="Ga2O3_VAC_TETRA_2x2x2.d12",
        label="test case",
        properties=("bands", "doss"),
        cphf_jobs=(CphfJob.dynamic("GaO_TETRA_VAC", start_nm=7500, stop_nm=12500, steps=10, damping=0.002),),
    )
    result = generate_crystal_inputs(
        input_dir="examples/crystal23/scf",
        out_dir=tmp_path / "runs",
        cores=128,
        cases=(case,),
    )
    case_dir = tmp_path / "runs" / "crystal23" / "GaO_TETRA_VAC"
    assert (case_dir / "GaO_TETRA_VAC.d12").exists()
    assert (case_dir / "CPHF_GaO_TETRA_VAC_7500NM_12500NM.d12").exists()
    assert (case_dir / "BANDS_GaO_TETRA_VAC.d3").exists()
    assert (case_dir / "DOSS_GaO_TETRA_VAC.d3").exists()
    assert (case_dir / "run_case.sh").exists()
    assert result.manifest.exists()
    manifest = result.manifest.read_text()
    assert "runPcry23 128 GaO_TETRA_VAC" in manifest
    assert "runPprop23 128 BANDS_GaO_TETRA_VAC GaO_TETRA_VAC" in manifest
