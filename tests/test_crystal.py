from pathlib import Path

from ga2o3dichroism.crystal import CphfSettings, cphf_block, insert_cphf_block, property_input


def test_band_input_has_expected_keywords() -> None:
    text = property_input("bands", "GaO_TETRA_VAC")
    assert text.startswith("BAND")
    assert "BANDS_GaO_TETRA_VAC" in text
    assert text.rstrip().endswith("END")


def test_dynamic_cphf_block_contains_wavelength_range() -> None:
    settings = CphfSettings.dynamic(start_nm=7500, stop_nm=12500, steps=10, damping=0.002)
    text = cphf_block(settings)
    assert "DYNAMIC" in text
    assert "7500" in text
    assert "12500" in text
    assert "DAMPING" in text


def test_insert_cphf_block_before_geometry_end() -> None:
    source = Path("examples/crystal23/scf/Ga2O3_VAC_TETRA_2x2x2.d12")
    text = insert_cphf_block(source.read_text(), CphfSettings.static())
    lines = [line.strip().upper() for line in text.splitlines()]
    assert lines.index("CPHF") < lines.index("END")
