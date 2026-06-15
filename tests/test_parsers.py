from pathlib import Path

from ga2o3dichroism.parsers import parse_cphf_output, read_band, read_doss


def test_read_doss_numeric_table(tmp_path: Path) -> None:
    path = tmp_path / "sample.DOSS"
    path.write_text("# header\n0.0 1.0 2.0\n0.1 1.5 2.5\n")
    df = read_doss(path, energy_unit="ev")
    assert list(df.columns) == ["energy", "y1", "y2", "energy_ev"]
    assert df.shape == (2, 4)


def test_read_band_long_form(tmp_path: Path) -> None:
    path = tmp_path / "sample.BAND"
    path.write_text("0.0 -1.0 2.0\n1.0 -0.5 2.5\n")
    df = read_band(path)
    assert set(df.columns) == {"k", "band", "energy_ev"}
    assert len(df) == 4


def test_parse_cphf_filename_wavelength(tmp_path: Path) -> None:
    path = tmp_path / "case_7500nm.out"
    path.write_text("ABSORPTION 1.0 2.0 3.0\n")
    row = parse_cphf_output(path)
    assert row["wavelength_nm"] == 7500
    assert row["abs_x"] == 1.0
