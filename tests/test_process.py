import logging
import zipfile
from pathlib import Path

import pandas as pd
import pytest

from influx_si_data_manager.__main__ import args_parse, process

ROOT = Path(__file__).resolve().parents[1]
TEST_DATA = ROOT / "influx_si_data_manager" / "test_data"
MTF_DATA = ROOT / "data" / "influx_si_data"
MAPPING = TEST_DATA / "mapping.txt"


@pytest.fixture(autouse=True)
def reset_logger():
    yield
    logger = logging.getLogger("root")
    for handler in logger.handlers[:]:
        handler.close()
        logger.removeHandler(handler)


@pytest.fixture
def physiofit_file(tmp_path):
    """PhysioFit summary restricted to the experiment present in the IsoCor test data"""
    summary = pd.read_csv(TEST_DATA / "summary.csv")
    path = tmp_path / "summary.csv"
    summary[summary["experiments"] == "40mM_Rep1"].to_csv(path, index=False)
    return path


def run(tmp_path, monkeypatch, *extra_args, netw=True):
    monkeypatch.chdir(tmp_path)
    cli_args = ["--isocor", str(TEST_DATA / "isocor_results.tabular")]
    if netw:
        cli_args += ["--netw", str(MTF_DATA / "e_coli.netw")]
    process(args_parse().parse_args(cli_args + list(extra_args)))


def read_member(archive, member):
    with zipfile.ZipFile(archive) as zf, zf.open(member) as f:
        return pd.read_csv(f, sep="\t")


def mflux_values(tmp_path):
    mflux = read_member(tmp_path / "40mM_Rep1.zip", "40mM_Rep1.mflux")
    return dict(zip(mflux["Flux"], mflux["Value"]))


def test_mapping_is_applied_to_physiofit_data(tmp_path, monkeypatch, physiofit_file):
    run(tmp_path, monkeypatch, "--physiofit", str(physiofit_file), "--mapping", str(MAPPING))
    assert {"BM", "Aceupt", "Glucupt"} <= set(mflux_values(tmp_path))


def test_normalization_by_parameter_name(tmp_path, monkeypatch, physiofit_file):
    run(tmp_path, monkeypatch, "--physiofit", str(physiofit_file), "--mapping", str(MAPPING),
        "--normalization", "BM")
    values = mflux_values(tmp_path)
    assert values["BM"] == pytest.approx(1.0)
    assert values["Aceupt"] == pytest.approx(-0.7078295328322147 / 0.493829563404029)


def test_normalization_by_value(tmp_path, monkeypatch, physiofit_file):
    run(tmp_path, monkeypatch, "--physiofit", str(physiofit_file), "--mapping", str(MAPPING),
        "--normalization", "2")
    assert mflux_values(tmp_path)["Aceupt"] == pytest.approx(-0.7078295328322147 / 2)


def test_normalization_requires_physiofit_data(tmp_path, monkeypatch):
    with pytest.raises(ValueError, match="PhysioFit"):
        run(tmp_path, monkeypatch, "--normalization", "2")


def test_physiofit_data_is_optional(tmp_path, monkeypatch):
    run(tmp_path, monkeypatch)
    with zipfile.ZipFile(tmp_path / "AA_2.zip") as zf:
        members = zf.namelist()
    assert "AA_2.miso" in members
    assert "AA_2.netw" in members
    assert not any(member.endswith(".mflux") for member in members)


def test_log_defaults_to_working_directory(tmp_path, monkeypatch):
    run(tmp_path, monkeypatch)
    assert (tmp_path / "log.txt").is_file()


def test_galaxy_none_strings_mean_no_input(tmp_path, monkeypatch, physiofit_file):
    # Galaxy passes the string "None" for optional inputs left empty
    run(tmp_path, monkeypatch, "--physiofit", str(physiofit_file), "--mapping", str(MAPPING),
        "--linp", "None", "--normalization", "None", "--log", str(tmp_path / "run.log"))
    with zipfile.ZipFile(tmp_path / "40mM_Rep1.zip") as zf:
        members = zf.namelist()
    assert "40mM_Rep1.mflux" in members
    assert not any(member.endswith(".linp") for member in members)


def test_missing_network_file_is_an_error(tmp_path, monkeypatch):
    with pytest.raises(ValueError, match=r"\.netw"):
        run(tmp_path, monkeypatch, netw=False)
