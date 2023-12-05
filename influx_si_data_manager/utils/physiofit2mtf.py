"""
Module for converting PhysioFit output into .mflux mtf file for influx_si
"""
import logging
from pathlib import Path

import pandas as pd

_logger = logging.getLogger("root.physiofit2mtf")

def physiofit2mtf(
        physiofit_res: str,
):
    """
    Generate dataframes with the .mflux structure (1 dataframe per experiment/mflux file)

    :param physiofit_res: path to physiofit summary file
    :return: list containing the dfs
    """

    # Get data
    _logger.info("Reading PhysioFit data...")
    data_path = Path(physiofit_res)
    data = pd.read_csv(data_path, sep=",")
    _logger.debug(f"PhysioFit data before indexing:\n{data}")
    data = data.loc[(~data["parameter name"].str.contains("_M0")) & (~data["parameter name"].str.contains("_0"))]
    _logger.debug(f"PhysioFit data after indexing:\n{data}")
    # build .mflux
    _logger.info("Building .mflux file...")
    mflux_file = pd.DataFrame(columns=["Id", "Comment", "Flux", "Value", "SD"])
    mflux_file["Id"] = data["experiments"]
    mflux_file["Flux"] = data["parameter name"]
    mflux_file["Value"] = data["optimal"]
    mflux_file["SD"] = data["sd"]
    mflux_dfs = [
        (exp_name, mflux_file[mflux_file["Id"] == exp_name].copy())
        for exp_name in sorted(mflux_file["Id"].unique())
    ]
    _logger.debug("List of PhysioFit experiments and associated dataframes:")
    for (exp, df) in mflux_dfs:
        _logger.debug(f"Experiment {exp:}\n{df}")
    return mflux_dfs

