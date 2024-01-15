"""
Module for converting isocor output into .miso mtf file for influx_si
"""

from pathlib import Path
import logging
import pandas as pd

_logger = logging.getLogger("root.isocor2mtf")


def isocor2mtf(
        isocor_res: str,
        sd: float = 0.02
):
    """
    Generate the list of dataframes in .miso format

    :param isocor_res: path to isocor results file
    :param sd: measurements standard deviation (defaults to 0.02)
    :return: None
    """

    # Get data
    _logger.debug(f"Selected default SD value: {sd}")
    _logger.info("Reading IsoCor data...")
    _logger.debug(f"IsoCor data path: {isocor_res}")
    data_path = Path(isocor_res)
    data = pd.read_csv(data_path, sep="\t")
    _logger.debug(f"IsoCor data:\n{data}")
    # Build dataframes with .miso structure
    useful_cols = ["sample", "metabolite", "isotopologue", "isotopologue_fraction"]
    _logger.info("Indexing data...")
    data = data[useful_cols]
    _logger.info("Generating fragments...")
    data = generate_fragments(data)
    _logger.debug(f"Data after fragment generation:\n{data}")
    _logger.info("Mapping new column names onto data...")
    column_name_map = {
        "isotopologue": "Isospecies",
        "isotopologue_fraction": "Value"
    }
    data = data.rename(
        mapper=column_name_map,
        axis=1
    )
    data["SD"] = sd
    data["Time"] = ""
    data["Id"] = ""
    data["Comment"] = ""
    data["Dataset"] = "MS-1"
    data = data[["sample", "Id", "Comment", "Specie", "Fragment", "Dataset", "Isospecies", "SD", "Time"]]
    _logger.info("Generating Isospecies...")
    data["Isospecies"] = "M" + data["Isospecies"].astype(str)
    miso_dfs = [
        (exp_name, data[data["sample"] == exp_name].copy())
        for exp_name in sorted(data["sample"].unique())
    ]
    _logger.debug("List of experiments and associated dataframes:")
    for (exp, df) in miso_dfs:
        _logger.debug(f"Experiment {exp:}\n{df}")

    return miso_dfs


def generate_fragments(data: pd.DataFrame):
    """
    Get the list of carbon numbers from MS and MS/MS experiments

    :param data: isocor res data
    :return: dataframe containing specie and fragment columns
    """

    # Handle MS/MS
    df_f = data[data.metabolite.str.contains("__f")].copy()
    df = data[~data.metabolite.str.contains("__f")].copy()
    if not df_f.empty:
        df_f[["Specie", "Fragment"]] = df_f.metabolite.str.split("__f", expand=True)
    tmps = []
    # Build Fragment and Specie columns
    for metabolite in df.metabolite.unique():
        tmp = df[df["metabolite"] == metabolite].copy()
        c = tmp["isotopologue"].max()
        s = ""
        for carbon in range(1, c + 1):
            s = s + f"{str(carbon)},"
        s = s[:-1]
        tmp["Fragment"] = s
        tmp["Specie"] = tmp.metabolite
        tmps.append(tmp)
    df = pd.concat(tmps)
    if not df_f.empty:
        data = pd.concat([df, df_f])
    else:
        data = df
    # Remove metabolite column (metabolite name is now in Specie)
    data = data.drop("metabolite", axis=1)
    return data

if __name__ == "__main__":

    mtf = isocor2mtf(r"C:\Users\legregam\PycharmProjects\influx_data_manager\test-data\isocor2mtf\Data_galaxy_test.tsv")
    print(mtf)
