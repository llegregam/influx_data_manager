import argparse
import logging
import zipfile
from pathlib import Path

import pandas as pd

from influx_si_data_manager.utils.isocor2mtf import isocor2mtf
from influx_si_data_manager.utils.physiofit2mtf import physiofit2mtf, normalize_data
from influx_si_data_manager.utils.map_data import map_data


def _init_logger(log_path, debug=False):
    """
    Initialize the root logger
    :return:
    """

    logger = logging.getLogger("root")
    logger.setLevel(logging.DEBUG)
    # stream_handler = logging.StreamHandler(sys.stdout)
    # stream_handler.setLevel(logging.DEBUG)
    file_handler = logging.FileHandler(log_path, mode="w")
    if debug:
        file_handler.setLevel(logging.DEBUG)
    else:
        file_handler.setLevel(logging.INFO)
    formatter = logging.Formatter(
        "%(levelname)s:%(name)s: %(message)s"
    )
    # stream_handler.setFormatter(formatter)
    file_handler.setFormatter(formatter)
    # logger.addHandler(stream_handler)
    logger.addHandler(file_handler)


def args_parse():
    """
    Parse arguments from Command-Line Interface
    :return: Argument Parser containing args
    """

    parser = argparse.ArgumentParser(
        "Influx_si data connector: enabling tool interoperability on Workflow4Metabolomics, see here: "
        "workflow4metabolomics.usegalaxy.org"
    )

    # parser.add_argument(
    #    "-v", "--version", action="version",
    # )
    parser.add_argument(
        "-m", "--mapping", type=str,
        help="Path to the mapping file"
    )

    # Get paths to data
    influx_files = parser.add_argument_group(
        "Influx input files",
        "List of all the input files that influx_si can accept. Some are variable and change for each experiment "
        "(.mflux or physiofit output, .miso or isocor output, .mmet, etc) and some are static for the same batch of "
        "experiments (.netw, .linp, etc). Check out the influx user documentation for more information."
    )

    influx_files.add_argument(
        "-p", "--physiofit", type=str,
        help="Path to physiofit summary output file"
    )
    influx_files.add_argument(
        "-i", "--isocor", type=str,
        help="Path to isocor results output file"
    )
    influx_files.add_argument(
        "-li", "--linp", type=str,
        help="Path to .linp file containing the label input information"
    )
    influx_files.add_argument(
        "-ne", "--netw", type=str,
        help="Path to .netw file containing the list of biochemical reactions & label transitions"
    )
    influx_files.add_argument(
        "-mm", "--mmet", type=str,
        help="Path to .mmet file containing the stationary specie information"
    )
    influx_files.add_argument(
        "-cn", "--cnstr", type=str,
        help="Path to .cnstr file containing constraints on fluxes and specie conentrations"
    )
    influx_files.add_argument(
        "-tv", "--tvar", type=str,
        help="Path to .tvar file containing the types of variables"
    )
    influx_files.add_argument(
        "-op", "--opt", type=str,
        help="Path to .opt file containing extra options to pass to influx_si"
    )
    parser.add_argument(
        "-l", "--log", type=str,
        help="Output path for log"
    )
    parser.add_argument(
        '-n', '--normalization', action='store', type=str, default='None',
        help='Normalize extracellular fluxes by a value in the input extracellular fluxes ' \
        'dataset (e.g. glc uptake) or by a specific user-defined value'
    )
    parser.add_argument(
        '-v', '--verbose', action='store_true',
        help='Activate debug mode'
    )

    return parser


def _is_given(value):
    """
    Check if an optional input was given. Galaxy passes the string 'None' for
    optional inputs left empty, the CLI passes None.
    """
    return value not in (None, "None", "")


def _parse_normalization(value):
    """
    Normalization is either a number or the name of a PhysioFit parameter
    """
    try:
        return float(value)
    except ValueError:
        return value


def process(args):

    # initialize root
    log_path = args.log if _is_given(args.log) else "./log.txt"
    _init_logger(str(Path(log_path)), args.verbose)

    # get logger
    _logger = logging.getLogger("root")

    _logger.debug("Run arguments:")
    for key, val in vars(args).items():
        _logger.debug(f"{key} : {val}")

    if not _is_given(args.netw):
        msg = 'Network file containing reaction and carbon transitions (.netw) is mandatory.'
        _logger.error(msg)
        raise ValueError(msg)

    _logger.info("Generating mflux and miso dataframes...")

    isocor_data = pd.read_csv(args.isocor, sep="\t")
    physiofit_data = None
    if _is_given(args.physiofit):
        _logger.info("Reading PhysioFit data...")
        physiofit_data = pd.read_csv(args.physiofit, sep=",")
    else:
        _logger.info("No PhysioFit data given, archives will not contain .mflux files")

    if _is_given(args.mapping):
        _logger.info(f"Mapping file detected: {args.mapping}")
        isocor_data = map_data(
            mapping_file=args.mapping,
            data=isocor_data,
            from_tool="isocor"
        )
        if physiofit_data is not None:
            physiofit_data = map_data(
                mapping_file=args.mapping,
                data=physiofit_data,
                from_tool="physiofit"
            )

    if _is_given(args.normalization):
        if physiofit_data is None:
            msg = "Normalization of extracellular fluxes requires PhysioFit data"
            _logger.error(msg)
            raise ValueError(msg)
        norm_value = _parse_normalization(args.normalization)
        _logger.info(f"Normalization of extracellular fluxes by {norm_value}")
        physiofit_data = normalize_data(
            physiofit_data=physiofit_data,
            norm_value=norm_value
        )
        _logger.info(f"Normalized PhysioFit data:\n{physiofit_data}")

    miso_dfs = isocor2mtf(isocor_res=isocor_data)
    miso_names = [exp[0] for exp in miso_dfs]

    mflux_by_experiment = {}
    if physiofit_data is not None:
        mflux_dfs = physiofit2mtf(data=physiofit_data)
        mflux_names = [exp[0] for exp in mflux_dfs]
        if mflux_names != miso_names:
            msg = (f"Sample names in miso files and mflux files are not the same:\nmflux names: {mflux_names}"
                   f"\nmiso names: {miso_names}")
            _logger.error(msg)
            raise ValueError(msg)
        mflux_by_experiment = dict(mflux_dfs)

    _logger.info(f"Experiment Names:\n{miso_names}")

    # List of files that should be static (non variable)
    non_var_files = [
        "linp",
        "netw",
        "mmet",
        "cnstr",
        "tvar",
        "opt"
    ]

    _logger.info("Building archives...")
    # Build archive & export for discovery in Galaxy workflow
    for experiment, miso in miso_dfs:
        with zipfile.ZipFile(f"{experiment}.zip", "w", compression=zipfile.ZIP_DEFLATED) as output_zip:
            _logger.info(f"Building archive for experiment {experiment}")

            # Handle isocor output (i.e corrected labelling) file:
            with output_zip.open(f"{experiment}.miso", "w") as miso_file:
                _logger.info(f'Adding {experiment}.miso')
                _logger.info(f"Data:\n{miso}")
                miso.to_csv(miso_file, index=False, sep="\t")

            # Handle physiofit output (i.e extracellular fluxes) file:
            if experiment in mflux_by_experiment:
                mflux = mflux_by_experiment[experiment]
                with output_zip.open(f"{experiment}.mflux", "w") as mflux_file:
                    _logger.info(f'Adding {experiment}.mflux')
                    _logger.info(f"Data:\n{mflux}")
                    mflux.to_csv(mflux_file, index=False, sep="\t")

            # Handle the other mtf files
            for nvf in non_var_files:
                nvf_file_path = vars(args)[nvf]
                _logger.debug(f"nvf file path: {nvf_file_path}")
                if not _is_given(nvf_file_path):
                    _logger.info(f'No {nvf} file given')
                    continue
                if nvf == "linp":
                    df = pd.read_csv(nvf_file_path, sep="\t", comment="#", dtype={"Isotopomer": str})
                else:
                    df = pd.read_csv(nvf_file_path, sep="\t", comment="#")
                with output_zip.open(f"{experiment}.{nvf}", "w") as nvf_file:
                    _logger.info(f'Adding {experiment}.{nvf}')
                    _logger.info(f"Data:\n{df}")
                    df.to_csv(nvf_file, index=False, sep="\t")


def main():
    parser = args_parse()
    args = parser.parse_args()
    process(args)


if __name__ == "__main__":
    main()
