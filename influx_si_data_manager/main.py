import argparse
import logging
import zipfile
from pathlib import Path

import pandas as pd

from influx_si_data_manager.utils.isocor2mtf import isocor2mtf
from influx_si_data_manager.utils.physiofit2mtf import physiofit2mtf


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
        "Influx_si data manager: handling data interoperability on Workflow4Metabolomics, see here: "
        "workflow4metabolomics.usegalaxy.org"
    )

    influx_files = parser.add_argument_group(
        "Influx input files",
        "List of all the input files that influx_si can accept. Some are variable and change for each experiment "
        "(.mflux or physiofit output, .miso or isocor output, .mmet, etc) and some are static for the same batch of "
        "experiments (.netw, .linp, etc). Check out the influx user documentation for more information."
    )

    # Get paths to data
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

    # # Give output collection paths
    # parser.add_argument(
    #     "-z", "--zip", action='store_true',
    #     help="Output path for zip containing all the files for influx launch"
    # )
    parser.add_argument(
        "-l", "--log", type=str,
        help="Output path for log"
    )
    parser.add_argument(
        '-v', '--verbose', action='store_true',
        help='Activate debug mode'
    )

    return parser


def process(args):
    # initialize root
    if hasattr(args, "log"):
        _init_logger(str(Path(args.log)), args.verbose)
    else:
        _init_logger("./log.txt", args.verbose)

    # get logger
    _logger = logging.getLogger("root")

    _logger.debug("Run arguments:")
    for key, val in vars(args).items():
        _logger.debug(f"{key} : {val}")

    _logger.info("Generating mflux and miso dataframes...")
    mflux_dfs = physiofit2mtf(
        physiofit_res=args.physiofit
    )
    miso_dfs = isocor2mtf(
        isocor_res=args.isocor
    )

    # Check experiment names
    mflux_names = [exp[0] for exp in mflux_dfs]
    miso_names = [exp[0] for exp in miso_dfs]

    if mflux_names != miso_names:
        msg = (f"Sample names in miso files and mflux files are not the same:\nmflux names: {mflux_names}"
               f"\nmiso names: {miso_names}")
        _logger.error(msg)
        raise ValueError(msg)

    _logger.info(f"Experiment Names:\n{mflux_names}")

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
    for mflux, miso in zip(mflux_dfs, miso_dfs):
        with zipfile.ZipFile(f"{mflux[0]}.zip", "w", compression=zipfile.ZIP_DEFLATED) as output_zip:
            _logger.info(f"Building archive for experiment {mflux[0]}")

            # Handle isocor output (i.e corrected labelling) file:
            with output_zip.open(f"{miso[0]}.miso", "w") as miso_file:
                _logger.info(f'Adding {miso[0]}.miso')
                _logger.info(f"Data:\n{miso[1]}")
                miso[1].to_csv(miso_file, index=False, sep="\t")

            # Handle physiofit output (i.e extracellular fluxes) file:
            if mflux[0] != 'None':
                with output_zip.open(f"{mflux[0]}.mflux", "w") as mflux_file:
                    _logger.info(f'Adding {mflux[0]}.mflux')
                    _logger.info(f"Data:\n{mflux[1]}")
                    mflux[1].to_csv(mflux_file, index=False, sep="\t")

            # Handle the other mtf files
            for nvf in non_var_files:
                nvf_file_path = vars(args)[nvf]
                nvf_file_type = type(nvf_file_path)
                if nvf == 'netw' and nvf_file_path == 'None':
                    msg = 'Network file containing reaction and carbon transitions (.netw) is mandatory.'
                    _logger.error(msg)
                    raise ValueError(msg)
                _logger.debug(f"nvf file path: {nvf_file_path}")
                _logger.debug(f'nvf file type: {nvf_file_type}')
                if nvf_file_path != 'None':
                    _logger.info(f'Adding {nvf_file_path}.{nvf}')
                    df = pd.read_csv(nvf_file_path, sep="\t", comment="#")
                    with output_zip.open(f"{mflux[0]}.{nvf}", "w") as nvf_file:
                        _logger.info(f'Adding {mflux[0]}.{nvf}')
                        _logger.info(f"Data:\n{df}")
                        df.to_csv(nvf_file, index=False, sep="\t")
                else:
                    _logger.info(f'No {nvf} file given')

    # mydir = Path(".").glob('**/*')
    # files = [x for x in mydir]
    # _logger.info(files)


def main():
    parser = args_parse()
    args = parser.parse_args()
    process(args)


if __name__ == "__main__":
    main()
