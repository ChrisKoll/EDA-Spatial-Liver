# Standard imports
import argparse
import datetime
import gzip
import logging
import os
import sys

# Third-party imports
import scanpy as sc
import pandas as pd
from scipy import io
from scipy.sparse import csr_matrix


def setup_logging() -> tuple[logging.Logger, logging.Formatter]:
    """Sets up logging for the application with a console handler and a specific log format.
    The logger is configured to output logs to the console (stdout) with a timestamp,
    log level, and message.

    This function returns the configured logger and formatter so they can be used
    in other parts of the application.

    :return: A tuple containing the configured logger instance and the formatter used.
    :rtype: tuple[logging.Logger, logging.Formatter]
    """

    # Create a logger
    logger = logging.getLogger("data_integration")
    logger.setLevel(logging.DEBUG)  # <- Important: Set the log level

    # Prevent duplicate handlers
    if logger.hasHandlers():
        logger.handlers.clear()

    # Create formatter and add it to handlers
    formatter = logging.Formatter(
        "%(asctime)s |:| LEVEL: %(levelname)s |:| %(message)s",
        datefmt="%d-%m-%Y %H:%M:%S",
    )

    # Create console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)

    # Add handlers to the logger
    logger.addHandler(console_handler)

    return logger, formatter


def validate_path(
    path: str, is_dir: bool = False, logger: logging.Logger = None
) -> str:
    """Validates whether the given path is a valid directory or file.

    If `is_dir` is set to `True`, the function checks if the path corresponds to an existing directory. If `is_dir` is `False`, it checks if the path corresponds to an existing file.

    If the path is invalid, a `FileNotFoundError` is raised, and an error message is logged if a logger is provided. If the path is valid, a success message is printed and logged (if a logger is provided).

    :param path: The path to the file or directory to be validated.
    :type path: str
    :param is_dir: Boolean flag indicating whether to validate a directory (`True`) or a file (`False`), defaults to `False`.
    :type is_dir: bool
    :param logger: Optional logger instance for logging validation messages. If not provided, no logging occurs.
    :type logger: Logger, optional

    :raises FileNotFoundError: If the path does not exist.

    :return: The validated path if the validation is successful.
    :rtype: str
    """

    if is_dir:  # Given path is directory
        if not os.path.isdir(path):  # If path is not a valid directory
            err_msg = f"Directory not found: {path}"

            if logger:  # If logger is provided, log error in log file
                logger.critical(err_msg)

            raise FileNotFoundError(err_msg)
        else:
            succ_msg = f"{'📁' if is_dir else '💾'} Valid {'directory' if is_dir else 'file'}: {path}"

            if logger:  # If logger is provided, log success in log file
                logger.info(succ_msg)

            return path
    else:  # Given path is file
        if not os.path.isfile(path):  # If path is not a valid file
            err_msg = f"File not found: {path}"

            if logger:  # If logger is provided, log error in log file
                logger.critical(err_msg)

            raise FileNotFoundError(err_msg)
        else:
            succ_msg = f"{'📁' if is_dir else '💾'} Valid {'directory' if is_dir else 'file'}: {path}"

            if logger:  # If logger is provided, log success in log file
                logger.info(succ_msg)

            return path


def convert_training_data():

    # |================|
    # | Setup Function |
    # |================|

    # Setup parser
    parser = argparse.ArgumentParser(
        description="Convert CellRanger mouse liver training data to AnnData"
    )

    # Add flags
    parser.add_argument(
        "--data", "-d", type=str, required=True, help="Path to the data file"
    )
    parser.add_argument(
        "--out_dir", "-o", type=str, required=True, help="Path to the output directory"
    )

    # Parse arguments
    args = parser.parse_args()

    # Setup logging
    # -> Formatter used to add file logging with givne output directory
    logger, formatter = setup_logging()

    logger.debug("Logger setup complete.")

    data_dir = validate_path(args.data, is_dir=True, logger=logger)
    out_dir = validate_path(args.out_dir, is_dir=True, logger=logger)

    # Create directory for results
    now = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M")
    results_dir = os.path.join(
        out_dir, f"training_data_integration_{now}"
    )  # Save to add content later
    os.makedirs(results_dir, exist_ok=True)

    # |===============|
    # | Extend Logger |
    # |===============|

    # Logs to given log file
    file_handler = logging.FileHandler(
        os.path.join(results_dir, f"training_data_integration_{now}.log")
    )
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    # |==============|
    # | Prepare Data |
    # |==============|

    logger.info("> Start data conversion")
    logger.info(f"Read data folder: {data_dir}")

    # CellRanger outputs in directory
    # Stored as `.mtx` format
    matrix = io.mmread(os.path.join(data_dir, "matrix.mtx.gz")).tocsc()

    # Load barcodes
    with gzip.open(os.path.join(data_dir, "barcodes.tsv.gz"), "rt") as f:
        barcodes = [line.strip() for line in f]

    # Load gene names
    with gzip.open(os.path.join(data_dir, "features.tsv.gz"), "rt") as f:
        genes = [line.strip() for line in f]

    # Construct AnnData manually
    adata = sc.AnnData(X=matrix.T)
    adata.var_names = genes
    adata.obs_names = barcodes

    logger.info("Data was successfully loaded as AnnData")
    logger.debug(adata)

    adata.write(
        os.path.join(results_dir, f"{os.path.basename(data_dir.rstrip('/'))}.h5ad")
    )

    logger.info(f"💾 Saved AnnData object to: {os.path.basename(data_dir.rstrip('/'))}")
    logger.info(f"> Finished!")

    # Log some quality metrics
    logger.info("==Quality Metrics==")

    logger.info(f"Shape of the AnnData object: {adata.shape}")
    logger.info(adata.obs.head())
    logger.info(adata.var.head())

    sc.pp.calculate_qc_metrics(adata, inplace=True)

    sc.settings.figdir = results_dir
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts"],
        jitter=0.4,
        multi_panel=True,
        save="mouse_liver_stst.svg",
    )


if __name__ == "__main__":
    convert_training_data()
