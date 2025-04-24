# Standard imports
import argparse
import datetime
import logging
import os
import sys

# Third-party imports
import scanpy as sc


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

    # Create directory for results
    now = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M")
    results_dir = os.path.join(
        args.out_dir, f"training_data_integration_{now}"
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
    logger.info(f"Read data file: {args.data}")

    # CellRanger outputs in directory
    # Stored as `.mtx` format
    adata = sc.read_10x_mtx(args.data, make_unique=True)

    logger.info("Data was successfully loaded")
    logger.debug(adata)

    adata.write(os.path.join(results_dir, os.path.basename(args.data)))

    logger.info(f"💾 Saved AnnData object to: {os.path.basename(args.data)}")
    logger.info(f"> Finished!")

    # Log some quality metrics
    logger.info("==Quality Metrics==")

    logger.info(f"Shape of the AnnData object: {adata.shape}")
    logger.info(adata.obs.head())
    logger.info(adata.var.head())

    sc.pp.calculate_qc_metrics(adata, inplace=True)
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts"],
        jitter=0.4,
        multi_panel=True,
        save="mla_stst.png",
    )


if __name__ == "__main__":
    convert_training_data()
