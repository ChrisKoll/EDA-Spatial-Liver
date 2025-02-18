# Standard imports
import argparse
import datetime
import logging
import os
import sys

# Third-party imports
import scanpy as sc
import scvi
import torch

# Set options
sc.settings.verbosity = 0
scvi.settings.dl_num_workers = 7
scvi.settings.verbosity = 10
torch.set_float32_matmul_precision("high")


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


def integrate_data():
    """Integrates single-cell RNA sequencing data using scVI and generates visualizations.

    This function integrates single-cell counts from different datasets stored in one `AnnData` object. The scVI model is trained and the results are extracted and stored in the `AnnData` object. Lastly the results are further visualized using UMAP embeddings and Leiden clustering.

    Raises:
        FileNotFoundError: If the specified input data file does not exist.
        ValueError: If invalid arguments are provided.
    """

    # |================|
    # | Setup Function |
    # |================|

    # Setup parser
    parser = argparse.ArgumentParser(description="Integrate data using scVI")

    # Add flags
    parser.add_argument(
        "--data", "-d", type=str, required=True, help="Path to the data file"
    )
    parser.add_argument(
        "--out_dir", "-o", type=str, required=True, help="Path to the output directory"
    )
    parser.add_argument(
        "--categories",
        "-c",
        nargs="+",
        type=str,
        required=True,
        help="List of categories to correct for",
    )

    # Parse arguments
    args = parser.parse_args()

    # Setup logging
    # -> Formatter used to add file logging with givne output directory
    logger, formatter = setup_logging()

    logger.debug("Logger setup complete.")

    data = validate_path(args.data, is_dir=False, logger=logger)
    out_dir = validate_path(args.out_dir, is_dir=True, logger=logger)

    # Create directory for results
    now = datetime.datetime.now().strftime("%d-%m-%Y_%H-%M")
    results_dir = os.path.join(
        out_dir, f"data_integration_{now}"
    )  # Save to add content later
    os.makedirs(results_dir, exist_ok=True)

    # |===============|
    # | Extend Logger |
    # |===============|

    # Logs to given log file
    file_handler = logging.FileHandler(
        os.path.join(results_dir, f"data_integration_{now}.log")
    )
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    # |==============|
    # | Prepare Data |
    # |==============|

    logger.info("> Start data integration")
    logger.info(f"Read data file: {data}")

    # Load data
    adata = scvi.data.read_h5ad(data)

    logger.info("Data was successfully loaded")
    logger.debug(adata)

    logger.info(f"> Correct for following categories: {args.categories}")

    # Retrieve raw counts
    bdata = adata.raw.to_adata()

    # Setup the data for processing
    # Model should correct for list of given categories
    scvi.model.SCVI.setup_anndata(bdata, categorical_covariate_keys=args.categories)

    # |=============|
    # | Train Model |
    # |=============|

    model_layers = 2
    latent_dim = 30

    logger.info(
        f"> Build model with: {model_layers} layers, {latent_dim} latent dimensions"
    )

    # Standard values from tutorial
    # https://docs.scvi-tools.org/en/stable/tutorials/notebooks/scrna/harmonization.html
    model = scvi.model.SCVI(
        bdata, n_layers=model_layers, n_latent=latent_dim, gene_likelihood="nb"
    )

    model.train()

    # Save model
    model_dir = os.path.join(results_dir, f"scvi_model_{now}")
    model.save(model_dir)

    # |=================|
    # | Extract Results |
    # |=================|

    # Get latent representation
    latent = model.get_latent_representation()
    adata.obsm["latent_scVI"] = latent

    logger.info(f"Extracted latent representation of shape {latent.shape}")

    # Save normalized values
    # -> Library size scales the values to common value
    adata.layers["scVI_normalized"] = model.get_normalized_expression(library_size=10e4)

    logger.info(f"Saved scVI normalized values")

    # |===================|
    # | Visualize Results |
    # |===================|

    # Saves all figures to given directory
    sc.settings.figdir = results_dir

    # Calculate UMAP for non-corrected data
    # Visualize afterwards
    sc.tl.pca(bdata)
    sc.pp.neighbors(bdata, n_pcs=30, n_neighbors=20)
    sc.tl.umap(bdata, min_dist=0.3)

    sc.pl.umap(
        bdata,
        color=["slide"],
        frameon=False,
        save="_slide_not_corrected.png",
        show=False,
    )
    sc.pl.umap(
        bdata,
        color=["condition"],
        frameon=False,
        save="_condition_not_corrected.png",
        show=False,
    )
    sc.pl.umap(
        bdata,
        color=["fov"],
        frameon=False,
        save="_fov_not_corrected.png",
        show=False,
    )

    # Calculate UMAP for corrected data
    # Visualize afterwards
    sc.pp.neighbors(adata, use_rep="latent_scVI")
    sc.tl.umap(adata, min_dist=0.3)

    sc.pl.umap(
        adata,
        color=["slide"],
        frameon=False,
        save="_slide_corrected.png",
        show=False,
    )
    sc.pl.umap(
        adata,
        color=["condition"],
        frameon=False,
        save="_condition_corrected.png",
        show=False,
    )
    sc.pl.umap(
        adata,
        color=["fov"],
        frameon=False,
        save="_fov_corrected.png",
        show=False,
    )

    # Cluster the latent space
    sc.tl.leiden(
        adata,
        key_added="leiden_scVI",
        resolution=0.5,
        flavor="igraph",
        n_iterations=2,
        directed=False,
    )

    sc.pl.umap(
        adata,
        color=["leiden_scVI"],
        frameon=False,
        save="_leiden_corrected.png",
        show=False,
    )

    adata.write(os.path.join(results_dir, os.path.basename(data)))

    logger.info(f"💾 Saved AnnData object to: {os.path.basename(data)}")
    logger.info(f"> Finished!")


if __name__ == "__main__":
    integrate_data()
