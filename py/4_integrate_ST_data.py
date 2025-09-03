"""Single-cell transcriptomics data integration using scVI.

This module provides functionality to integrate single-cell RNA sequencing datasets
using the scVI (single-cell Variational Inference) model for batch correction and
dimensionality reduction.
"""

# Standard imports
import argparse
import datetime

# Third-party imports
import scvi

# Local imports
from utils import utils


__license__ = "MIT"
__copyright__ = "Copyright (c) 2025 Christian Kolland"


# Constants
TIMESTAMP_FORMAT: str = "%Y-%m-%d_%H-%M"
LOGFILE_NAME: str = "scVI_run.log"
DEFAULT_LAYER: str = "raw"
DEFAULT_N_LAYERS: int = 1
DEFAULT_HIDDEN_SIZE: int = 128
DEFAULT_LATENT_SIZE: int = 10
GENE_LIKELIHOOD: str = "nb"  # negative binomial


def main():
    """
    Integrates single-cell transcriptomics datasets using the `scVI` model.

    This function **parses command-line arguments** to configure the integration process,
    **loads an `AnnData` object**, **sets up the `scVI` model** for batch correction and dimensionality
    reduction, **trains the model**, and **saves the results** including the latent representation,
    normalized expression values.

    Command-line arguments:
        - `--adata` / `-d` (str): Path to the input .h5ad file containing the `AnnData` object.
        - --`layer` / `-l` (str, optional): Name of the layer holding raw counts. Defaults to "raw".
        - --`out_dir` / `-o` (str): Path to the output directory for results.
        - --`batch-effects` / `-b` (list of str): List of batch effect categories to correct for.
        - --`n-layers` / `-n` (int, optional): Number of hidden layers in the scVI model. Defaults to 1.
        - --`size-hidden` / `-s` (int, optional): Number of nodes per hidden layer. Defaults to 128.
        - --`size-latent` / `-t` (int, optional): Dimensionality of the latent space. Defaults to 10.
    """

    # --------------------
    #     Setup parser
    # --------------------

    parser = argparse.ArgumentParser(
        description="Integrates different single-cell datasets using the scVI model"
    )

    parser.add_argument(
        "--adata",
        "-d",
        type=str,
        required=True,
        help="Path to the input file (.h5ad) containing the AnnData object",
    )
    parser.add_argument(
        "--layer",
        "-l",
        type=str,
        default=DEFAULT_LAYER,
        help=f"Name of layer that holds raw counts (default: {DEFAULT_LAYER})",
    )
    parser.add_argument(
        "--out_dir",
        "-o",
        type=str,
        required=True,
        help="Path to the output directory for saving results",
    )
    parser.add_argument(
        "--batch-effects",
        "-b",
        nargs="+",
        type=str,
        required=True,
        help="List of batch effect categories to correct for (e.g., 'batch', 'donor')",
    )
    parser.add_argument(
        "--n-layers",
        "-n",
        type=int,
        default=DEFAULT_N_LAYERS,
        help=f"Number of hidden layers in the scVI model (default: {DEFAULT_N_LAYERS})",
    )
    parser.add_argument(
        "--size-hidden",
        "-s",
        type=int,
        default=DEFAULT_HIDDEN_SIZE,
        help=f"Number of nodes per hidden layer (default: {DEFAULT_HIDDEN_SIZE})",
    )
    parser.add_argument(
        "--size-latent",
        "-t",
        type=int,
        default=DEFAULT_LATENT_SIZE,
        help=f"Dimensionality of the latent space (default: {DEFAULT_LATENT_SIZE})",
    )

    args = parser.parse_args()

    # Validate and resolve file paths
    data_file = utils.assert_path(args.adata, assert_dir=False)
    out_dir = utils.assert_path(args.out_dir, assert_dir=True)

    # Create run directory
    now = datetime.datetime.now().strftime(TIMESTAMP_FORMAT)
    run_dir = out_dir / f"scVI_run_{now}"
    run_dir.mkdir(parents=True, exist_ok=True)  # Create run directory

    # ---------------------
    #     Setup logging
    # ---------------------

    log_file = run_dir / LOGFILE_NAME
    logger = utils.setup_logging(log_file)
    logger.info("✅ Setup complete.")
    logger.info("----")

    # ------------------------
    #     Data integration
    # ------------------------

    logger.info("Start integration process...")

    # Load the AnnData object
    adata = scvi.data.read_h5ad(data_file)
    logger.info(f"Successfully loaded dataset: '{data_file}'")
    logger.debug(adata)

    # Create backup of unprocessed data for safety
    backup_file = run_dir / f"unprocessed_{data_file.name}"
    adata.write(backup_file, compression="gzip")
    logger.debug(f"Saved unprocessed AnnData object as fallback: '{backup_file}'")

    # scVI needs raw counts as input!
    # Setup AnnData for scVI processing
    scvi.model.SCVI.setup_anndata(
        adata,
        layer=args.layer,
        categorical_covariate_keys=args.batch_effects,
    )
    logger.info(
        f"Configured scVI to correct for the following batch effects: {args.batch_effects}"
    )

    # Initialize the scVI model
    # Gene likelihood is set to `nb` (negative binomial)
    # -> Our CosMx data has no 0 expression genes
    model = scvi.model.SCVI(
        adata,
        n_layers=args.n_layers,
        n_hidden=args.size_hidden,
        n_latent=args.size_latent,
        gene_likelihood=GENE_LIKELIHOOD,
    )

    logger.info(
        f"Built scVI model with:\n"
        f"- {args.n_layers} hidden layers\n"
        f"- {args.size_hidden} nodes per hidden layer\n"
        f"- {args.size_latent} latent dimensions\n"
        f"- Gene likelihood: {GENE_LIKELIHOOD}"
    )

    # Train the model
    model.train()

    # Save model
    model_path = run_dir / f"scVI_trained.pt"
    model.save(model_path)
    logger.info(f"Saved trained model to: '{model_path}'")

    # -----------------------
    #     Extract results
    # -----------------------

    # Extract latent representation (low-dimensional embedding)
    # This captures the main biological variation while removing batch effects
    adata.obsm["X_scVI"] = model.get_latent_representation()
    logger.info(f"Saved latent representation to AnnData object.")

    # Extract normalized expression values
    # These are batch-corrected expression values suitable for downstream analysis
    adata.layers["scVI_normalized"] = model.get_normalized_expression()
    logger.info(f"Saved scVI normalized values to AnnData object.")

    # Save the updated AnnData object with all results
    adata_file = run_dir / data_file.name
    adata.write(adata_file, compression="gzip")
    logger.info(f"💾 Saved AnnData object to: '{adata_file}'")

    logger.info("🎉 Finished!")


if __name__ == "__main__":
    main()
