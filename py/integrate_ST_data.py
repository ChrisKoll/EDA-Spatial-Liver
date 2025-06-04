# Standard imports
import argparse
import datetime
import os
from pathlib import Path

# Third-party imports
from anndata import AnnData
import scanpy as sc
import scvi
import torch

# Local imports
import utils


def visualize_umap(adata: AnnData, layer: str, batch_effects: list[str]) -> None:
    """
    Visualizes UMAP embeddings for a given AnnData object using a specified data layer and highlights batch effects.

    :param adata: The annotated data matrix (AnnData object) containing the data to visualize.
    :type adata: AnnData

    :param layer: The key of the layer in `adata.layers` to use for UMAP calculation.
    :type layer: str

    :param batch_effect: A list of column names in `adata.obs` representing batch effect variables to visualize on the UMAP.
    :type batch_effect: list[str]

    :returns: The function saves UMAP plots to disk and does not return any value.
    :rtype: None
    """

    # Calculate UMAP for layer of choice
    sc.pp.neighbors(adata, use_rep=layer)
    sc.tl.umap(adata)

    # Plot UMAP for each batch effect to see if it is corrected
    for effect in batch_effects:
        sc.pl.umap(
            adata,
            color=effect,
            frameon=False,
            show=False,
            save=f"_{effect}_corrected.png",
        )


def main():
    """
    Integrates single-cell transcriptomics datasets using the scVI model.

    This function parses command-line arguments to configure the integration process,
    loads an AnnData object, sets up the scVI model for batch correction and dimensionality
    reduction, trains the model, and saves the results including the latent representation,
    normalized expression values, and visualizations.

    Command-line arguments:
        - --adata / -d (str): Path to the input .h5ad file containing the AnnData object.
        - --layer / -l (str, optional): Name of the layer holding raw counts. Defaults to "raw".
        - --out_dir / -o (str): Path to the output directory for results.
        - --batch-effects / -b (list of str): List of batch effect categories to correct for.
        - --n-layers / -n (int, optional): Number of hidden layers in the scVI model. Defaults to 1.
        - --size-hidden / -s (int, optional): Number of nodes per hidden layer. Defaults to 128.
        - --size-latent / -t (int, optional): Dimensionality of the latent space. Defaults to 10.
    """

    # --------------------
    #     Setup parser
    # --------------------

    parser = argparse.ArgumentParser(
        description="Integrates different datasets using the scVI model"
    )

    parser.add_argument(
        "--adata", "-d", type=str, required=True, help="Path to the .h5ad file"
    )
    parser.add_argument(
        "--layer",
        "-l",
        type=str,
        default="raw",
        help="Name of layer that holds raw counts",
    )
    parser.add_argument(
        "--out_dir", "-o", type=str, required=True, help="Path to the output directory"
    )
    parser.add_argument(
        "--batch-effects",
        "-b",
        nargs="+",
        type=str,
        required=True,
        help="List of batch effects to correct for",
    )
    parser.add_argument(
        "--n-layers",
        "-n",
        type=int,
        default=1,
        help="Number of hidden layers in the model",
    )
    parser.add_argument(
        "--size-hidden",
        "-s",
        type=int,
        default=128,
        help="Number of nodes per hidden layer",
    )
    parser.add_argument(
        "--size-latent",
        "-t",
        type=int,
        default=10,
        help="Dimensionality of the latent space",
    )

    args = parser.parse_args()

    # Verify paths
    data_file = utils.verify_path(args.data, is_dir=False)
    out_dir = utils.verify_path(args.out_dir, is_dir=True)

    # Create run directory
    now = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M")
    run_dir = os.path.join(out_dir, f"run_{now}")
    os.makedirs(run_dir, exist_ok=True)  # Create run directory

    # ---------------------
    #     Setup logging
    # ---------------------

    logger = utils.setup_logging(out_dir)
    logger.info("Setup complete.")
    logger.info("Start integration process.")

    # ------------------------
    #     Data integration
    # ------------------------

    # Load data
    adata = scvi.data.read_h5ad(data_file)
    logger.info(f"Dataset loaded from `{data_file}`")
    logger.debug(adata)

    # scVI needs raw counts as input!
    # Setup data for processing
    scvi.model.SCVI.setup_anndata(
        adata,
        layer=args.layer,
        categorical_covariate_keys=args.categories,
    )
    logger.info(f"Correct for the following categories: {args.categories}")

    # Setup model
    # Gene likelihood is set to `nb` (negative binomial)
    # -> CosMx data has no 0 expression genes
    model = scvi.model.SCVI(
        adata,
        n_layers=args.n_layers,
        n_hidden=args.size_hidden,
        n_latent=args.size_latent,
        gene_likelihood="nb",
    )
    logger.info(
        f"Build model:\n- {args.n_layers} layers\n- {args.size_hidden} nodes per hidden layer\n- {args.size_latent} latent dimensions"
    )

    # Train the model
    model.train()

    # Save model
    model_path = os.path.join(run_dir, f"scVI_{now}")
    model.save(model_path)

    # -----------------------
    #     Extract results
    # -----------------------

    # Save latent representation
    adata.obsm["X_scVI"] = model.get_latent_representation()
    logger.info(f"Saved latent representation to .h5ad")

    # Save normalized values
    adata.layers["scVI_normalized"] = model.get_normalized_expression()
    logger.info(f"Saved scVI normalized values to .h5ad")

    # Save newly generated data
    adata_path = Path(os.path.join(run_dir, os.path.basename(data_file)))
    adata.write(adata_path)
    logger.info(f"💾 Saved AnnData object to: {os.path.basename(data_file)}")

    # -------------------------
    #     Visualize results
    # -------------------------

    # Save figures in run directory
    sc.settings.figdir = run_dir

    # Visualize raw counts
    visualize_umap(adata, "raw", args.batch_effects)

    # Visualize corrected counts
    visualize_umap(adata, "X_scVI", args.batch_effects)


if __name__ == "__main__":
    main()
