# Standard imports
from datetime import datetime
import os
import sys

# Third-party imports
import scanpy as sc
import scvi
import seaborn as sns
import torch

# Set options
scvi.settings.verbosity = 10
torch.set_float32_matmul_precision("high")
sc.set_figure_params(figsize=(6, 6), frameon=False)
sns.set_theme()
save_dir = "./"


def main(path_to_data: str):
    """
    Runs the scVI model on the spatial liver dataset

    Args:
        path_to_data: Path to data source.
    """
    if not path_to_data:
        raise ValueError("Path to data must be provided")

    print(">>> Start")

    # Load data
    adata = scvi.data.read_h5ad(path_to_data)
    print("AnnData loaded successfully")
    print(adata)

    # Setup the data for processing
    # Model should correct for slide and condition
    scvi.model.SCVI.setup_anndata(
        adata, categorical_covariate_keys=["Slide_name", "condition"]
    )

    # Standard values from tutorial
    # https://docs.scvi-tools.org/en/stable/tutorials/notebooks/scrna/harmonization.html
    model = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
    print(model)

    model.train()

    # Save model
    model_dir = os.path.join(
        save_dir, f"scvi_model_{datetime.now().strftime('%d-%m-%Y_%H:%M')}"
    )
    model.save(model_dir, overwrite=True)

    # Save latent representation
    SCVI_LATENT_KEY = "X_scVI"

    latent = model.get_latent_representation()
    adata.obsm[SCVI_LATENT_KEY] = latent
    print(latent.shape)

    # Calculate UMAP for non-corrected data
    # Visualize afterwards
    sc.tl.pca(adata)
    sc.pp.neighbors(adata, n_pcs=30, n_neighbors=20)
    sc.tl.umap(adata, min_dist=0.3)

    sc.pl.umap(
        adata,
        color=["Slide_name"],
        frameon=False,
        save="UMAP_slide_name_not_corrected.svg",
        show=False,
    )
    sc.pl.umap(
        adata,
        color=["condition"],
        frameon=False,
        save="UMAP_condition_not_corrected.svg",
        show=False,
    )

    # Calculate UMAP for corrected data
    # Visualize afterwards
    sc.pp.neighbors(adata, use_rep=SCVI_LATENT_KEY)
    sc.tl.umap(adata, min_dist=0.3)

    sc.pl.umap(
        adata,
        color=["Slide_name"],
        frameon=False,
        save="UMAP_slide_name_corrected.svg",
        show=False,
    )
    sc.pl.umap(
        adata,
        color=["condition"],
        frameon=False,
        save="UMAP_condition_corrected.svg",
        show=False,
    )

    print(">>> Done!")


if __name__ == "__main__":
    main(sys.argv[1])
