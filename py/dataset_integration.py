# Standard imports
import os
import sys
import tempfile

# Third-party imports
import scanpy as sc
import scvi
import seaborn as sns
import torch

sc.set_figure_params(figsize=(6, 6), frameon=False)
sns.set_theme()
scvi.settings.verbosity = 10
torch.set_float32_matmul_precision("high")
save_dir = tempfile.TemporaryDirectory()


def main(path_to_data: str):
    """
    Runs the scVI model on the spatial liver dataset

    Args:
        path_to_data: Path to data source.
    """
    if not path_to_data:
        raise ValueError("Path to data must be provided")

    # Load data
    adata = scvi.data.read_h5ad(path_to_data)
    print("AnnData loaded successfully")
    print(adata)

    scvi.model.SCVI.setup_anndata(adata, batch_key=["Slide_name"])
    # Standard values from tutorial
    # https://docs.scvi-tools.org/en/stable/tutorials/notebooks/scrna/harmonization.html
    # Gene likelihood might need to change
    model = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
    model.train()

    model_dir = os.path.join(save_dir.name, "scvi_model")
    model.save(model_dir, overwrite=True)

    SCVI_LATENT_KEY = "X_scVI"
    adata.obsm[SCVI_LATENT_KEY] = model.get_latent_representation()

    # Cluster dataset and visualize latent space
    sc.pp.neighbors(adata, use_rep=SCVI_LATENT_KEY)
    sc.tl.leiden(adata)

    sc.tl.umap(adata, min_dist=0.3)
    sc.pl.umap(
        adata,
        color=["Slide_name", "leiden"],
        frameon=False,
        ncols=1,
    )


if __name__ == "__main__":
    main(sys.argv[1])
