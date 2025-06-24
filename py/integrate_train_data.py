import argparse
from datetime import datetime
import gzip

import scanpy as sc
import pandas as pd
from scipy import io

from utils import utils


def convert_training_data():

    # -------------
    #     Setup
    # -------------

    # Setup parser
    parser = argparse.ArgumentParser(
        description="Convert CellRanger mouse liver training data to AnnData."
    )

    parser.add_argument(
        "--data", "-d", type=str, required=True, help="Path to the data directory."
    )
    parser.add_argument(
        "--annotation",
        "-a",
        type=str,
        required=True,
        help="Path to the data annotation file. Must be '.csv'.",
    )
    parser.add_argument(
        "--output", "-o", type=str, required=True, help="Path to the output directory."
    )

    # Parse arguments
    args = parser.parse_args()

    data_dir = utils.assert_path(args.data)
    anno_file = utils.assert_path(args.annotation, assert_dir=False)
    out_dir = utils.assert_path(args.output)

    # Create directory for results
    now = datetime.now().strftime("%Y-%m-%d_%H-%M")
    results_dir = out_dir / now
    results_dir.mkdir(parents=True, exist_ok=True)

    # Setup logging
    logging_path = results_dir / "train_data_integration.log"
    logger = utils.setup_logging(logging_path)
    logger.info("Setup complete.")

    # --------------------
    #     Prepare data
    # --------------------

    logger.info("Start integration of training data...")

    # CellRanger outputs in directory
    # Stored as `.mtx` format
    matrix = io.mmread(
        data_dir / "matrix.mtx.gz"
    ).tocsc()  # Converst to csr when transposed
    logger.info(f"Loaded data matrix of shape: {matrix.shape}")

    # Load barcodes
    with gzip.open(data_dir / "barcodes.tsv.gz", "rt") as f:
        barcodes = [line.strip() for line in f]
    logger.info(f"Successfully loaded {len(barcodes)} observations.")

    # Load gene names
    with gzip.open(data_dir / "features.tsv.gz", "rt") as f:
        genes = [line.strip() for line in f]
    logger.info(f"Successfully loaded {len(genes)} feature names.")

    # Construct AnnData manually
    adata = sc.AnnData(X=matrix.T)
    adata.var_names = genes
    adata.obs_names = barcodes
    logger.info(f"Data was successfully converted to AnnData.\n{adata}")

    # ----------------------
    #     Add annotation
    # ----------------------

    # For this special case: All values not found in the annotation did not pass the QC
    # Remove them for better training results
    obs_anno = pd.read_csv(anno_file)
    logger.debug(f"Annotation matrix loaded successfully.\n{obs_anno}")

    # Set cell column as index that can be alligned
    obs_anno = obs_anno.set_index("cell")

    # Reorder annotation to match adata.obs_names
    obs_anno = obs_anno.reindex(adata.obs_names)

    # Assign annotation as observations to adata
    adata.obs = obs_anno
    logger.info(f"Annotation matrix added as observations to AnnData object.\n{adata}")

    # Save number of obs for user ouput
    init_obs = adata.n_obs
    # Drop all samples that did not pass the QC
    obs_mask = ~adata.obs.isna().all(axis=1)
    adata = adata[obs_mask].copy()
    logger.info(
        f"Removed low quality samples. Reduced number from {init_obs} to {adata.n_obs}."
    )

    # ------------
    #     Save
    # ------------

    adata_path = results_dir / f"{results_dir.parts[-1]}.h5ad"
    adata.write(adata_path)
    logger.info(f"💾 Saved AnnData object to: {adata_path}")
    logger.info(f"Finished!")

    # ----------
    #     QC
    # ----------

    # Computes QC metrics for each cell and each gene
    sc.pp.calculate_qc_metrics(adata, inplace=True)

    sc.settings.figdir = results_dir
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts"],
        jitter=0.4,
        multi_panel=True,
        save="quality_metrics.png",
    )


if __name__ == "__main__":
    convert_training_data()
