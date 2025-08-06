"""Convert CellRanger output to AnnData format with preprocessing and quality control.

This module provides functionality to convert CellRanger output files (matrix.mtx.gz,
barcodes.tsv.gz, features.tsv.gz) to AnnData format, preprocess the data, and generate
quality control metrics.
"""

import argparse
from datetime import datetime
import gzip
from logging import Logger
from pathlib import Path

from anndata import AnnData
import scanpy as sc
import pandas as pd
from scipy import io

from utils import utils


# Constants
DATETIME_FORMAT = "%Y-%m-%d_%H-%M"
LOG_FILENAME = "train_data_integration.log"
MATRIX_FILENAME = "matrix.mtx.gz"
BARCODES_FILENAME = "barcodes.tsv.gz"
FEATURES_FILENAME = "features.tsv.gz"
CELL_INDEX_COL = "cell"
RAW_LAYER_NAME = "raw"
LOG1P_LAYER_NAME = "log1p"
COUNTS_PER_SAMPLE_COL = "counts_per_sample"
SPARSITY_COL = "sparsity"
COUNTS_PER_GENE_COL = "counts_per_gene"
ANNOTATION_COL = "annot"
UMAP_1_COL = "UMAP_1"
UMAP_2_COL = "UMAP_2"
CLUSTER_COL = "cluster"


def convert_CellRanger_to_h5ad(
    data_dir: Path,
    anno_file: Path | None,
    logger: Logger,
) -> AnnData:
    """Convert CellRanger output files to AnnData format.

    Reads CellRanger output files (matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz)
    and converts them to an AnnData object. Optionally adds cell annotations and
    filters out low-quality samples that did not pass QC.

    :param data_dir: Path to directory containing CellRanger output files
    :type data_dir: Path
    :param anno_file: Path to CSV file containing cell annotations, or None if no annotations
    :type anno_file: Path | None
    :param logger: Logger instance for logging progress and information
    :type logger: Logger

    :return: AnnData object containing the converted single-cell data
    :rtype: AnnData

    :raises FileNotFoundError: If required CellRanger files are not found in data_dir
    :raises ValueError: If matrix dimensions don't match barcode/feature counts
    """
    logger.info("Start integration of training data...")

    # CellRanger outputs in directory
    # Matrix is stored in Matrix Market format (.mtx)
    matrix = io.mmread(
        data_dir / MATRIX_FILENAME
    ).tocsc()  # Converst to csr when transposed
    logger.info(f"📁 Loaded data matrix of shape: {matrix.shape}")

    # Load cell barcodes (observation names)
    with gzip.open(data_dir / BARCODES_FILENAME, "rt") as f:
        barcodes = [line.strip() for line in f]
    logger.info(f"📁 Loaded {len(barcodes)} observations.")

    # Load gene/feature names (variable names)
    with gzip.open(data_dir / FEATURES_FILENAME, "rt") as f:
        genes = [line.strip() for line in f]
    logger.info(f"📁 Loaded {len(genes)} feature names.")

    # Construct AnnData object
    # Note: Matrix needs to be transposed because CellRanger stores genes x cells
    # but AnnData expects cells x genes
    adata = sc.AnnData(X=matrix.T)
    adata.var_names = genes
    adata.obs_names = barcodes
    logger.info(f"Data was successfully converted to AnnData.\n{adata}")

    if anno_file is not None:
        # For this special case: All values not found in the annotation did not pass the QC
        # Remove them for better training results
        obs_anno = pd.read_csv(anno_file)
        logger.debug(f"📁 Annotation matrix loaded successfully.\n{obs_anno}")

        # Set cell column as index that can be alligned
        obs_anno = obs_anno.set_index(CELL_INDEX_COL)

        # Reorder annotation to match adata.obs_names
        obs_anno = obs_anno.reindex(adata.obs_names)

        # Assign annotation as observations to adata
        adata.obs = obs_anno
        logger.info(
            f"Annotation matrix added as observations to AnnData object.\n{adata}"
        )

        # Save number of obs for user ouput
        init_obs = adata.n_obs
        # Drop all samples that did not pass the QC
        # Assumption: Cells not in annotation file failed QC and should be removed
        obs_mask = ~adata.obs.isna().all(axis=1)
        adata = adata[obs_mask].copy()
        logger.info(
            f"Removed low quality samples. Reduced number from {init_obs} to {adata.n_obs}."
        )

    return adata


def preprocess_data(adata: AnnData, results_dir: Path, logger: Logger) -> AnnData:
    """Preprocess single-cell data and calculate quality metrics.

    Adds different data layers (raw counts, log1p transformed), calculates
    quality control metrics for both cells and genes, and saves the processed
    AnnData object to disk.

    :param adata: AnnData object containing single-cell data
    :type adata: AnnData
    :param results_dir: Directory path where processed data will be saved
    :type results_dir: Path
    :param logger: Logger instance for logging progress
    :type logger: Logger

    :return: Preprocessed AnnData object with added layers and metrics
    :rtype: AnnData
    """
    # Store original raw counts in a separate layer for future reference
    adata.layers[RAW_LAYER_NAME] = adata.X.copy()
    logger.info("Added original counts to layer 'raw'.")

    # Add log1p transformed counts (log(x + 1)) which is commonly used in scRNA-seq
    adata.layers[LOG1P_LAYER_NAME] = sc.pp.log1p(adata.X, copy=True)
    logger.info("Added log1p scaled counts to layer 'log1p'.")

    # Calculate cell-level (observation) quality metrics
    # Total read counts per cell
    counts_per_sample = adata.X.sum(axis=1)
    adata.obs["counts_per_sample"] = counts_per_sample.A1

    # Sparsity: fraction of genes with zero expression in each cell
    sample_sparsity = 1 - (adata.X.getnnz(axis=1) / adata.n_vars)
    adata.obs["sparsity"] = sample_sparsity
    logger.info("Added sample specific metrics to 'obs'.")

    # Calculate gene-level (variable) quality metrics
    # Convert to CSC format for efficient column operations
    csc_matrix = adata.X.tocsc()

    # Total counts per gene across all cells
    counts_per_gene = csc_matrix.sum(axis=0)
    adata.var["counts_per_gene"] = counts_per_gene.A1

    # Gene sparsity: fraction of cells with zero expression for each gene
    gene_sparsity = 1 - (csc_matrix.getnnz(axis=0) / adata.n_obs)
    adata.var["sparsity"] = gene_sparsity
    logger.info("Added gene specific metrics to 'var'.")
    logger.info(f"{adata}")

    # Save processed AnnData object to disk
    # Use directory name as filename for consistency
    adata_path = results_dir / f"{results_dir.parts[-1]}.h5ad"
    adata.write(adata_path, compression="gzip")
    logger.info(f"💾 Saved AnnData object to: '{adata_path}'.")

    return adata


def run_qc(adata: AnnData, results_dir: Path, logger: Logger):
    """Generate and save quality control metrics and reports.

    Creates separate files containing cell-level and gene-level quality metrics
    in Feather format for efficient storage and retrieval. Also logs basic
    statistics about the data normalization.

    :param adata: AnnData object containing processed single-cell data
    :type adata: AnnData
    :param results_dir: Directory path where QC files will be saved
    :type results_dir: Path
    :param logger: Logger instance for logging QC information
    :type logger: Logger
    """
    # Create and save cell-level metrics DataFrame
    sample_metrics = pd.DataFrame(
        {
            "counts": adata.obs[COUNTS_PER_SAMPLE_COL],
            "sparsity": adata.obs[SPARSITY_COL],
            "annotation": adata.obs[ANNOTATION_COL],
            "umap_1": adata.obs[UMAP_1_COL],
            "umap_2": adata.obs[UMAP_2_COL],
            "cluster": adata.obs[CLUSTER_COL],
        },
        index=adata.obs_names,
    )

    save_path = results_dir / f"{results_dir.parts[-1]}_sample_metrics.feather"
    sample_metrics.to_feather(save_path)
    logger.info(f"💾 Saved sample specific metrics in: '{save_path}'.")

    # Create and save gene-level metrics DataFrame
    gene_metrics = pd.DataFrame(
        {
            "counts": adata.var[COUNTS_PER_GENE_COL],
            "sparsity": adata.var[SPARSITY_COL],
        },
        index=adata.var_names,
    )

    save_path = results_dir / f"{results_dir.parts[-1]}_gene_metrics.feather"
    gene_metrics.to_feather(save_path)
    logger.info(f"💾 Saved gene specific metrics in: '{save_path}'.")

    # Log basic data statistics for quality assessment
    logger.info("---")
    logger.info(f"Range of counts: {adata.X.min()} - {adata.X.max()}.")
    logger.info(f"Mean expression: {adata.X.mean()}.")
    logger.info("---")

    logger.info(f"Finished!")


def main():
    """Main function to orchestrate the CellRanger to AnnData conversion pipeline.

    Parses command line arguments, sets up logging, and executes the three main
    pipeline steps: conversion, preprocessing, and quality control.
    """
    # Setup argument parser
    parser = argparse.ArgumentParser(
        description="Convert CellRanger mouse liver training data to AnnData."
    )

    parser.add_argument(
        "--data",
        "-d",
        type=str,
        required=True,
        help="Path to the data directory containing CellRanger output files.",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        required=True,
        help="Path to the output directory where results will be saved.",
    )
    parser.add_argument(
        "--annotation",
        "-a",
        type=str,
        help="Path to the data annotation file. Must be '.csv' format.",
    )

    # Parse command line arguments
    args = parser.parse_args()

    # Validate input paths
    data_dir = utils.assert_path(args.data)
    out_dir = utils.assert_path(args.output)

    if args.annotation is not None:
        anno_file = utils.assert_path(args.annotation, assert_dir=False)
    else:
        anno_file = None

    # Create timestamped results directory
    now = datetime.now().strftime(DATETIME_FORMAT)
    results_dir = out_dir / f"run_{now}"
    results_dir.mkdir(parents=True, exist_ok=True)

    # Setup logging
    logging_path = results_dir / LOG_FILENAME
    logger = utils.setup_logging(logging_path)
    logger.info("✅ Setup complete.")
    logger.info("----")

    # Execute pipeline steps
    adata = convert_CellRanger_to_h5ad(data_dir, anno_file, logger)
    adata = preprocess_data(adata, results_dir, logger)
    run_qc(adata, results_dir, logger)


if __name__ == "__main__":
    main()
