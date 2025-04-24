import scanpy as sc

# CellRanger liver data in mtx format
mouse_liver_all_samples = "../data/training_data/countTable_mouseStSt"

adata = sc.read_10x_mtx(mouse_liver_all_samples, make_unique=True)

with open("integrate_training_data.log", "w") as f:
    f.write(f"Shape of the AnnData object: {adata.shape}\n")
    f.write(adata.obs.head())
    f.write("\n")
    f.write(adata.var.head())
    f.write("\n")

sc.pp.calculate_qc_metrics(adata, inplace=True)
sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts"],
    jitter=0.4,
    multi_panel=True,
    save="mla_stst.png",
)

adata.write("mouse_liver_atlas_stst.h5ad")
