**Master’s Thesis Project – Bioinformatics**  
Goethe University Frankfurt

**Thesis Title:**  
*Neural Network-Based Gene Expression Imputation for Subcellular Spatial Transcriptomics Data*

For the companion repository focusing on model training and evaluation, see: [Liver-Expression-Imputation](https://github.com/CKolland/Liver-Expression-Imputation).

# EDA-Spatial-Liver

[![Python Version](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![R version](https://img.shields.io/badge/R-%3E%3D4.4.2-blue)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## 📂 Project Structure

The project workflow is organized into two main parts:

1. Exploratory Data Analysis (EDA) in *R*

    - Initial inspection and quality control of spatial transcriptomics datasets
    - Merging of two *Seurat* objects for integrated analysis
    
2. Transition to *Python* & Advanced Analysis

    - Conversion of merged *Seurat* objects into *AnnData* format for use in *Python*
    - Preprocessing and preparation for neural network–based model training
    - Integration with *scVI-tools* for representation learning and harmonization
    - Clustering and cell type annotation within the scverse framework

This hybrid approach leverages both *Seurat* (*R*) and *scverse* (*Python*), combining their strengths for high-resolution analysis of spatial transcriptomics data.

>⚠️ Note: Large raw and intermediate data files have been excluded from this repository due to size limitations. They can be provided upon request.


---

```
EDA-Spatial-Liver/
├── R/
│   ├── 1_IDE_mouse_liver_ST.qmd
│   ├── 2_merge_ST_Seurat_objects.qmd
│   ├── 3_EDA_mouse_liver_ST.qmd
│   ├── 4_extract_from_Seurat.qmd
│   └── helpers/
│       ├── render_plots.R
│       └── render_tables.R
└── py/
    ├── 1_assemble_h5ad.ipynb
    ├── 2_integrate_train_data.py
    ├── 2_visualize_train_integration.ipynb
    ├── 3_update_h5ad.ipynb
    ├── 4_integrate_ST_data.py
    ├── 4_visualize_ST_integration.ipynb
    ├── 5_cell_type_annotation.ipynb
    └── 6_export_results.ipynb
```

---

## ⚙️ Installation

```bash
git clone https://github.com/CKolland/EDA-Spatial-Liver.git
cd EDA-Spatial-Liver
```

### Setup R Environment

Open project in *R* with `EDA-Spatial-Liver.Rproj`

```R
library(renv)
renv::restore()
```

### Setup Python Environment

```bash
mamba env create -n <env_name> -f environment.yml
mamba activate <env_name>
```

## 📜 License

This project is licensed under the MIT License.   
See [LICENSE](LICENSE) for details.

## 📬 Contact

*Maintainer:* Christian Kolland ([Schulz Lab](https://schulzlab.github.io/))

For questions, requests (including access to large data files), or feedback, please contact the maintainer directly.
