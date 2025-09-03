# EDA-Spatial-Liver

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Exploratory data analysis of a spatial liver dataset.

---

## Project Structure Analysis

**Main Programming Language(s) and Frameworks:**  
- Python (Jupyter Notebooks, scanpy, anndata, pandas, numpy, matplotlib)

**Project Type:**  
- Data analysis toolkit for spatial transcriptomics, focused on assembling and exploring `.h5ad` (AnnData) files.

**Dependencies and Package Managers:**  
- Dependencies managed via `pip` (see `requirements.txt`).
- Jupyter Notebook for interactive analysis.

**Build Tools and Configuration Files:**  
- No build tools detected (not required for notebooks).
- Configuration via Python scripts and notebook cells.

---

## 1. Project Title and Description

**EDA-Spatial-Liver**  
A toolkit for assembling, processing, and exploring spatial transcriptomics data from liver tissue.  
It streamlines the creation of AnnData (`.h5ad`) objects from raw data, enabling efficient downstream analysis and visualization.

---

## 2. Features

- Assemble AnnData objects from raw spatial transcriptomics files
- Integrate count matrices, spatial coordinates, and metadata
- Quality control and preprocessing workflows
- Interactive exploratory data analysis in Jupyter Notebooks
- Visualization of spatial gene expression patterns
- Modular Python scripts for reproducibility

---

## 3. Installation

### Prerequisites

- Python 3.8 or higher
- pip
- Jupyter Notebook

### Step-by-Step Instructions

```bash
git clone https://github.com/yourusername/EDA-Spatial-Liver.git
cd EDA-Spatial-Liver
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
jupyter notebook
```

---

## 4. Usage

### Assembling AnnData Objects

Open `1_assemble_h5ad.ipynb` and follow the notebook instructions to assemble AnnData objects.

```python
import scanpy as sc
import anndata

# Load count matrix, spatial coordinates, and metadata
# ...see notebook for details...

adata = anndata.AnnData(X=counts, obs=metadata, obsm={"spatial": coordinates})
adata.write("liver_spatial.h5ad")
```

### Exploratory Data Analysis

Use the provided notebooks to visualize and analyze spatial gene expression.

```python
sc.pl.spatial(adata, color="GeneA")
```

---

## 5. Configuration

- `requirements.txt` for dependencies
- Paths and settings are set within notebook cells or Python scripts

---

## 6. Development

### Setting Up

- Follow installation steps above.
- Open notebooks in Jupyter for interactive development.

### Running Tests

If unit tests are present, run:

```bash
pytest
```

### Building the Project

- No build step required; notebooks and scripts are ready to use.

### Contributing

- Fork the repository and submit pull requests.
- Follow PEP8 style guidelines.
- Add docstrings and comments to new code.

---

## 7. Additional Sections

### License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

### Acknowledgments

- Built with scanpy, anndata, and the Python scientific stack.
- Inspired by open-source spatial transcriptomics
