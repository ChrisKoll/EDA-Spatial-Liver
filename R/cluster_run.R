# Installation of SeuratDisk
# if (!requireNamespace("remotes", quietly = TRUE)) {
#   install.packages("remotes")
# }
# remotes::install_github("mojaveazure/seurat-disk")

# Load packages
library(dplyr)
library(Seurat)
library(SeuratDisk)
library(SeuratObject)
library(tibble)

# Load data
sobj_liver1 <- SeuratObject::LoadSeuratRds("/projects/spatial_liver/work/liver1_seurat_object.Rds")
sobj_liver2 <- SeuratObject::LoadSeuratRds("/projects/spatial_liver/work/liver2_seurat_object.Rds")

# ---------------------- #
# Add metadata to liver1 #
# ---------------------- #

# Disease condition annotation for liver1
healthy_l1 <- c(1, 2, 15, 16, 17, 18, 19, 20, 21, 22, 23, 45)
cirrhosis_l1 <- c(5, 6, 24, 25, 26, 30, 31, 34, 35, 40, 41, 44)
steatosis_l1 <- c(3, 4, 9, 10, 13, 14, 27, 28, 29, 38, 39)
aclf_l1 <- c(7, 8, 11, 12, 32, 33, 36, 37, 42, 43)

# Create condition tibble to add to metadata
conditions_l1 <- sobj_liver1@meta.data$fov |>
  tibble::as_tibble() |>
  # Sets condition according to look up vectors
  dplyr::mutate(
    condition = dplyr::case_when(
      value %in% healthy_l1 ~ "healthy",
      value %in% cirrhosis_l1 ~ "cirrhosis",
      value %in% steatosis_l1 ~ "steatosis",
      value %in% aclf_l1 ~ "aclf",
      TRUE ~ "Other"
    )
  ) |>
  # Do not need value column
  dplyr::select(-value)

# Add to metadata
sobj_liver1 <- Seurat::AddMetaData(
  sobj_liver1,
  metadata = conditions_l1
)

# ---------------------- #
# Add metadata to liver2 #
# ---------------------- #

cirrhosis_l2 <- c(1, 2, 3, 4, 5, 6, 7, 8, 13, 14, 17, 18, 19, 20, 21, 22, 23, 24, 25, 36, 37, 44, 45)
aclf_l2 <- c(9, 10, 11, 12, 15, 16, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 38, 39, 40, 41, 42, 43)

conditions_l2 <- sobj_liver2@meta.data$fov |>
  tibble::as_tibble() |>
  dplyr::mutate(
    condition = dplyr::case_when(
      value %in% cirrhosis_l2 ~ "cirrhosis",
      value %in% aclf_l2 ~ "aclf",
      TRUE ~ "Other"
    )
  ) |>
  dplyr::select(-value)

sobj_liver2 <- Seurat::AddMetaData(
  sobj_liver2,
  metadata = conditions_l2
)

# -------------------- #
# Merge Seurat objects #
# -------------------- #

sobj_liver_merged <- merge(
  sobj_liver1,
  y = sobj_liver2,
  add.cell.ids = c("liver1", "liver2")
)

SeuratDisk::SaveH5Seurat(sobj_liver_merged, "../data/liver_seurat_object.h5Seurat")
SeuratDisk::Convert("../data/liver_seurat_object.h5Seurat", dest = "h5ad")
