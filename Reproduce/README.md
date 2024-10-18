# Data Processing and Visualization

This repository contains scripts for processing and visualizing various omics of datasets: Brain5k, Oihane, Gutierrez, Quake, and Spleen datasets.

## Contents

1. `process_datasets.py`: Python script for processing multiple types of datasets.
2. `Reproduce_Brain5k.ipynb`: Jupyter notebook for reproducing Brain5k dataset figures.
3. `Reproduce_Gutierrez.ipynb`: Jupyter notebook for reproducing Gutierrez dataset figures.
4. `Reproduce_hbrc.ipynb`: Jupyter notebook for reproducing HBRC dataset figures.
5. `Reproduce_Oihane.ipynb`: Jupyter notebook for reproducing Oihane dataset figures.
6. `Reproduce_Quake.ipynb`: Jupyter notebook for reproducing Quake dataset figures.
7. `Reproduce_Spleen.ipynb`: Jupyter notebook for reproducing Spleen dataset figures.

## Requirements

- Python 3.7+
- pandas
- anndata
- scanpy
- numpy
- matplotlib
- torch
- SpatialGlue

## Data

Users can either:
1. Use the pre-processed data in the "../dataset" folder directly, or
2. Download the original data from [link] and run `process_datasets.py` to process and clean the data.

The script assumes that the following data files are present in the "../dataset" directory:

- `10x-ATAC-Brain5k.h5ad`
- `Oihane.h5ad`
- `Gutierrez.h5ad`
- `Quake_Smart-seq2_Lung.h5ad`
- `adata_RNA.h5ad`
- `adata_Pro.h5ad`

## Usage

### Data Processing

To process all datasets, run the `process_datasets.py` script, and this script will process the following datasets:

- Brain5k
- Oihane
- Gutierrez
- Quake
- Spleen

The processed data will be saved as CSV files in the "../dataset" directory.

### Figure Reproduction

To reproduce the figures from the article, open and run each corresponding Jupyter notebook:

- `Reproduce_Brain5k.ipynb`
- `Reproduce_Gutierrez.ipynb`
- `Reproduce_hbrc.ipynb`
- `Reproduce_Oihane.ipynb`
- `Reproduce_Quake.ipynb`
- `Reproduce_Spleen.ipynb`
