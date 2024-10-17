# ADMReproduce

This R script (main.R) performs various dimensionality reduction and visualization methods on single-cell RNA-seq data, followed by ensemble and meta-analysis techniques. It also calculates evaluation metrics for each method.

## Requirements

Before running this script, ensure you have the following:

- R (version 4.3.1 or later)
- The ADM R package (available at [Seven595/ADM](https://github.com/Seven595/ADM))

- Other related packages

  | Package      | Version    |
  | ------------ | ---------- |
  | gridExtra    | 2.3        |
  | stats        | 4.3.1      |
  | MASS         | 7.3.60.0.1 |
  | rARPACK      | 0.11.0     |
  | utils        | 4.3.1      |
  | grDevices    | 4.3.1      |
  | graphics     | 4.3.1      |
  | methods      | 4.3.1      |
  | phateR       | 1.0.7      |
  | ggplot2      | 3.5.1      |
  | mclust       | 6.1.1      |
  | cluster      | 2.1.6      |
  | magrittr     | 2.0.3      |
  | tidyr        | 1.3.1      |
  | dplyr        | 1.1.4      |
  | ggrepel      | 0.9.6      |
  | BiocParallel | 1.36.0     |
  | uwot         | 0.2.2      |
  | aricode      | 1.0.3      |
  | DDoutlier    | 0.1.0      |
  | diffudist    | 1.0.1      |
  | fitdistrplus | 1.2.1      |
  | igraph       | 2.0.3      |
  | Rtsne        | 0.17       |
  | dimRed       | 0.2.6      |

## Usage

1. Set the dataset name in the `dataset` variable. Available options are:
   - Gutierrez
   - Oihane
   - Quake
   - Brain5k
   - mir
   - Spleen
   - metabolism
   - gene
   
2. Run the script in R or RStudio.

## Script Overview

The script performs the following steps:

1. Loads necessary libraries and custom functions
2. Defines a color list for visualization
3. Loads the specified dataset
4. Executes various candidate visualization methods:
   - PCA
   - MDS
   - iMDS
   - Sammon
   - HLLE
   - Isomap
   - kPCA
   - LEIM
   - UMAP
   - tSNE (with perplexities 10 and 30)
   - PHATE
   - KEF
5. Executes the meta-spec and ADM visualization
7. Processes and visualizes meta-method results
8. Visualizes individual method results
9. Outputs a summary of results, including ARI, NMI, and Silhouette coefficient for each method

## Output

The script generates visualizations for each method and calculates evaluation metrics. Results can be accessed through the `ind_result` list, where each element contains:

- `plot`: A ggplot object of the visualization
- `ari`: Adjusted Rand Index
- `nmi`: Normalized Mutual Information
- `silhouette`: Silhouette coefficient

A summary of these metrics is printed to the console for easy comparison.

## Repository Structure

- `datasets/`: Contains partial datasets used in our study, all dataset used in our study can be downloaded from https://drive.google.com/drive/folders/1s9CtUIZYY5iwkHN-vOlrfcqDvCa4m8p1?usp=drive_link
- `Reproduce/`: Jupyter notebooks for reproducing figures and analyses.

  
