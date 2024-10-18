# Import necessary libraries
import pandas as pd
import anndata as ad
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import os
import torch
import SpatialGlue
from SpatialGlue.preprocess import fix_seed, clr_normalize_each_cell, pca, construct_neighbor_graph
from SpatialGlue.SpatialGlue_pyG import Train_SpatialGlue

# Set device for PyTorch
device = torch.device('cuda:1' if torch.cuda.is_available() else 'cpu')

# Function to process Brain5k dataset
def process_brain5k():
    # Load the dataset
    adata = ad.read_h5ad("../dataset/10x-ATAC-Brain5k.h5ad")
    
    # Preprocess the data
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat')
    
    # Filter for highly variable genes
    adata_hvg = adata[:, adata.var.highly_variable].copy()
    
    # Convert to pandas DataFrame and save
    pbmc_5k_data = pd.DataFrame(adata_hvg.X.toarray(), index=adata_hvg.obs.index, columns=adata_hvg.var.index)
    pbmc_5k_data.to_csv("../dataset/pbmc_5k_data.csv")
    
    # Save metadata
    metadata = adata_hvg.obs
    metadata.to_csv("../dataset/pbmc_5k_metadata.csv")

# Function to process Oihane dataset
def process_oihane():
    # Load the dataset
    adata = ad.read_h5ad("../dataset/Oihane.h5ad")
    
    # Save raw data
    Oihane = pd.DataFrame(adata.X.toarray(), index=adata.obs.index, columns=adata.var.index)
    Oihane.to_csv("../dataset/Oihane_raw.csv")
    
    # Filter genes
    sc.pp.filter_genes(adata, min_cells=int(adata.n_obs * 0.2))
    
    # Process and save filtered data
    Oihane_filtered = pd.DataFrame(adata.X.toarray(), index=adata.obs.index, columns=adata.var.index)
    Oihane_filtered_log = np.log10(Oihane_filtered + 1)
    Oihane_filtered_log.to_csv("../dataset/Oihane_data.csv")
    
    # Save metadata
    adata.obs.to_csv("../dataset/Oihane_metadata.csv")

# Function to process Gutierrez dataset
def process_gutierrez():
    # Load the dataset
    adata = ad.read_h5ad("../dataset/Gutierrez.h5ad")
    
    # Save raw data
    Gutierrez = pd.DataFrame(adata.X.toarray(), index=adata.obs.index, columns=adata.var.index)
    Gutierrez.to_csv("../dataset/Gutierrez_raw.csv")
    
    # Filter genes
    sc.pp.filter_genes(adata, min_cells=int(adata.n_obs * 0.30))
    
    # Process and save filtered data
    Gutierrez_filtered = pd.DataFrame(adata.X.toarray(), index=adata.obs.index, columns=adata.var.index)
    Gutierrez_filtered_log = np.log10(Gutierrez_filtered + 1)
    Gutierrez_filtered_log.to_csv("../dataset/Gutierrez_data.csv")
    
    # Save metadata
    adata.obs.to_csv("../dataset/Gutierrez_metadata.csv")

# Function to process Quake dataset
def process_quake():
    # Load the dataset
    adata = ad.read_h5ad("../dataset/Quake_Smart-seq2_Lung.h5ad")
    
    # Save raw data
    Quake = pd.DataFrame(adata.X.toarray(), index=adata.obs.index, columns=adata.var.index)
    Quake.to_csv("../dataset/Quake_raw.csv")
    
    # Filter genes
    sc.pp.filter_genes(adata, min_cells=int(adata.n_obs * 0.30))
    
    # Process and save filtered data
    Quake_filtered = pd.DataFrame(adata.X.toarray(), index=adata.obs.index, columns=adata.var.index)
    Quake_filtered_log = np.log10(Quake_filtered + 1)
    Quake_filtered_log.to_csv("../dataset/Quake_data.csv")
    
    # Save metadata
    adata.obs.to_csv("../dataset/Quake_metadata.csv")

# Function to process Spleen dataset
def process_spleen():
    # Set random seed
    random_seed = 2022
    fix_seed(random_seed)
    
    # Load datasets
    adata_omics1 = sc.read_h5ad('../dataset/adata_RNA.h5ad')
    adata_omics2 = sc.read_h5ad('../dataset/adata_Pro.h5ad')
    
    # Make var names unique
    adata_omics1.var_names_make_unique()
    adata_omics2.var_names_make_unique()
    
    # Set data type
    data_type = 'SPOTS'
    
    # Preprocess adata_omics1 (RNA)
    sc.pp.filter_genes(adata_omics1, min_cells=10)
    sc.pp.highly_variable_genes(adata_omics1, flavor="seurat_v3", n_top_genes=3000)
    sc.pp.normalize_total(adata_omics1, target_sum=1e4)
    sc.pp.log1p(adata_omics1)
    sc.pp.scale(adata_omics1)
    
    # Calculate PCA for adata_omics1
    adata_omics1_high = adata_omics1[:, adata_omics1.var['highly_variable']]
    adata_omics1.obsm['feat'] = pca(adata_omics1_high, n_comps=adata_omics2.n_vars-1)
    
    # Preprocess adata_omics2 (Protein)
    adata_omics2 = clr_normalize_each_cell(adata_omics2)
    sc.pp.scale(adata_omics2)
    adata_omics2.obsm['feat'] = pca(adata_omics2, n_comps=adata_omics2.n_vars-1)
    
    # Construct neighbor graph
    data = construct_neighbor_graph(adata_omics1, adata_omics2, datatype=data_type)
    
    # Train SpatialGlue model
    model = Train_SpatialGlue(data, datatype=data_type, device=device)
    output = model.train()
    
    # Process output
    adata = adata_omics1.copy()
    adata.obsm['emb_latent_omics1'] = output['emb_latent_omics1'].copy()
    adata.obsm['emb_latent_omics2'] = output['emb_latent_omics2'].copy()
    adata.obsm['SpatialGlue'] = output['SpatialGlue'].copy()
    adata.obsm['alpha'] = output['alpha']
    adata.obsm['alpha_omics1'] = output['alpha_omics1']
    adata.obsm['alpha_omics2'] = output['alpha_omics2']
    adata.obsm["spatialglue_copy"] = adata.obsm["SpatialGlue"]
    
    # Visualization
    fig, ax_list = plt.subplots(1, 2, figsize=(7, 3))
    sc.pp.neighbors(adata, use_rep='SpatialGlue', n_neighbors=10)
    sc.tl.umap(adata)
    sc.pl.umap(adata, color='SpatialGlue', ax=ax_list[0], title='SpatialGlue', s=20, show=False)
    sc.pl.embedding(adata, basis='spatial', color='SpatialGlue', ax=ax_list[1], title='SpatialGlue', s=25, show=False)
    plt.tight_layout(w_pad=0.3)
    
    # Annotation
    adata.obs['SpatialGlue_number'] = adata.obs['SpatialGlue'].copy()
    adata.obs['SpatialGlue'] = adata.obs['SpatialGlue'].cat.rename_categories({
        1: 'MMMØ', 2: 'T cell', 3: 'B cell', 4: 'MZMØ', 5: 'RpMØ'
    })
    list_ = ['MZMØ','MMMØ','RpMØ','B cell', 'T cell']
    adata.obs['SpatialGlue'] = pd.Categorical(adata.obs['SpatialGlue'], categories=list_, ordered=True)
    
    # Save results
    spleen_label = pd.DataFrame(adata.obs["SpatialGlue"])
    spleen_pro = pd.DataFrame(adata_omics2.X, index=adata_omics2.obs.index, columns=adata_omics2.var.index)
    spleen_pro.to_csv("../dataset/Spleen_pro_data.csv")
    spleen_label.to_csv("../dataset/Spleen_pro_label.csv")

# Main execution
if __name__ == "__main__":
    process_brain5k()
    process_oihane()
    process_gutierrez()
    process_quake()
    process_spleen()