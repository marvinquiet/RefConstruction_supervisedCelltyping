import os, sys
import pandas as pd
import numpy as np
import scanpy.api as sc
import anndata
from numpy.random import seed

# plot
import seaborn as sns
import matplotlib.pyplot as plt

RANDOM_SEED=0
seed(RANDOM_SEED)
np.random.seed(RANDOM_SEED)

## --- analyze mousebrain FC 
data_dir = "/home/wma36/gpu/data"

from preprocess import load_mousebrain_data
adata = load_mousebrain_data.load_brain_adata(data_dir, region="FC", ind=None, celltype_gran=0)

from preprocess import process_train_test_data
filtered_adata = process_train_test_data.process_adata(adata)

## use Seurat to select genes
sc.pp.highly_variable_genes(filtered_adata, n_top_genes=1000, subset=True)

## do PCA
sc.tl.pca(filtered_adata, random_state=RANDOM_SEED)
## do tSNE on PCA result
sc.tl.tsne(filtered_adata, n_pcs=50, use_rep="X_pca",
        learning_rate=300,
        perplexity=30,
        n_jobs=4,
        random_state=RANDOM_SEED)
sc.pl.tsne(filtered_adata, color=['ind', 'cell.type', 'cluster'], legend_fontsize="xx-small")
plt.savefig(data_dir+os.sep+'Mousebrain'+os.sep+'tSNE.png')

## do UMAP
sc.pp.neighbors(filtered_adata, n_neighbors=20, use_rep="X_pca", random_state=RANDOM_SEED) 
sc.tl.umap(filtered_adata, random_state=RANDOM_SEED)
sc.pl.umap(filtered_adata, color=['ind', 'cell.type', 'cluster'], legend_fontsize="xx-small")
plt.savefig(data_dir+os.sep+'Mousebrain'+os.sep+'umap.png')

