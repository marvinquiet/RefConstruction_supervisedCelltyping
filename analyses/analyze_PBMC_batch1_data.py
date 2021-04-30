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

## ---  analyze PBMC batch 1
data_dir = "/home/wma36/gpu/data"
result_dir = "/home/wma36/data_analysis/PBMC_batch1_analysis"

from preprocess import load_PBMC_data
PBMC_batch1_df = load_PBMC_data.load_PBMC_batch1_df(data_dir)
## plot on celltype
plt.clf()
fig, ax = plt.subplots(figsize=(6, 6))
sns_plot = sns.scatterplot(data=PBMC_batch1_df,
        x="tsne1", y="tsne2",
        hue="cell.type", alpha=0.5, s=3)
sns_plot.set_title("PBMC batch 1 (by celltype) tSNE")
sns_plot.set_xlabel("tSNE1")
sns_plot.set_ylabel("tSNE2")
lgd = plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1),ncol=4)

fig.tight_layout()
fig.savefig(result_dir+os.sep+"PBMC_batch1_celltype_tSNE.png",
        bbox_extra_artists=(lgd,),
        bbox_inches='tight', dpi=300)
## plot on batch
plt.clf()
fig, ax = plt.subplots(figsize=(6, 6))
sns_plot = sns.scatterplot(data=PBMC_batch1_df,
        x="tsne1", y="tsne2",
        hue="batch", alpha=0.5, s=3)
sns_plot.set_title("PBMC batch 1 (by batch) tSNE")
sns_plot.set_xlabel("tSNE1")
sns_plot.set_ylabel("tSNE2")
lgd = plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1),ncol=4)

fig.tight_layout()
fig.savefig(result_dir+os.sep+"PBMC_batch1_ABC_tSNE.png",
        bbox_extra_artists=(lgd,),
        bbox_inches='tight', dpi=300)

## plot on S1(1154) and S5 (1085)
plt.clf()
PBMC_batch1_df['ind'] = PBMC_batch1_df['ind'].astype('category')
fig, ax = plt.subplots(figsize=(6, 6))
sns_plot = sns.scatterplot(data=PBMC_batch1_df,
        x="tsne1", y="tsne2",
        hue="ind", alpha=0.5, s=3)
sns_plot.set_title("PBMC batch 1 (by individual) tSNE")
sns_plot.set_xlabel("tSNE1")
sns_plot.set_ylabel("tSNE2")
lgd = plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1),ncol=4)

fig.tight_layout()
fig.savefig(result_dir+os.sep+"PBMC_batch1_ind_tSNE.png",
        bbox_extra_artists=(lgd,),
        bbox_inches='tight', dpi=300)

## PBMC data extract out S1 and S5 to do tSNE
PBMC_batch1_adata = load_PBMC_data.load_PBMC_batch1_data(data_dir+os.sep+"PBMC_batch1")
PBMC_batch1_df['barcode'] = PBMC_batch1_df['batch'] + '_' + PBMC_batch1_df.index

S1_df = PBMC_batch1_df[PBMC_batch1_df["ind"] == 1154]
S1_df = S1_df[S1_df["batch"] != "B"] # 1551 cells

S2_df = PBMC_batch1_df[PBMC_batch1_df["ind"] == 1085]
S2_df = S2_df[S2_df["batch"] != "A"] # 1583 cells

S1S2_barcodes = S1_df["barcode"].tolist() + S2_df["barcode"].tolist() # 3134 cells
S1S2_adata = PBMC_batch1_adata[S1S2_barcodes]
# tSNE on data
tmp_adata = sc.AnnData(S1S2_adata.X.toarray())
tmp_adata.obs = S1S2_adata.obs
tmp_adata.obs.index.name = None
obs_df = tmp_adata.obs.merge(PBMC_batch1_df, on="barcode", left_index=True)
obs_df.index = obs_df["barcode"]
tmp_adata.obs = obs_df
## do PCA first
sc.tl.pca(tmp_adata, 
        random_state=RANDOM_SEED)

## do tSNE on PCA result
sc.tl.tsne(tmp_adata, n_pcs=50, use_rep="X_pca",
        learning_rate=300,
        perplexity=30,
        n_jobs=4,
        random_state=RANDOM_SEED)
obsm_data = pd.DataFrame(tmp_adata.obsm["X_tsne"])
obsm_data.columns = ['tsne1', 'tsne2']
obsm_data.index = tmp_adata.obs.index
obsm_data['cell.type'] = tmp_adata.obs['cell.type']
obsm_data['ind'] = tmp_adata.obs['ind']
obsm_data.to_csv(result_dir+os.sep+"tSNE_on_S1_S5.csv")
## plot tSNE on cell type
plt.clf()
fig, ax = plt.subplots(figsize=(6, 6))
sns_plot = sns.scatterplot(data=obsm_data,
        x="tsne1", y="tsne2",
        hue="cell.type", alpha=0.5, s=3)
sns_plot.set_title("S1 and S5 t-SNE (celltype)")
sns_plot.set_xlabel("tSNE1")
sns_plot.set_ylabel("tSNE2")
lgd = plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1),ncol=4)

fig.tight_layout()
fig.savefig(result_dir+os.sep+"tSNE_on_S1_S5_celltype.png",
        bbox_extra_artists=(lgd,),
        bbox_inches='tight', dpi=300)

## plot tSNE on cell type
plt.clf()
obsm_data["ind"] = obsm_data["ind"].astype("category")
fig, ax = plt.subplots(figsize=(6, 6))
sns_plot = sns.scatterplot(data=obsm_data,
        x="tsne1", y="tsne2",
        hue="ind", alpha=0.5, s=3)
sns_plot.set_title("S1 and S5 t-SNE (individual)")
sns_plot.set_xlabel("tSNE1")
sns_plot.set_ylabel("tSNE2")
lgd = plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1),ncol=4)

fig.tight_layout()
fig.savefig(result_dir+os.sep+"tSNE_on_S1_S5_ind.png",
        bbox_extra_artists=(lgd,),
        bbox_inches='tight', dpi=300)
