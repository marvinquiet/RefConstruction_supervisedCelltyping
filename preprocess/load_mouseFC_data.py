import os
import anndata
import numpy as np
import pandas as pd

def load_FC_adata(data_dir, devstage=None, ind=None, treatment=None,
        celltype_gran=0, curate=False):
    '''Loading data from mouse frontal cortex along with metadata information

    @devstage: P21 or Adult
    @ind: individual sample prediction: P21Sample1-3, PFCSample1-12
    @treatment: Saline or Cocaine
    @celltype_gran: major cell types or sub-cell types
    '''
    ## load data; set genes and cells
    adata = anndata.read_csv(data_dir+os.sep+"MouseFC_GSE124952/GSE124952_expression_matrix.csv").T

    adata.var["gene_symbols"] = adata.var.index
    adata.var_names_make_unique(join="-")

    adata.obs['barcode'] = adata.obs.index
    adata.obs_names_make_unique(join="-")

    ## load metadata
    meta_df = pd.read_csv(data_dir+os.sep+"MouseFC_GSE124952/GSE124952_meta_data.csv", index_col=0)
    adata.obs = adata.obs.merge(meta_df, left_index=True, right_index=True)

    if 0 == celltype_gran:
        adata.obs.rename(columns={"CellType": "cell.type"}, inplace=True)  ## change cell type names

        if curate: ## if curation is needed for cross dataset prediction
            adata_obs = adata.obs
            #adata_obs["cell.type"].replace(['Oligo', 'NF Oligo'], 'Oligodendrocytes', inplace=True)
            #adata_obs["cell.type"].replace(['OPC'], 'Polydendrocytes', inplace=True)
            adata_obs["cell.type"].replace(['Oligo'], 'Oligodendrocytes', inplace=True)
            adata_obs["cell.type"].replace(['OPC', 'NF Oligo'], 'Polydendrocytes', inplace=True)
            adata_obs["cell.type"].replace(['Astro'], 'Astrocytes', inplace=True)
            adata_obs["cell.type"].replace(['Excitatory'], 'Neuron', inplace=True)
            adata_obs["cell.type"].replace(['Inhibitory'], 'Interneuron', inplace=True)
            adata_obs["cell.type"].replace(['Endo'], 'Endothelial', inplace=True)
            ## Microglia stays as microglia
            adata.obs = adata_obs

    elif 1 == celltype_gran:
        adata.obs.rename(columns={"L2_clusters": "cell.type"}, inplace=True)

    ## subset individuals
    if ind is not None:
        ind_cells = adata.obs[adata.obs['Sample'].isin(ind.split('_'))].index
        adata = adata[ind_cells]

    if devstage is not None:
        dev_cells = adata.obs[adata.obs["DevStage"] == devstage].index
        adata = adata[dev_cells]

    if treatment is not None:
        treat_cells = adata.obs[adata.obs["treatment"] == treatment].index
        adata = adata[treat_cells]

    return adata

