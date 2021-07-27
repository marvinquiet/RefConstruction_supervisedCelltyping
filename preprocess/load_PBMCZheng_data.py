import os
import anndata
import numpy as np
import pandas as pd

def combine_celltypes(data_dir):
    '''Merge the downloaded FACS-sorted data into one
    '''
    celltypes = ["b_cells", "cd14_monocytes", "cd34", "cd56_nk",
            "cd4_t_helper", "naive_t", "memory_t", "regulatory_t", ## CD4+ T cells
            "cytotoxic_t", "naive_cytotoxic"] ## CD8+ T cells
    adata_list = []
    for celltype in celltypes:
        celltype_dir = data_dir+os.sep+'PBMC_Zheng_FACS'+os.sep+celltype
        adata = anndata.read_mtx(celltype_dir+os.sep+'matrix.mtx').T
        ## load genes
        genes = pd.read_csv(celltype_dir+os.sep+'genes.tsv',
                header=None, sep='\t')
        adata.var['gene_symbols'] = genes[1].values
        adata.var_names = adata.var['gene_symbols']
        adata.var_names_make_unique(join="-")
        ## load cells
        cells = pd.read_csv(celltype_dir+os.sep+'barcodes.tsv',
                header=None, sep='\t')
        adata.obs['barcode'] = cells[0].values
        adata.obs_names = cells[0]
        adata.obs_names_make_unique(join="-")
        ## append adata
        adata_list.append(adata)

    final_adata = anndata.AnnData.concatenate(*adata_list,
            join='inner', batch_key="cell.type", batch_categories=celltypes) #inner
    final_adata.var.index.name = None
    final_adata.obs.index.name = None
    final_adata.write(data_dir+os.sep+'PBMC_Zheng_FACS/FACS_adata.h5ad')
    return final_adata

def load_PBMCZheng_data(data_dir, curate=False):
    '''Load PBMC Zheng FACS-sorted data
    '''
    adata = anndata.read_h5ad(data_dir+os.sep+'PBMC_Zheng_FACS/FACS_adata.h5ad')
    if curate:
        adata = adata[adata.obs["cell.type"] != "cd34"] ## remove cd34 cells
        adata_obs = adata.obs
        adata_obs["cell.type"].replace(['b_cells'], 'B cells', inplace=True)
        adata_obs["cell.type"].replace(['cd14_monocytes'], 'CD14+ Monocytes', inplace=True)
        adata_obs["cell.type"].replace(['cd56_nk'], 'NK cells', inplace=True)
        adata_obs["cell.type"].replace(['cd4_t_helper', 'naive_t', 'memory_t', 'regulatory_t'], 'CD4 T cells', inplace=True)
        adata_obs["cell.type"].replace(['cytotoxic_t', 'naive_cytotoxic'], 'CD8 T cells', inplace=True)
        adata.obs = adata_obs
    return adata

