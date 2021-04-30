import os
import anndata
import numpy as np
import pandas as pd

def write_brain_adata(data_dir, region="FC"):
    '''Loading data from different brain region and make it an anndata, store it
    @region: FC/HC

    FC: 29463 genes * 194027 cells -> 71,445 cells
    HC: 27953 genes * 134430 cells -> 53,204 cells

    Note: the mtx and genes/barcodes come from RDS file using writeMM and write
    Because the file is too large, I generated h5ad and deleted original data

    data <- readRDS("Hippocampus.RDS")
    writeMM(data, "mouse_HC.mtx")
    write(rownames(data), "mouse_HC_genes.tsv")
    write(colnames(data), "mouse_HC_barcodes.tsv")
    '''
    ## load data as anndata
    adata = anndata.read_mtx(data_dir+os.sep+'Mousebrain/mouse_'+region+'.mtx').T

    ## load cells and genes
    genes = pd.read_csv(data_dir+os.sep+'Mousebrain/mouse_'+region+'_genes.tsv',
            header=None, sep='\t')
    adata.var['gene_symbols'] = genes[0].values
    adata.var_names = adata.var['gene_symbols']
    adata.var_names_make_unique(join="-")

    cells = pd.read_csv(data_dir+os.sep+'Mousebrain/mouse_'+region+'_barcodes.tsv',
            header=None, sep='\t')
    adata.obs['barcode'] = cells[0].values
    adata.obs_names = cells[0]
    adata.obs_names_make_unique(join="-")

    ## load metadata information
    df = pd.read_csv(data_dir+os.sep+"Mousebrain/Mousebrain_metadata.csv", index_col=0)
    df = df[df["mouse_celltypes"] != "unknown"] # remove unknown cell types
    common_barcodes = set(df["barcodes"]).intersection(set(adata.obs["barcode"])) # 53,204 cells
    adata = adata[list(common_barcodes)]

    adata.obs = adata.obs.merge(df, left_on="barcode", right_on="barcodes")
    adata.obs.index = adata.obs["barcode"]
    adata.obs.index.name = None
    adata.var.index.name=None

    adata.write(data_dir+os.sep+'Mousebrain/'+region+'_adata.h5ad')

def load_brain_adata(data_dir, region="FC", ind="P60HippoRep1", 
        celltype=None, celltype_gran=0, curate=False):
    '''Return certain individual in the specific region
    
    @region: different brain regions
    @ind: individuals to return; if None, return all
        - FC: P60FCAldh1l1Rep1/P60FCCx3cr1Rep1/P60FCRep1,2,3,4,6
        - HC: P60HippoRep1-6
        can use _ as separator to input multiple
    @celltype_gran: granularity of celltype, 0: major cell types, 1: sub-celltypes
    '''
    adata = anndata.read_h5ad(data_dir+os.sep+'Mousebrain/'+region+'_adata.h5ad')
    if 0 == celltype_gran:
        adata.obs.rename(columns={"mouse_celltypes": "cell.type"}, inplace=True) # change cell types column name

        ## curate cell types for cross-region prediction
        if curate:
            if "FC" == region:
                adata_obs = adata.obs
                adata_obs["cell.type"].replace(['Interneuron_CGE', 'Interneuron_MGE'], 'Interneuron', inplace=True)
                adata_obs["cell.type"].replace(['Neuron_Claustrum', 'Neuron_L2/3', 
                    'Neuron_L5', 'Neuron_L5b', 'Neuron_L6'], 'Neuron', inplace=True)
                adata.obs = adata_obs
            elif "HC" == region:
                adata_obs = adata.obs
                adata_new_obs = adata_obs[(adata_obs['cell.type'] != "Choroid_Plexus") & (adata_obs['cell.type'] != "Ependyma") &
                        (adata_obs['cell.type'] != "Neurongenesis_Mitosis") ]
                adata = adata[adata_new_obs.index]
 
    if 1 == celltype_gran:
        adata.obs["cell.type"] = adata.obs["cluster"].astype(str)+'.'+adata.obs["mouse_celltypes"].astype(str)

    if ind is not None: 
        ind_cells = adata.obs[adata.obs["ind"].isin(ind.split('_'))].index
        adata = adata[ind_cells]

    return adata
