import os
import anndata
import numpy as np
import pandas as pd

## --- curate PBMC demultiplex cell types to two granularity
def curate_PBMC_demulx_celltypes(adata, celltype_gran=1, 
        celltype_label="cell.type"):
    '''Curate PBMC demulx cell types into two-starge prediction

    @adata: anndata object
    @celltype_gran: default 1 because the given cell.type are sub-cell types;
        if 0, then curate to major cell types
    @celltype_label: column indicator for cell type in adata.obs
    '''

    celltype_categories = {'B cells': ['B cells'],
            'T cells': ['CD4 T cells', 'CD8 T cells'],
            'NK cells': ['NK cells'],
            'Dendritic cells': ['Dendritic cells'],
            'Monocytes': ['CD14+ Monocytes', 'FCGR3A+ Monocytes'],
            'Megakaryocytes': ['Megakaryocytes']
            }
    major_celltypes = []
    for celltype in adata.obs["cell.type"].tolist():
        flag = False
        for major_celltype, sub_celltypes in celltype_categories.items():
            if celltype in sub_celltypes:
                major_celltypes.append(major_celltype)

                ## if found
                flag = True
                break
        if False == flag:
            major_celltypes.append('unknown')
    adata.obs['majortypes'] = major_celltypes

    if 0 == celltype_gran:
        adata.obs.rename(columns={celltype_label: "subtypes", 
                                    "majortypes": celltype_label}, 
                        inplace=True)
    return adata

# --------------
# load PBMC liger data functions
# --------------
def load_PBMC_liger_data(data_dir):
    '''Loading PBMC liger 6k cells data
    '''
    ## load data as anndata
    adata = anndata.read_mtx(data_dir+os.sep+'PBMC_demuxlet/matrix.mtx').T

    ## load cells and genes
    genes = pd.read_csv(data_dir+os.sep+'PBMC_demuxlet/genes.tsv', 
            header=None, sep='\t')
    adata.var['gene_symbols'] = genes[0].values
    adata.var_names = adata.var['gene_symbols']
    # Make sure the gene names are unique
    adata.var_names_make_unique(join="-")

    cells = pd.read_csv(data_dir+os.sep+'PBMC_demuxlet/barcodes.tsv', 
            header=None, sep='\t')
    adata.obs['barcode'] = cells[0].values
    adata.obs_names = cells[0]
    # Make sure the cell names are unique
    adata.obs_names_make_unique(join="-")
    return adata

# --------------
# load PBMC batch2 data functions
# --------------
def load_PBMC_batch2_data(data_dir, condition=None, ind=None):
    '''Loading PBMC batch2 data

    @condition: ctrl/stim
    @ind: 101 1015 1016 1039 107 1244 1256 1488
    '''
    ## load genes
    genes = pd.read_csv(data_dir+os.sep+'PBMC_demuxlet/GSE96583_batch2.genes.tsv',
                        header=None, sep='\t')
 
    ## load control data
    ctrl_adata = anndata.read_mtx(data_dir+os.sep+'PBMC_demuxlet/GSM2560248_2.1.mtx').T
    ctrl_adata.var['gene_symbols'] = genes[1].values
    ctrl_adata.var_names = ctrl_adata.var['gene_symbols']
    # Make sure the gene names are unique
    ctrl_adata.var_names_make_unique(join="-")

    cells = pd.read_csv(data_dir+os.sep+'PBMC_demuxlet/GSM2560248_barcodes.tsv',
                        header=None, sep='\t')
    ctrl_adata.obs['barcode'] = 'ctrl' + cells[0].values
    ctrl_adata.obs_names = cells[0]
    # Make sure the cell names are unique
    ctrl_adata.obs_names_make_unique(join="-")

    ## load stim data
    stim_adata = anndata.read_mtx(data_dir+os.sep+'PBMC_demuxlet/GSM2560249_2.2.mtx').T
    stim_adata.var['gene_symbols'] = genes[1].values
    stim_adata.var_names = stim_adata.var['gene_symbols']
    # Make sure the gene names are unique
    stim_adata.var_names_make_unique(join="-")

    cells = pd.read_csv(data_dir+os.sep+'PBMC_demuxlet/GSM2560249_barcodes.tsv',
                        header=None, sep='\t')
    stim_adata.obs['barcode'] = 'stim' + cells[0].values
    stim_adata.obs_names = cells[0]
    # Make sure the cell names are unique
    stim_adata.obs_names_make_unique(join="-")

    ## combine control/stimulated data together
    adata = ctrl_adata.concatenate(stim_adata, 
            batch_key="condition", batch_categories=['control', 'stimulated'])
    adata.obs.index = adata.obs["barcode"]
    adata.obs.index.name = None

    ## load meta data information
    PBMC_batch2_df = load_PBMC_batch2_df(data_dir)
    common_barcodes = set(PBMC_batch2_df['barcode']).intersection(set(adata.obs['barcode']))
    adata = adata[list(common_barcodes)]

    adata.obs = adata.obs.merge(PBMC_batch2_df, left_on="barcode", right_on="barcode")
    adata.obs.index = adata.obs["barcode"]
    adata.obs.index.name = None
    adata.var.index.name = None

    adata.obs.rename(columns={'cell': 'cell.type'}, inplace=True)

    if condition is not None:
        cond_cells = adata.obs[adata.obs["condition"] == condition].index
        adata = adata[cond_cells]

    if ind is not None:
        ind_cells = adata.obs[adata.obs["ind"].isin(ind.split('_'))].index
        adata = adata[ind_cells]
    return adata


def load_PBMC_batch2_df(data_dir):
    '''load PBMC metadata information
    filter the cells with singlet only and add barcode column as unique identifier
    '''

    ## load meta data information
    PBMC_df = pd.read_csv(data_dir+os.sep+"PBMC_demuxlet/GSE96583_batch2.total.tsne.df.tsv", 
            header=0, sep="\t")
    ## find singlet, 24,673 cells
    PBMC_df = PBMC_df[(PBMC_df["multiplets"] == "singlet") & (~PBMC_df["cell"].isnull())]
    PBMC_df['barcode'] = PBMC_df['stim'] + PBMC_df.index
    return PBMC_df

# --------------
# load PBMC batch1 data functions
# --------------
def load_PBMC_batch1_data(data_dir, batch=None, ind=None):
    '''Loading PBMC batch1 S1 dataset from W1-3

    @batch: A/B/C
    @ind: 1043 1079 1154 1249 1493 1511 1598 1085
    '''
    ## load batch1 genes
    genes = pd.read_csv(data_dir+os.sep+'PBMC_demuxlet/GSE96583_batch1.genes.tsv',
                        header=None, sep='\t')
 
    ## load matrix A data
    A_adata = anndata.read_mtx(data_dir+os.sep+'PBMC_demuxlet/GSM2560245_A.mat').T # 3639 inviduals
    ## load cells
    A_adata.var['gene_symbols'] = genes[1].values
    A_adata.var_names = A_adata.var['gene_symbols']
    A_adata.var_names_make_unique(join="-") # make unique

    A_cells = pd.read_csv(data_dir+os.sep+'PBMC_demuxlet/GSM2560245_barcodes.tsv',
                        header=None, sep='\t')
    A_adata.obs['barcode'] = 'A_' + A_cells[0].values
    A_adata.obs_names = A_cells[0]
    A_adata.obs_names_make_unique(join="-") # make unique

    ## load matrix B data
    B_adata = anndata.read_mtx(data_dir+os.sep+'PBMC_demuxlet/GSM2560246_B.mat').T # 4246 inviduals
    ## load cells
    B_adata.var['gene_symbols'] = genes[1].values
    B_adata.var_names = B_adata.var['gene_symbols']
    B_adata.var_names_make_unique(join="-") # make unique

    B_cells = pd.read_csv(data_dir+os.sep+'PBMC_demuxlet/GSM2560246_barcodes.tsv',
                        header=None, sep='\t')
    B_adata.obs['barcode'] = 'B_' + B_cells[0].values
    B_adata.obs_names = B_cells[0]
    B_adata.obs_names_make_unique(join="-") # make unique

    ## load matrix C data
    C_adata = anndata.read_mtx(data_dir+os.sep+'PBMC_demuxlet/GSM2560247_C.mat').T # 6145 inviduals
    ## load cells
    C_adata.var['gene_symbols'] = genes[1].values
    C_adata.var_names = C_adata.var['gene_symbols']
    C_adata.var_names_make_unique(join="-") # make unique

    C_cells = pd.read_csv(data_dir+os.sep+'PBMC_demuxlet/GSM2560247_barcodes.tsv',
                        header=None, sep='\t')
    C_adata.obs['barcode'] = 'C_' + C_cells[0].values
    C_adata.obs_names = C_cells[0]
    C_adata.obs_names_make_unique(join="-") # make unique

    # combine data together
    adata = A_adata.concatenate(B_adata, C_adata,
            batch_key="batch", batch_categories=['A', 'B', 'C'])
    adata.obs.index = adata.obs["barcode"]
    adata.obs.index.name = None

    ## load meta data information
    PBMC_batch1_df = load_PBMC_batch1_df(data_dir)
    common_barcodes = set(PBMC_batch1_df['barcode']).intersection(set(adata.obs['barcode']))
    adata = adata[list(common_barcodes)]

    adata.obs = adata.obs.merge(PBMC_batch1_df, on=["barcode", "batch"], how="left")
    adata.obs.index = adata.obs["barcode"]
    adata.obs.index.name = None
    adata.var.index.name = None

    if batch is not None:
        batch_cells = adata.obs[adata.obs['batch'] == batch].index
        adata = adata[batch_cells]

    if ind is not None:
        ind_list = [int(x) for x in ind.split('_')]
        ind_cells = adata.obs[adata.obs["ind"].isin(ind_list)].index
        adata = adata[ind_cells]
    return adata

def load_PBMC_batch1_df(data_dir):
    '''load PBMC metadata information
    filter the cells with singlet only and add barcode column as unique identifier
    '''

    ## load meta data information
    PBMC_df = pd.read_csv(data_dir+os.sep+"PBMC_demuxlet/GSE96583_batch1.total.tsne.df.tsv", 
            header=0, sep="\t")
    ## find singlet, 11,815 cells
    PBMC_df = PBMC_df[(PBMC_df["multiplets"] == "singlet") & (~PBMC_df["cell.type"].isnull())]
    PBMC_df['barcode'] = PBMC_df['batch'] + '_' + PBMC_df.index
    return PBMC_df
