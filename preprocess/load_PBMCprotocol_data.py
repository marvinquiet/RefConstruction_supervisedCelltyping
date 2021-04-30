import os
import anndata
import numpy as np
import pandas as pd

# --------------
# load PBMC protocols function
# --------------
def load_PBMC_protocols_data(data_dir):
    '''load PBMC from different protocols

    After comparing between cells.read.new.txt and cells.umi.new.txt, the cell names
    remain as cells.umi.new.txt. Therefore, compare the count matrix and check reads.
    Turn out to be umi and read count are not the same..

    - Extract out plate-based cells from counts.read.txt
    - Extract out droplet-based cells from counts.umi.txt
    - Concatenate to anndata

    @ protocol: plate-based (Smart-Seq2/CEL-Seq2), droplet-based (10X v2, 10X v3,
    Drop-seq, Seq-Well, inDrops)
    '''
    plate_protocols = ["CEL-Seq2", "Smart-seq2"]

    metadata_df = pd.read_csv(data_dir+os.sep+"PBMC_protocols/metadata.txt", 
            header=0, sep="\t")

    ## plate-based data
    read_adata = anndata.read_mtx(data_dir+os.sep+"PBMC_protocols/counts.read.txt").T
    read_cells = pd.read_csv(data_dir+os.sep+"PBMC_protocols/cells.read.new.txt",
            header=None)
    read_genes = pd.read_csv(data_dir+os.sep+"PBMC_protocols/genes.read.txt",
            header=None)

    read_adata.var['gene_symbols'] = [x.split('_')[1] for x in read_genes[0].values]
    read_adata.var_names = read_adata.var['gene_symbols']
    read_adata.var_names_make_unique(join="-") # make unique
    read_adata.var_names.name = None

    read_adata.obs['barcode'] = read_cells[0].values
    read_adata.obs_names = read_adata.obs['barcode']
    read_adata.obs_names_make_unique(join="-") ## make unique
    read_adata.obs_names.name = None

    plate_metadata = metadata_df[metadata_df['Method'].isin(plate_protocols)]
    common_cells = set(plate_metadata['NAME']).intersection(set(read_adata.obs_names))
    common_cells = list(common_cells)
    read_adata = read_adata[common_cells] # 1052 cells

    obs_df = read_adata.obs.merge(plate_metadata, how='left', 
            left_index=True, right_on='NAME')
    obs_df.index = obs_df['barcode'].values
    read_adata.obs = obs_df

    ## umi-based data
    umi_adata = anndata.read_mtx(data_dir+os.sep+"PBMC_protocols/counts.umi.txt").T
    umi_cells = pd.read_csv(data_dir+os.sep+"PBMC_protocols/cells.umi.new.txt",
            header=None)
    umi_genes = pd.read_csv(data_dir+os.sep+"PBMC_protocols/genes.umi.txt",
            header=None)

    umi_adata.var['gene_symbols'] = [x.split('_')[1] for x in umi_genes[0].values]
    umi_adata.var_names = umi_adata.var['gene_symbols']
    umi_adata.var_names_make_unique(join="-") # make unique
    umi_adata.var_names.name = None

    umi_adata.obs['barcode'] = umi_cells[0].values
    umi_adata.obs_names = umi_adata.obs['barcode']
    umi_adata.obs_names_make_unique(join="-") ## make unique
    umi_adata.obs_names.name = None

    droplet_metadata = metadata_df[~metadata_df['Method'].isin(plate_protocols)]
    common_cells = set(droplet_metadata['NAME']).intersection(set(umi_adata.obs_names))
    common_cells = list(common_cells)
    umi_adata = umi_adata[common_cells] # 29969 cells

    obs_df = umi_adata.obs.merge(droplet_metadata, how='left', 
            left_index=True, right_on='NAME')
    obs_df.index = obs_df['barcode'].values
    umi_adata.obs = obs_df

    ## concatenate adata together
    adata = read_adata.concatenate(umi_adata, batch_key="protocol_type", 
            batch_categories=['plate', 'droplet'])
    adata.obs.rename(columns={'CellType': 'cell.type'}, inplace=True)

    adata_obs = adata.obs
    adata_obs['Method'].replace(['10x Chromium (v2)', '10x Chromium (v2) A',
        '10x Chromium (v2) B'], '10x-v2', inplace=True)
    adata_obs['Method'].replace(['10x Chromium (v3)'], '10x-v3', inplace=True)
    adata.obs = adata_obs

    return adata 
