import os
import anndata
import numpy as np
import pandas as pd

def load_mouseprotocol_adata(data_dir, exp=None, protocol=None, curate=False):
    '''load Mouse cortex data from different protocols

    After comparing between count.umis.txt and count.reads.txt, although the
    same dimension, there are entries different with each other

    - Extract out plate-based cells from count.reads.txt
    - Extract out droplet-based cells from counts.umi.txt
    - Concatenate to anndata

    @exp: cortex1, cortex2
    @protocol: plate-based (Smart-seq2), droplet-based(DroNc-seq, sci-RNA-seq,
    10x Chromium)
    @curate: whether to curate for cell types
    '''

    plate_protocols = ["Smart-seq2"]

    metadata_df = pd.read_csv(data_dir+os.sep+"Mousecortex_protocols/metadata.txt",
            header=0, sep="\t")
    metadata_df = metadata_df[metadata_df["CellType"] != "Unassigned"]  ## remove unassigned

    ## they used the same cells and genes indicator
    cells = pd.read_csv(data_dir+os.sep+"Mousecortex_protocols/cell.names.new.txt", header=None)
    genes = pd.read_csv(data_dir+os.sep+"Mousecortex_protocols/genes.counts.txt", header=None)

    ## plate-based data
    read_adata = anndata.read_mtx(data_dir+os.sep+"Mousecortex_protocols/count.reads.txt").T
    read_adata.var['gene_symbols'] = [x.split('_')[1] for x in genes[0].values]
    read_adata.var_names = read_adata.var['gene_symbols']
    read_adata.var_names_make_unique(join="-") # make unique
    read_adata.var_names.name = None

    read_adata.obs['barcode'] = cells[0].values
    read_adata.obs_names = read_adata.obs['barcode']
    read_adata.obs_names_make_unique(join="-") ## make unique
    read_adata.obs_names.name = None

    plate_metadata = metadata_df[metadata_df['Method'].isin(plate_protocols)]
    common_cells = set(plate_metadata['NAME']).intersection(set(read_adata.obs_names))
    common_cells = list(common_cells)
    read_adata = read_adata[common_cells]

    obs_df = read_adata.obs.merge(plate_metadata, how='left', 
            left_index=True, right_on='NAME')
    obs_df.index = obs_df['barcode'].values
    read_adata.obs = obs_df

    ## umi-based data
    umi_adata = anndata.read_mtx(data_dir+os.sep+"Mousecortex_protocols/count.umis.txt").T
    umi_adata.var['gene_symbols'] = [x.split('_')[1] for x in genes[0].values]
    umi_adata.var_names = umi_adata.var['gene_symbols']
    umi_adata.var_names_make_unique(join="-") # make unique
    umi_adata.var_names.name = None

    umi_adata.obs['barcode'] = cells[0].values
    umi_adata.obs_names = umi_adata.obs['barcode']
    umi_adata.obs_names_make_unique(join="-") ## make unique
    umi_adata.obs_names.name = None

    droplet_metadata = metadata_df[~metadata_df['Method'].isin(plate_protocols)]
    common_cells = set(droplet_metadata['NAME']).intersection(set(umi_adata.obs_names))
    common_cells = list(common_cells)
    umi_adata = umi_adata[common_cells] 

    obs_df = umi_adata.obs.merge(droplet_metadata, how='left', 
            left_index=True, right_on='NAME')
    obs_df.index = obs_df['barcode'].values
    umi_adata.obs = obs_df

    ## concatenate adata together
    adata = read_adata.concatenate(umi_adata, batch_key="protocol_type", 
            batch_categories=['plate', 'droplet'])
    adata.obs.rename(columns={'CellType': 'cell.type'}, inplace=True)

    adata_obs = adata.obs
    adata_obs['Method'].replace(['10x Chromium'], '10x', inplace=True)
    adata.obs = adata_obs

    if exp is not None:
        exp_cells = adata.obs[adata.obs['Experiment'] == exp].index
        adata = adata[exp_cells]

    if protocol is not None:
        proc_cells = adata.obs[adata.obs['Method'] == protocol].index
        adata = adata[proc_cells]

    if curate:
        adata_obs = adata.obs
        adata_obs["cell.type"].replace(['Astrocyte'], 'Astrocytes', inplace=True)
        adata_obs["cell.type"].replace(['Excitatory neuron'], 'Neuron', inplace=True)
        adata_obs["cell.type"].replace(['Inhibitory neuron'], 'Interneuron', inplace=True)
        adata_obs["cell.type"].replace(['Oligodendrocyte'], 'Oligodendrocytes', inplace=True)
        adata_obs["cell.type"].replace(['OPC'], 'Polydendrocytes', inplace=True)
        #adata_obs["cell.type"].replace(['Pericyte'], 'Mural', inplace=True) ## seems not the same cell types
        ## Endothelial and Microglia as it is
        adata.obs = adata_obs

    return adata
