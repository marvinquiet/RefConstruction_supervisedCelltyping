import os
import anndata
import numpy as np
import pandas as pd

def filter_mouseAllen_10X_data(data_dir):
    '''Filter out frontal cortex data (ACA, PL, ILA, ORB (the ORB cannot be separated, therefore included)) from mouse Allen brain 10X dataset
    https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-whole-cortex-and-hippocampus-10x
    '''
    metadata = pd.read_csv(data_dir+os.sep+'AllenBrain/sub_metadata.csv', header=0, sep=',')
    idx = ((metadata['region_label'] == 'ACA') | (metadata['region_label'] == 'PL;ILA;ORB')) & (pd.notna(metadata['cell_type_alias_label'])) ## extract out cells from ACA, PL, ILA, ORB -> as frontal cortex
    reg_metadata = metadata[idx]
    reg_metadata['cell.type'] = [x.split('_')[1] for x in reg_metadata['cell_type_alias_label']]
    reg_metadata.to_csv(data_dir+os.sep+'AllenBrain/filtered_metadata.csv')
    reg_cells = reg_metadata['sample_name'].tolist() ## 106,892 cells

    files = [f for f in os.listdir(data_dir+os.sep+'AllenBrain') if 'smaller_matrix' in f]
    fout = open(data_dir+os.sep+'AllenBrain'+os.sep+'filtered_matrix.csv', 'w')
    ## write header
    data_header_file = open(data_dir+os.sep+'AllenBrain/matrix_header.csv', 'r')
    data_header = data_header_file.readline().strip()
    data_header_file.close()
    fout.write(data_header+'\n')
    for f in files:
        print("Reading {}...".format(f))
        fpath = data_dir+os.sep+'AllenBrain'+os.sep+f
        fop = open(fpath, 'r')
        for line in fop:
            sample_id = line.strip().split(',')[0]
            if sample_id in reg_cells:
                fout.write(line)
        fop.close()
    fout.close()

def load_mouseAllen_10X_data(data_dir, ind=None, celltype_gran=0, curate=False):
    '''Load mouse Allen cortex data from 10X reads

    @ind: individual sample -> 49 mice originally, 3 remained after filtering on male 
    @celltype_gran: major cell types (0) or sub-cell types (1)
    @curate: whether to curate to the same name as mouse FC, mouse pFC, mouse protocol
    '''
    adata = anndata.read_csv(data_dir+os.sep+'AllenBrain/filtered_matrix.csv')
    adata.var["gene_symboles"] = adata.var.index
    adata.var_names_make_unique(join="-")
    adata.obs['barcode'] = adata.obs.index
    adata.obs_names_make_unique(join="-")

    ## load metadata
    meta_df = pd.read_csv(data_dir+os.sep+'AllenBrain/filtered_metadata.csv', index_col=0)
    meta_df['external_donor_name_label'] = meta_df['external_donor_name_label'].astype(str)
    sample_info = pd.read_csv(data_dir+os.sep+'AllenBrain/sample_info.csv', index_col=0) ## use to filter gender
    merged_meta_df = meta_df.merge(sample_info, how='left', left_on="external_donor_name_label", right_on="external_donor_name")
    adata.obs = adata.obs.merge(merged_meta_df, left_index=True, right_on="sample_name", how='left')  ## merge anndata with metadata
    adata.obs.index = adata.obs["barcode"]
    adata.obs.index.name = None
    adata = adata[adata.obs['sex_label'] == 'M']  ## filter males out

    if 0 == celltype_gran:
        if curate:
            adata_obs = adata.obs
            adata_obs['cell.type'].replace(['L2 IT HATA', 'L2 IT RSP-ACA', 'L2/3 IT APr', 
                'L2/3 IT CTX', 'L2/3 IT PPP', 'L3 RSP-ACA', 'L4/5 IT CTX', 'L5 IT CTX', 
                'L5 NP CTX', 'L5 PT CTX', 'L5 PT RSP-ACA', 'L5/6 IT CTX', 'L6 CT CTX',
                'L6 Car3', 'L6 IT CTX', 'L6 NP CT CTX', 'L6b CTX','NP PPP', 'NP SUB', 
                'CA2-IG-FC', 'L6b RHP', 'L6 IT RHP', 'L2/3 IT ENTl', 'L2 IT RSPv', 
                'L5 IT TPE-ENT', 'L2/3 IT TPE', 'L2 IT ProS'], 'Neuron', inplace=True) ## 27 types
            adata_obs['cell.type'].replace(['Sncg', 'Pvalb', 'Pvalb Vipr2', 'Sst',
                'Sst Chodl', 'Vip', 'Lamp5', 'Lamp5 Lhx6', 'Pax6', 'Ndnf HPF'], 'Interneuron', inplace=True) ## 10 types
            adata_obs['cell.type'].replace(['Astro'], 'Astrocytes', inplace=True)
            adata_obs['cell.type'].replace(['Endo'], 'Endothelial', inplace=True)
            adata_obs['cell.type'].replace(['Micro'], 'Microglia', inplace=True)
            adata_obs['cell.type'].replace(['OPC'], 'Polydendrocytes', inplace=True)
            adata_obs['cell.type'].replace(['Oligo'], 'Oligodendrocytes', inplace=True)
            adata_obs['cell.type'].replace(['Peri'], 'Pericytes', inplace=True)
            adata.obs = adata_obs
            ## 'CR', 'PVM', 'SMC', 'VLMC' celltype are out of all other database, therefore not included
            remained_idx = ~adata.obs['cell.type'].isin(['CR', 'PVM', 'SMC', 'VLMC'])
            adata = adata[remained_idx]

    if 1 == celltype_gran: ## for sub-celltypes
        adata.obs.rename(columns={"cell.type": "major_celltype"}, inplace=True)
        adata.obs.rename(columns={"cell_type_alias_label": "cell.type"}, inplace=True)

    if ind is not None:
        ind_cells = adata.obs[adata.obs["external_donor_name_label"].isin(ind.split('_'))].index
        adata = adata[ind_cells]

    return adata

def filter_mouseAllen_SS_data(data_dir):
    '''Filter out frontal cortex data with Smart-seq2 sequenced
    Downloaded from: http://data.nemoarchive.org/biccn/grant/u19_zeng/zeng/transcriptome/scell/SSv4/mouse/processed/YaoHippo2020/

    Important Note: in our original thought, we wanted to separate frontal cortex region and select males as our data. However,
        it is impossible because even for the whole dataset, some cell types are lacked for most individuals due to the FACS-sorted
        and capture of cells for Smart-Seqv4.
    -> Therefore, for SSv4, we will only filter out males and then do cross-validation without filtering out specific brain regions.
    '''
    import h5py
    import scipy.sparse as ss

    fname = data_dir+os.sep+'AllenBrain/smrt.h5'
    h5f = h5py.File(fname,'r')
    exons = h5f['/data/exon/'] ## extract exons reads of genes
    x = exons['x']
    i = exons['i']
    p = exons['p']
    dims = exons['dims']
    sparse_mat = ss.csc_matrix((x[0:x.len()],
                               i[0:i.len()],
                               p[0:p.len()]),
                               shape = (dims[0],dims[1]))
    ## get gene names and sample names
    gene_names = h5f['/gene_names']
    sample_names = h5f['/sample_names']
    genes = [x.decode('utf-8') for x in gene_names]
    samples = [x.decode('utf-8') for x in sample_names]

    ## create anndata object
    adata = anndata.AnnData(sparse_mat)
    adata.var['gene_symbols'] = genes
    adata.var_names = genes
    adata.var_names_make_unique(join="-")
    adata.obs['barcode'] = samples
    adata.obs_names = samples
    adata.obs_names_make_unique(join="-")

    metadata = pd.read_csv(data_dir+os.sep+"AllenBrain/CTX_Hip_anno_SSv4.csv", index_col=0)
    metadata['donor_label'] = metadata['donor_label'].astype(str)
    adata.obs = adata.obs.merge(metadata, how="left", left_on="barcode", right_on="exp_component_name")
    adata.obs.index = adata.obs["barcode"]
    adata.obs.index.name = None

    ## filter out brain region 
    #idx = ((adata.obs['region_label'] == 'ACA') | (adata.obs['region_label'] == 'PL-ILA')) & \
    #        (adata.obs['sex_label'] == 'M') & (pd.notna(adata.obs['label'])) ## extract out cells from ACA, PL, ILA, ORB -> as frontal cortex
    idx = (adata.obs['sex_label'] == 'M') & (pd.notna(adata.obs['label'])) ## extract out cells from ACA, PL, ILA, ORB -> as frontal cortex
    filtered_adata = adata[idx]

    ## write out metadata and matrix
    filtered_adata.obs.to_csv(data_dir+os.sep+'AllenBrain/filtered_SS_metadata.csv')
    mat_df = pd.DataFrame(filtered_adata.X.toarray().astype(int), index=filtered_adata.obs_names,
            columns=filtered_adata.var_names) ## shape: (45626, 45768)
    mat_df.to_csv(data_dir+os.sep+'AllenBrain/filtered_SS_matrix.csv')

def load_mouseAllen_SS_data(data_dir, ind=None, celltype_gran=0, curate=False):
    '''Load mouse Allen cortex data from Smart-seq exon reads

    @ind: individual sample -> 529 mice originally, 19 mice after filtering
    @celltype_gran: major cell types (0) or sub-cell types (1)
    @curate: whether to curate to the same name as mouse FC, mouse pFC, mouse protocol
    '''
    adata = anndata.read_csv(data_dir+os.sep+'AllenBrain/filtered_SS_matrix.csv')
    adata.var["gene_symboles"] = adata.var.index
    adata.var_names_make_unique(join="-")
    adata.obs['barcode'] = adata.obs.index
    adata.obs_names_make_unique(join="-")

    ## load metadata
    meta_df = pd.read_csv(data_dir+os.sep+'AllenBrain/filtered_SS_metadata.csv', index_col=0)
    adata.obs = adata.obs.merge(meta_df, left_on='barcode', right_on='barcode', how='left')  ## merge anndata with metadata
    adata.obs.index = adata.obs["barcode"]
    adata.obs.index.name = None
    adata.obs = adata.obs[['barcode', 'donor_label', 'region_label', 'joint_region_label', 'label', 'cluster_label', 'ss_cluster_label']]

    if 0 == celltype_gran:
        adata.obs.rename(columns={'label':'cell.type'}, inplace=True)
        if curate:
            adata_obs = adata.obs
            adata_obs['cell.type'].replace(['L2/3 IT CTX', 'L2/3 IT RSP', 'L2/3 IT RSPv',
                'L2 IT ENTl', 'L4/5 IT CTX', 'L4 RSP-ACA', 'L5/6 IT CTX', 'L5 IT CTX',
                'L5 IT TPE-ENT', 'L5 NP CT CTX ', 'L5 NP CTX', 'L5 PT CTX', 'L6b CTX', 
                'L6 CT CTX', 'L6 IT CTX', 'CA1', 'CA1-ProS', 'CA2-IG-FC', 'CA3',
                'CT SUB', 'Car3', 'DG', 'L2 IT ENTm', 'L2 IT PAR', 'L2/3 IT APr',
                'L2/3 IT ENTl', 'L2/3 IT HATA', 'L2/3 IT PPP', 'L3 IT ENTl', 'L3 IT ENTm',
                'L5 PPP', 'L6 IT ENTl', 'L6b/CT ENT', 'NP PPP', 'NP SUB', 'SUB', 'SUB-ProS'], 'Neuron', inplace=True)
            adata_obs['cell.type'].replace(['Lamp5', 'Lamp5 Lhx6', 'Pax6', 'Pvalb',
                'Sncg', 'Sst', 'Sst Chodl', 'Vip'], 'Interneuron', inplace=True)
            adata_obs['cell.type'].replace(['Astro'], 'Astrocytes', inplace=True)
            adata_obs['cell.type'].replace(['Endo'], 'Endothelial', inplace=True)
            adata_obs['cell.type'].replace(['Micro-PVM'], 'Microglia', inplace=True)
            adata_obs['cell.type'].replace(['Oligo'], 'Oligodendrocytes', inplace=True)
            adata_obs['cell.type'].replace(['SMC-Peri'], 'Pericytes', inplace=True)
            adata.obs = adata_obs
            ## 'CR', 'Meis2', 'Meis2 HPF', 'Ndnf HPF', 'VLMC' celltype are out of all other database, therefore not included
            ## The polydendrocytes (OPCs) are missing
            remained_idx = ~adata.obs['cell.type'].isin(['CR', 'Meis2', 'Meis2 HPF', 'Ndnf HPF', 'VLMC'])
            adata = adata[remained_idx]

    if 1 == celltype_gran: ## for sub-celltypes
        adata.obs.rename(columns={"cluster_label": "cell.type"}, inplace=True)

    if ind is not None:
        ind_cells = adata.obs[adata.obs["donor_label"].isin(ind.split('_'))].index
        adata = adata[ind_cells]

    return adata
