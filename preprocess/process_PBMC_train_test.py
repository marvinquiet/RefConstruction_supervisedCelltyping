import os, sys
import numpy as np
import pandas as pd
import scipy
import scanpy.api as sc
import anndata
import random

from preprocess.process_train_test_data import *

from preprocess import load_PBMC_data
from preprocess import load_PBMCprotocol_data


def process_batch1_ind(data_dir, result_dir, gene_no=1000, dr_seed=0,
        ind1="1154", ind2="1085", select_on="test", select_method="Seurat",
        purify=False, purify_method="distance", celltype_gran=1):
    ''' Process individual data of PBMC batch1

    @data_dir: where PBMC batch1 data stroes
    @result_dir: where to store PCA/tSNE/UMAP result
    @gene_no: number of HVGs called using Seurat
    @dr_seed: dimension reduction seed for reproducibility
    @celltype_gran: granularity of cell types, 0: major cell types, 1:sub-celltypes
        According to the dataset, the give cell types are sub-cell types
    '''
    load_ind, train_adata, test_adata = load_adata(result_dir)
    if load_ind:
        return train_adata, test_adata

    ## add batch info
    train_adata = load_PBMC_data.load_PBMC_batch1_data(data_dir, ind=ind1)
    test_adata = load_PBMC_data.load_PBMC_batch1_data(data_dir, ind=ind2)
    train_adata.obs.index.name = None
    test_adata.obs.index.name = None

    ## curate given sub-cell types to major cell types
    train_adata = load_PBMC_data.curate_PBMC_demulx_celltypes(train_adata, celltype_gran)
    test_adata = load_PBMC_data.curate_PBMC_demulx_celltypes(test_adata, celltype_gran)

    ## feature selection
    train_adata, test_adata = feature_selection_train_test(train_adata, test_adata,
            result_dir, gene_no, select_on, select_method)

    ## scale and analze
    train_adata, test_adata = scale_and_visualize(train_adata, test_adata, result_dir,
            purify=purify, purify_method=purify_method)
    save_adata(train_adata, test_adata, result_dir)

    return train_adata, test_adata

def process_batch1_multiinds(data_dir, result_dir, sample_ind, pred_ind,
        gene_no=1000, dr_seed=0, sample_seed=0, select_on="test",
        select_method="Seurat", celltype_gran=1, sample=False):
    '''Process multi individuals to predict one individual

    @sample_ind: downsample to a certain individual number
    @pred_ind: the predicted individual
    '''
    load_ind, train_adata, test_adata = load_adata(result_dir)
    if load_ind:
        return train_adata, test_adata

    batch1_inds = [1043, 1079, 1154, 1249, 1493, 1511, 1598, 1085]

    exclude_list = batch1_inds
    exclude_list.remove(int(pred_ind))

    train_adata = load_PBMC_data.load_PBMC_batch1_data(data_dir, ind=exclude_list)
    if sample: ## sample to average number
        avg_number = train_adata.shape[0] // len(exclude_list)

        print("=== Downsample to number:", avg_number)
        random.seed(sample_seed)
        sampled_cells = random.sample(list(train_adata.obs_names), k=avg_number)
        train_adata = train_adata[sampled_cells]

    test_adata = load_PBMC_data.load_PBMC_batch1_data(data_dir, ind=pred_ind)
    train_adata.obs.index.name = None
    test_adata.obs.index.name = None

    ## feature selection
    train_adata, test_adata = feature_selection_train_test(train_adata, test_adata,
            result_dir, gene_no, select_on, select_method)

    ## scale and analze
    train_adata, test_adata = scale_and_visualize(train_adata, test_adata, result_dir)
    save_adata(train_adata, test_adata, result_dir)

    return train_adata, test_adata


def process_batch1_ABC(data_dir, result_dir, gene_no=1000, dr_seed=0,
        batch1="A", batch2="B", select_on="test", select_method="Seurat",
        purify=False, purify_method="distance", celltype_gran=1):
    ''' Process two batches for GEDFN

    @data_dir: where PBMC batch1/batch2 data stroes
    @result_dir: where to store PCA/tSNE/UMAP result
    @gene_no: number of HVGs called using Seurat
    @dr_seed: dimension reduction seed for reproducibility
    @celltype_gran: default 1, given as sub-cell types
    '''
    load_ind, train_adata, test_adata = load_adata(result_dir)
    if load_ind:
        return train_adata, test_adata

    ## add batch info
    train_adata = load_PBMC_data.load_PBMC_batch1_data(data_dir, batch = batch1)
    test_adata = load_PBMC_data.load_PBMC_batch1_data(data_dir, batch = batch2)
    train_adata.obs.index.name = None
    test_adata.obs.index.name = None

    ## curate given sub-cell types to major cell types
    train_adata = load_PBMC_data.curate_PBMC_demulx_celltypes(train_adata, celltype_gran)
    test_adata = load_PBMC_data.curate_PBMC_demulx_celltypes(test_adata, celltype_gran)

    ## feature selection
    train_adata, test_adata = feature_selection_train_test(train_adata, test_adata,
            result_dir, gene_no, select_on, select_method)

    ## scale and analze
    train_adata, test_adata = scale_and_visualize(train_adata, test_adata, result_dir, 
            purify=purify, purify_method=purify_method)
    save_adata(train_adata, test_adata, result_dir)

    return train_adata, test_adata


def process_batch1_batchtoind(data_dir, result_dir, gene_no=1000, dr_seed=0,
        batch="A", ind="1085", select_on="test", select_method="Seurat", 
        sample=True, celltype_gran=1):
    '''Process batch1 A -> S5, as opposed to S1 -> S5

    @batch: batch A as training datasets
    @ind: ind S5 as predictor
    @sample: whether to downsample or not
    '''
    load_ind, train_adata, test_adata = load_adata(result_dir)
    if load_ind:
        return train_adata, test_adata

    ## get training dataset
    train_adata = load_PBMC_data.load_PBMC_batch1_data(data_dir, batch=batch)
    train_adata.obs.index.name = None

    ## get test dataset
    test_adata = load_PBMC_data.load_PBMC_batch1_data(data_dir, inds=ind)
    test_adata.obs.index.name = None

    ## downsample training dataset to S1 = 1,551 cells
    if sample:
        random.seed(dr_seed)
        sampled_cells = random.sample(list(train_adata.obs_names), k=1551)
        train_adata = train_adata[sampled_cells]

    ## curate given sub-cell types to major cell types
    train_adata = load_PBMC_data.curate_PBMC_demulx_celltypes(train_adata, celltype_gran)
    test_adata = load_PBMC_data.curate_PBMC_demulx_celltypes(test_adata, celltype_gran)

    ## feature selection
    train_adata, test_adata = feature_selection_train_test(train_adata, test_adata,
            result_dir, gene_no, select_on, select_method)

    ## scale and analze
    train_adata, test_adata = scale_and_visualize(train_adata, test_adata, result_dir)
    save_adata(train_adata, test_adata, result_dir)

    return train_adata, test_adata


# ---- process batch2 control and stimulated
def process_batch2_ctrl_stim(data_dir, result_dir, gene_no=1000, dr_seed=0, 
        cond1="control", cond2="stimulated", select_on="test", select_method="Seurat",
        purify=False, purify_method="distance", celltype_gran=1):
    '''Extract control as train adata, stimulated as test adata

    @train/test: control or stimultaed
    @select_on: "test" if selecting features from test; "train" if selecting features
                 from train.
    @select_method: Seurat/FEAST/None
    '''
    load_ind, train_adata, test_adata = load_adata(result_dir)
    if load_ind:
        return train_adata, test_adata

    train_adata = load_PBMC_data.load_PBMC_batch2_data(data_dir, condition=cond1)
    test_adata = load_PBMC_data.load_PBMC_batch2_data(data_dir, condition=cond2)

    print(train_adata)

    ## curate given sub-cell types to major cell types
    train_adata = load_PBMC_data.curate_PBMC_demulx_celltypes(train_adata, celltype_gran)
    test_adata = load_PBMC_data.curate_PBMC_demulx_celltypes(test_adata, celltype_gran)

    ## feature selection
    train_adata, test_adata = feature_selection_train_test(train_adata, test_adata,
            result_dir, gene_no, select_on, select_method)

    ## scale and analze
    train_adata, test_adata = scale_and_visualize(train_adata, test_adata, result_dir,
            purify=purify, purify_method=purify_method)
    save_adata(train_adata, test_adata, result_dir)

    return train_adata, test_adata

# ---- process PBMC batch1/2 individuals
def process_batch1_batch2_ind(data_dir, result_dir, input1, input2, 
        gene_no=1000, dr_seed=0, select_on="test", select_method="Seurat", 
        celltype_gran=1):
    '''Use individuals from one batch to predict another batch

    @input1/input2: can be batch1_indID/batc2_indID
    @celltype_gran: 0 major; 1 sub
    '''

    load_ind, train_adata, test_adata = load_adata(result_dir)
    if load_ind:
        return train_adata, test_adata

    ## split input1
    input1_list = input1.split('_')
    input1_batch = input1_list[0]
    if len(input1_list) > 1:
        input1_inds = '_'.join(input1_list[1:])
    else:
        input1_inds = None

    input2_list = input2.split('_')
    input2_batch = input2_list[0]
    if len(input2_list) > 1:
        input2_inds = '_'.join(input2_list[1:])
    else:
        input2_inds = None

    ## extract train and test adata according to batch information
    if input1_batch == "batch1":
        train_adata = load_PBMC_data.load_PBMC_batch1_data(data_dir, ind=input1_inds)
    elif "batch2" in input1_batch:
        cond = input1_batch.replace("batch2", "")
        if cond == "":
            train_adata = load_PBMC_data.load_PBMC_batch2_data(data_dir, ind=input1_inds)
        else:
            train_adata = load_PBMC_data.load_PBMC_batch2_data(data_dir, ind=input1_inds,
                condition=cond)

    if input2_batch == "batch1":
        test_adata = load_PBMC_data.load_PBMC_batch1_data(data_dir, ind=input2_inds)
    elif "batch2" in input2_batch:
        cond = input2_batch.replace("batch2", "")
        if cond == "":
            test_adata = load_PBMC_data.load_PBMC_batch2_data(data_dir, ind=input2_inds)
        else:
            test_adata = load_PBMC_data.load_PBMC_batch2_data(data_dir, ind=input2_inds,
                condition=cond)

    ## curate given sub-cell types to major cell types
    train_adata = load_PBMC_data.curate_PBMC_demulx_celltypes(train_adata, celltype_gran)
    test_adata = load_PBMC_data.curate_PBMC_demulx_celltypes(test_adata, celltype_gran)

    ## feature selection
    train_adata, test_adata = feature_selection_train_test(train_adata, test_adata,
            result_dir, gene_no, select_on, select_method)

    ## scale and analze
    train_adata, test_adata = scale_and_visualize(train_adata, test_adata, result_dir)
    save_adata(train_adata, test_adata, result_dir)

    return train_adata, test_adata


# ---- process PBMC 7 protocols
def process_PBMC_protocols_type(data_dir, result_dir, gene_no=1000, dr_seed=0,
        protocols1="Smart-seq2", protocols2="CEL-Seq2", exp="pbmc1",
        select_on="test", select_method="Seurat"):
    '''Compare protocols using PBMC1/PBMC2/mixture

    @ protocols1/protocols2: 10x Chromium (v2)/10x Chromium (v2) A/10x Chromium (v2) B/
        10x Chromium (v3)/CEL-Seq2/Drop-seq/Seq-Well/Smart-seq2/inDrops
    @ exp: pbmc1/pbmc2/pbmc1_pbmc2
    '''
    load_ind, train_adata, test_adata = load_adata(result_dir)
    if load_ind:
        return train_adata, test_adata

    ## load PBMC protocols adata
    adata = load_PBMC_data.load_PBMC_protocols_data(data_dir)
    exp_cells = adata.obs[adata.obs['Experiment'].isin(exp.split('_'))].index
    exp_adata = adata[exp_cells]

    train_cells = exp_adata.obs[exp_adata.obs['Method'].isin(protocols1.split('_'))].index
    train_adata = exp_adata[train_cells]
    test_cells = exp_adata.obs[exp_adata.obs['Method'].isin(protocols2.split('_'))].index
    test_adata = exp_adata[test_cells]

    # curate for common cell types
    common_celltypes = set(train_adata.obs["cell.type"]).intersection(set(test_adata.obs["cell.type"]))
    train_cells = train_adata.obs.loc[train_adata.obs["cell.type"].isin(common_celltypes)].index
    test_cells = test_adata.obs.loc[test_adata.obs["cell.type"].isin(common_celltypes)].index
    train_adata = train_adata[train_cells]
    test_adata = test_adata[test_cells]

    ## feature selection
    train_adata, test_adata = feature_selection_train_test(train_adata, test_adata,
            result_dir, gene_no, select_on, select_method)

    ## scale and analze
    train_adata, test_adata = scale_and_visualize(train_adata, test_adata, result_dir)
    save_adata(train_adata, test_adata, result_dir)

    return train_adata, test_adata

# ---- process PBMC 7 protocols batch
def process_PBMC_protocols_batch(data_dir, result_dir, gene_no=1000, dr_seed=0,
        batch1="pbmc1", batch2="pbmc2", protocol="Smart-seq2",
        select_on="test", select_method="Seurat"):
    '''Compare same protocol between batch1 and batch2

    @ protocol: 10x Chromium (v2)/10x Chromium (v2) A/10x Chromium (v2) B/
        10x Chromium (v3)/CEL-Seq2/Drop-seq/Seq-Well/Smart-seq2/inDrops
    @ batch1/batch2: pbmc1/pbmc2
    '''
    load_ind, train_adata, test_adata = load_adata(result_dir)
    if load_ind:
        return train_adata, test_adata

    ## load PBMC protocols adata
    adata = load_PBMCprotocol_data.load_PBMC_protocols_data(data_dir)
    pro_cells = adata.obs[adata.obs['Method'].isin(protocol.split('_'))].index
    pro_adata = adata[pro_cells]

    train_cells = pro_adata.obs[pro_adata.obs['Experiment'].isin([batch1])].index
    train_adata = pro_adata[train_cells]
    test_cells = pro_adata.obs[pro_adata.obs['Experiment'].isin([batch2])].index
    test_adata = pro_adata[test_cells]

    # curate for common cell types
    common_celltypes = set(train_adata.obs["cell.type"]).intersection(set(test_adata.obs["cell.type"]))
    train_cells = train_adata.obs.loc[train_adata.obs["cell.type"].isin(common_celltypes)].index
    test_cells = test_adata.obs.loc[test_adata.obs["cell.type"].isin(common_celltypes)].index
    train_adata = train_adata[train_cells]
    test_adata = test_adata[test_cells]

    ## feature selection
    train_adata, test_adata = feature_selection_train_test(train_adata, test_adata,
            result_dir, gene_no, select_on, select_method)

    ## scale and analze
    train_adata, test_adata = scale_and_visualize(train_adata, test_adata, result_dir)
    save_adata(train_adata, test_adata, result_dir)

    return train_adata, test_adata
