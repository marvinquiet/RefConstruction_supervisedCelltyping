import os, sys
import numpy as np
import pandas as pd
import scipy
import scanpy.api as sc
import anndata
import random

from preprocess.process_train_test_data import *

from preprocess import load_mousebrain_data
from preprocess import load_mouseFC_data
from preprocess import load_mouseprotocol_data

# ==== process mouse data
def process_mousebrain_data(data_dir, result_dir, ind1, ind2, region="FC",
        gene_no=1000, dr_seed=0, select_on="test", select_method="Seurat",
        celltype_gran=0, purify=False, purify_method="distance", curate=False):
    '''Process Frontal Cortex from one ind to another
    @ind: P60FCAldh1l1Rep1/P60FCCx3cr1Rep1/P60FCRep1/P60FCRep2/P60FCRep3/P60FCRep4/P60FCRep6
    @celltype_gran: granularity of cell types. 0: major cell types, 1: sub-celltypes
    @purify: whether to purify the training dataset
    @purify_method: purify using distance or SVM
    @curate: whether to curate mouse brain data by combining different neurons/interneurons
    '''
    load_ind, train_adata, test_adata = load_adata(result_dir)
    if load_ind:
        return train_adata, test_adata

    train_adata = load_mousebrain_data.load_brain_adata(data_dir, region=region,
            ind=ind1, celltype_gran=celltype_gran, curate=curate)
    test_adata = load_mousebrain_data.load_brain_adata(data_dir, region=region,
            ind=ind2, celltype_gran=celltype_gran, curate=curate)

    ## feature selection
    train_adata, test_adata = feature_selection_train_test(train_adata, test_adata,
            result_dir, gene_no, select_on, select_method)

    ## scale and analze
    train_adata, test_adata = scale_and_visualize(train_adata, test_adata, result_dir,
            purify=purify, purify_method=purify_method)
    save_adata(train_adata, test_adata, result_dir)

    return train_adata, test_adata

def process_mousebrain_multiinds_data(data_dir, result_dir, sample_ind, pred_ind, 
        region="FC", gene_no=1000, dr_seed=0, sample_seed=0, select_on="test", select_method="Seurat",
        celltype_gran=0, sample=False):
    '''Process mousebrain data from multiple individuals to another

    @sample_ind: downsample to a certain individual number
    @pred_ind: the predicted individual
    '''
    load_ind, train_adata, test_adata = load_adata(result_dir)
    if load_ind:
        return train_adata, test_adata

    ## individuals for brain regions
    FC_inds = ["P60FCAldh1l1Rep1", "P60FCCx3cr1Rep1", "P60FCRep1", "P60FCRep2",
            "P60FCRep3", "P60FCRep4", "P60FCRep6"]
    HC_inds = ["P60HippoRep1", "P60HippoRep2", "P60HippoRep3",
            "P60HippoRep4", "P60HippoRep5", "P60HippoRep6"]

    if region == "FC":
        exclude_list = FC_inds
    elif region == "HC":
        exclude_list = HC_inds

    exclude_list.remove(pred_ind)
    exclude_inds = '_'.join(exclude_list)

    train_adata = load_mousebrain_data.load_brain_adata(data_dir, region=region,
            ind=exclude_inds, celltype_gran=celltype_gran)

    ## sample to number of cells when compared to the sample_id
    if sample:
        avg_number = train_adata.shape[0] // len(exclude_list)

        print("=== Downsample to number:", avg_number)
        random.seed(sample_seed)
        sampled_cells = random.sample(list(train_adata.obs_names), k=avg_number)
        train_adata = train_adata[sampled_cells]

    test_adata = load_mousebrain_data.load_brain_adata(data_dir, region=region,
            ind=pred_ind, celltype_gran=celltype_gran)

    ## feature selection
    train_adata, test_adata = feature_selection_train_test(train_adata, test_adata,
            result_dir, gene_no, select_on, select_method)

    ## scale and analze
    train_adata, test_adata = scale_and_visualize(train_adata, test_adata, result_dir)
    save_adata(train_adata, test_adata, result_dir)

    return train_adata, test_adata

def process_mousebrain_regions(data_dir, result_dir, region1, region2,
        gene_no=1000, dr_seed=0, select_on="test", select_method="Seurat",
        celltype_gran=0):
    ''' Use one region to predict another
    can allow region_individual input

    @celltype_gran: with sub-cell types and check ARI
    '''
    load_ind, train_adata, test_adata = load_adata(result_dir)
    if load_ind:
        return train_adata, test_adata

    ## get region information and individuals
    regions1 = region1.split('_')
    regions2 = region2.split('_')

    region1 = regions1[0]
    if len(regions1) > 1:
        regions1_inds = '_'.join(regions1[1:])
    else:
        regions1_inds = None

    region2 = regions2[0]
    if len(regions2) > 1:
        regions2_inds = '_'.join(regions2[1:])
    else:
        regions2_inds = None

    ## get train/test anndata
    train_adata = load_mousebrain_data.load_brain_adata(data_dir, region=region1,
            ind=regions1_inds, celltype_gran=celltype_gran, curate=True)
    test_adata = load_mousebrain_data.load_brain_adata(data_dir, region=region2,
            ind=regions2_inds, celltype_gran=celltype_gran, curate=True)
 
    ## feature selection
    train_adata, test_adata = feature_selection_train_test(train_adata, test_adata,
            result_dir, gene_no, select_on, select_method)

    ## scale and analze
    train_adata, test_adata = scale_and_visualize(train_adata, test_adata, result_dir)
    save_adata(train_adata, test_adata, result_dir)

    return train_adata, test_adata


## === process mouse FC data
def process_mousebrain_FC_stage(data_dir, result_dir, stage1, stage2,
        gene_no=1000, dr_seed=0, select_on="test", select_method="Seurat",
        celltype_gran=0):
    '''Use one stage to predict another stage
    '''
    load_ind, train_adata, test_adata = load_adata(result_dir)
    if load_ind:
        return train_adata, test_adata
 
    train_adata = load_mouseFC_data.load_FC_adata(data_dir, devstage=stage1,
            celltype_gran=celltype_gran)
    test_adata = load_mouseFC_data.load_FC_adata(data_dir, devstage=stage2,
            celltype_gran=celltype_gran)

    if 1 == celltype_gran:
        ## find common sub-celltypes
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

def process_mousebrain_FC_datasets(data_dir, result_dir, dataset1_input, 
        dataset2_input, gene_no=1000, dr_seed=0, sample_seed=0, 
        select_on="test", select_method="Seurat", 
        celltype_gran=0, sample=False):
    '''Use FC from one dataset to predict another

    @dataset1_input: "FC_inds", if sample == True, 
        the inds indicate the number of cells to downsample
    @dataset2_input: "Adult_inds"
    the above two parameters can be interchangable
    '''
    load_ind, train_adata, test_adata = load_adata(result_dir)
    if load_ind:
        return train_adata, test_adata

    dataset1_split = dataset1_input.split('_')
    dataset2_split = dataset2_input.split('_')

    ## identify where the dataset points to
    regions = ["FC", "HC", "CB", "STR", "SN", "PC"]
    stages = ["Adult", "P21"]

    dataset1_input = dataset1_split[0]
    if len(dataset1_split) > 1:
        dataset1_inds = dataset1_split[1:]
        dataset1_inds = '_'.join(dataset1_inds)
    else:
        dataset1_inds = None

    dataset2_input = dataset2_split[0]
    if len(dataset2_split) > 1:
        dataset2_inds = dataset2_split[1:]
        dataset2_inds = '_'.join(dataset2_inds)
    else:
        dataset2_inds = None

    if dataset1_input in regions:
        train_adata = load_mousebrain_data.load_brain_adata(data_dir, 
                region=dataset1_input, ind=dataset1_inds, curate=True)
    else:
        info = dataset1_input.split('-')
        if info[0] in stages:
            if len(info) > 1:
                train_adata = load_mouseFC_data.load_FC_adata(data_dir, devstage=info[0],
                    treatment=info[1], curate=True)
            else:
                train_adata = load_mouseFC_data.load_FC_adata(data_dir, devstage=info[0],
                    ind=dataset1_inds, curate=True)

    if dataset2_input in regions:
        test_adata = load_mousebrain_data.load_brain_adata(data_dir,
                region=dataset2_input, ind=dataset2_inds, curate=True)
    else:
        info = dataset2_input.split('-')
        if info[0] in stages:
            if len(info) > 1:
                train_adata = load_mouseFC_data.load_FC_adata(data_dir, devstage=info[0],
                    treatment=info[1], curate=True)
            else:
                train_adata = load_mouseFC_data.load_FC_adata(data_dir, devstage=info[0],
                    ind=dataset2_inds, curate=True)
        
    if sample:
        if dataset1_input == "FC":  ## 7 individuals in mousebrian FC
            inds_no = 7
        
        if dataset1_input == "Adult": ## 6 individuals in mousebrain pFC
            inds_no = 6

        avg_number = train_adata.shape[0] // inds_no
        random.seed(sample_seed)
        sampled_cells = random.sample(list(train_adata.obs_names), k=avg_number)
        train_adata = train_adata[sampled_cells]

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

## Process Mouse cortex protocols datasets
def process_mousecortex_protocols_mouseFC_data(data_dir, result_dir, 
        input1, input2, gene_no=1000, dr_seed=0,
        select_on="test", select_method="Seurat", 
        curate=False):
    '''process both mousecortex protocols and mouse FC data to predict

    Note here we can only do major cell types because Cortex data does not 
    provide sub-cell types

    @input1/input2: MouseProtocol+exp+protocols_inds, mouseFC_inds
    @curate: whether to curate the datasets
    '''

    load_ind, train_adata, test_adata = load_adata(result_dir)
    if load_ind:
        return train_adata, test_adata

    ## load train and test adata
    input1_list = input1.split('_')
    input1_dataset = input1_list[0]
    if len(input1_list) > 1:
        input1_inds = '_'.join(input1_list[1:])
    else:
        input1_inds = None

    input2_list = input2.split('_')
    input2_dataset = input2_list[0]
    if len(input2_list) > 1:
        input2_inds = '_'.join(input2_list[1:])
    else:
        input2_inds = None


    ## extract train and test adata accordingly
    if "mouseFC" in input1_dataset:
        train_adata = load_mousebrain_data.load_brain_adata(data_dir, region="FC",
                ind=input1_inds, celltype_gran=0, curate=True)
    elif "MouseProtocol" in input1_dataset:
        info = input1_dataset.split('+')
        exp = info[1]
        protocol = info[2]

        if exp == "Both":
            train_adata = load_mouseprotocol_data.load_mouseprotocol_adata(data_dir,
                    exp=None, protocol=protocol, curate=True)
        else:
            train_adata = load_mouseprotocol_data.load_mouseprotocol_adata(data_dir,
                    exp=exp, protocol=protocol, curate=True)

    if "mouseFC" in input2_dataset:
        test_adata = load_mousebrain_data.load_brain_adata(data_dir, region="FC",
                ind=input2_inds, celltype_gran=0, curate=True)
    elif "MouseProtocol" in input2_dataset:
        info = input2_dataset.split('+')
        exp = info[1]
        protocol = info[2]

        if exp == "Both":
            test_adata = load_mouseprotocol_data.load_mouseprotocol_adata(data_dir,
                    exp=None, protocol=protocol, curate=True)
        else:
            test_adata = load_mouseprotocol_data.load_mouseprotocol_adata(data_dir,
                    exp=exp, protocol=protocol, curate=True)

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

def process_FCdatasets_mouseFC_data(data_dir, result_dir, datasets_ref, mousebrain_tgt,
        gene_no=1000, dr_seed=0, select_on="test", select_method="Seurat",
        celltype_gran=0, curate=False):
    '''Combine mouse brain protocol + pFC datasets together to predict the 
    mouse FC target individual

    @datasets_ref: MouseProtocol+Both+DroNc-seq_pFC+Adult+Saline
    @mousebrain_tgt: mousebrain dataset with certain targets
    '''

    load_ind, train_adata, test_adata = load_adata(result_dir)
    if load_ind:
        return train_adata, test_adata

    datasets_list = datasets_ref.split('_')
    adata_list = []
    for dataset in datasets_list:
        infos = dataset.split('+')
        dataset_name = infos[0]
        if dataset_name == "MouseProtocol":
            exp = infos[1]
            protocol = infos[2]
            if exp == "Both":
                adata = load_mouseprotocol_data.load_mouseprotocol_adata(data_dir,
                        exp=None, protocol=protocol, curate=curate)
            else:
                adata = load_mouseprotocol_data.load_mouseprotocol_adata(data_dur,
                        exp=exp, protocol=protocol, curate=curate)
        elif dataset_name == "pFC":
            devstage = infos[1]
            treatment = infos[2]
            adata = load_mouseFC_data.load_FC_adata(data_dir, devstage=devstage,
                    treatment=treatment, celltype_gran=celltype_gran, curate=curate)
        adata_list.append(adata)

    ## find common genes and concatenate adata
    common_columns = ['barcode', 'nGene', 'nUMI', 'percent.mito', 'cell.type']
    for adata in adata_list:
        adata_obs = adata.obs[common_columns]
        adata.obs = adata_obs
    train_adata = anndata.AnnData.concatenate(*adata_list, join='inner')
    del adata_list  ## release space

    ## get target individuals list
    tgt_list = mousebrain_tgt.split('_')
    region = tgt_list[0]
    if len(tgt_list) > 1:
        tgt_inds = '_'.join(tgt_list[1:])
    else:
        tgt_inds = None

    test_adata = load_mousebrain_data.load_brain_adata(data_dir, region=region,
            ind=tgt_inds, celltype_gran=celltype_gran, curate=curate)

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
