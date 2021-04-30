import os, sys
import numpy as np
import pandas as pd
import scipy
import scanpy.api as sc
import anndata

from preprocess.process_train_test_data import *

from preprocess import load_pancreatic_data

## --- process pancreatic data
def process_pancreas(data_dir, result_dir, dataset1="seg", dataset2="xin", 
        gene_no=1000, dr_seed=0, select_on="test", select_method="Seurat"):
    '''Use one dataset to predict another
    '''

    load_ind, train_adata, test_adata = load_adata(result_dir)
    if load_ind:
        return train_adata, test_adata

    train_adata = load_pancreatic_data.load_panc_data(data_dir, source=dataset1)
    test_adata = load_pancreatic_data.load_panc_data(data_dir, source=dataset2)

    ## feature selection
    train_adata, test_adata = feature_selection_train_test(train_adata, test_adata,
            result_dir, gene_no, select_on, select_method)

    ## scale and analze
    train_adata, test_adata = scale_and_visualize(train_adata, test_adata, result_dir)
    save_adata(train_adata, test_adata, result_dir)

    return train_adata, test_adata

def process_pancreas_seg_cond(data_dir, result_dir, cond1="Healthy", cond2="T2D",
        gene_no=1000, dr_seed=0, select_on="test", select_method="Seurat"):
    '''Process seg data with conditions
    '''
    load_ind, train_adata, test_adata = load_adata(result_dir)
    if load_ind:
        return train_adata, test_adata

    train_adata = load_pancreatic_data.load_seg_data(data_dir, condition=cond1)
    test_adata = load_pancreatic_data.load_seg_data(data_dir, condition=cond2)

    ## feature selection
    train_adata, test_adata = feature_selection_train_test(train_adata, test_adata,
            result_dir, gene_no, select_on, select_method)

    ## scale and analze
    train_adata, test_adata = scale_and_visualize(train_adata, test_adata, result_dir)
    save_adata(train_adata, test_adata, result_dir)

    return train_adata, test_adata


def process_pancreas_custom(data_dir, result_dir, seg_cond="Healthy", dataset="Muraro",
        gene_no=1000, dr_seed=0, select_on="test", select_method="Seurat"):
    '''Process seg Healthy and muraro
    '''
    load_ind, train_adata, test_adata = load_adata(result_dir)
    if load_ind:
        return train_adata, test_adata

    train_adata = load_pancreatic_data.load_seg_data(data_dir, condition=seg_cond)
    test_adata = load_pancreatic_data.load_panc_data(data_dir, source=dataset)

    ## feature selection
    train_adata, test_adata = feature_selection_train_test(train_adata, test_adata,
            result_dir, gene_no, select_on, select_method)

    ## scale and analze
    train_adata, test_adata = scale_and_visualize(train_adata, test_adata, result_dir)
    save_adata(train_adata, test_adata, result_dir)

    return train_adata, test_adata

def process_pancreas_seg_mix(data_dir, result_dir, main_cond="Healthy", pred_cond="T2D",
        gene_no=1000, dr_seed=0, select_on="test", select_method="Seurat"):
    '''Processing seg mix to predict another
    '''
    load_ind, train_adata, test_adata = load_adata(result_dir)
    if load_ind:
        return train_adata, test_adata

    adata = load_pancreatic_data.load_panc_data(data_dir, source="seg")
    adata.obs['condition'] = ["T2D" if "T2D" in cell else "Healthy" for cell in adata.obs.index]
    adata.obs['ind'] = [cell.split('_')[0] for cell in adata.obs.index]
    pred_inds = adata.obs[adata.obs['condition'] == pred_cond]["ind"].tolist()
    pred_inds.sort() ## sort for reproducibility in R code

    pred_ind = pred_inds[0] ## T2D ind: HP1504101T2D; Healthy ind: AZ
    print("===In seg mix, predicting:", pred_ind)

    main_cells = adata.obs[adata.obs["ind"] != pred_ind].index
    pred_cells = adata.obs[adata.obs["ind"] == pred_ind].index

    train_adata = adata[main_cells]
    test_adata = adata[pred_cells]

    ## feature selection
    train_adata, test_adata = feature_selection_train_test(train_adata, test_adata,
            result_dir, gene_no, select_on, select_method)

    ## scale and analze
    train_adata, test_adata = scale_and_visualize(train_adata, test_adata, result_dir)
    save_adata(train_adata, test_adata, result_dir)

    return train_adata, test_adata

def process_pancreas_multi_to_multi(data_dir, result_dir, cond1, cond2,
        gene_no=1000, dr_seed=0, select_on="test", select_method="Seurat"):
    '''Processing multi datasets combinations to predict other 
    '''
    load_ind, train_adata, test_adata = load_adata(result_dir)
    if load_ind:
        return train_adata, test_adata

    cond1_adata,cond2_adata = [], []
    for cond in cond1.split('_'):
        cond1_adata.append(load_pancreatic_data.load_panc_data(data_dir, source=cond))
    for cond in cond2.split('_'):
        cond2_adata.append(load_pancreatic_data.load_panc_data(data_dir, source=cond))

    train_adata = anndata.AnnData.concatenate(*cond1_adata,join='inner', 
            batch_key="data_source", batch_categories= cond1.split('_'))
    test_adata = anndata.AnnData.concatenate(*cond2_adata,join='inner', 
            batch_key="data_source", batch_categories= cond2.split('_'))

    ## feature selection
    train_adata, test_adata = feature_selection_train_test(train_adata, test_adata,
            result_dir, gene_no, select_on, select_method)

    ## scale and analze
    train_adata, test_adata = scale_and_visualize(train_adata, test_adata, result_dir)
    save_adata(train_adata, test_adata, result_dir)

    return train_adata, test_adata

