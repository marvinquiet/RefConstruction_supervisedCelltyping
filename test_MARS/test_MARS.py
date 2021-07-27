## should use mars environment

import torch
import os, sys, time
import numpy as np
import pandas as pd
import scanpy.api as sc
from anndata import AnnData
import anndata

from matplotlib import pyplot as plt
import matplotlib as mpl

sys.path.append("../")
sys.path.append('/homelocal/wma36/celltyping_refConstruct/test_MARS/mars')
from args_parser import get_parser
from model.mars import MARS
from model.experiment_dataset import ExperimentDataset
import warnings
warnings.filterwarnings('ignore')

## set seed
torch.manual_seed(0)
import random
random.seed(0)
np.random.seed(0)

## pytorch params
params, unknown = get_parser().parse_known_args()
params.cuda = True
device = 'cuda:0' if torch.cuda.is_available() and params.cuda else 'cpu'
params.device = device

def run_MARS(train_adata, test_adata, result_dir, prefix="MARS", 
        celltype_col="cell.type"):
    '''Run MARS pipeline
    '''

    adata = train_adata.concatenate(test_adata)

    from sklearn import preprocessing
    le = preprocessing.LabelEncoder()
    le.fit(train_adata.obs[celltype_col])
    y_train = np.array(le.transform(train_adata.obs[celltype_col]), dtype=np.int64)
    annotated = ExperimentDataset(train_adata.X, train_adata.obs_names, 
            train_adata.var_names, '0', y_train)
    
    y_test = np.array(le.transform(test_adata.obs[celltype_col]), dtype=np.int64)
    unannnotated = ExperimentDataset(test_adata.X, test_adata.obs_names, 
            test_adata.var_names, '1', y_test)
    
    pretrain_data = ExperimentDataset(test_adata.X, test_adata.obs_names, 
            test_adata.var_names, '1')
    n_clusters = len(np.unique(unannnotated.y))
    
    #params.epochs = 100 ## change to 50/100 and see the difference
    start = time.time()
    mars = MARS(n_clusters, params, [annotated], unannnotated, pretrain_data, hid_dim_1=1000, hid_dim_2=100)
    mars_adata, landmarks, scores = mars.train(evaluation_mode=True, save_all_embeddings=True) # evaluation mode
    end = time.time()
    with open(result_dir+os.sep+prefix+'_metrics.txt', 'w') as f:
        f.write("Acc:%s\n" % str(scores['accuracy']))
        f.write("ARI:%s\n" % str(scores['adj_rand']))
        f.write("macroF1:%s\n" % str(scores['f1_score']))
        f.write("runtime:%s\n" % str(end-start))
    
    ## give names to cell types (Aborted, very inaccurate)
    #keys = le.transform(le.classes_)
    #values = le.classes_
    #cell_type_name_map = dict()
    #for index, key in enumerate(keys):
    #    cell_type_name_map[key] = values[index]
    #interp_names = mars.name_cell_types(mars_adata, landmarks, cell_type_name_map)

    mars_adata.obs = mars_adata.obs.merge(adata.obs, 
            left_index=True, right_index=True, how="left")
    mars_adata = mars_adata[mars_adata.obs['dataset_batch'] == "test"]

    adata_mars = AnnData(mars_adata.obsm['MARS_embedding'])
    adata_mars.obs['MARS_labels'] = pd.Categorical(mars_adata.obs['MARS_labels'])
    adata_mars.obs['ground_truth'] = pd.Categorical(mars_adata.obs[celltype_col])
    adata_mars.obs.to_csv(result_dir+os.sep+prefix+"_predicted_obs.csv")

    sc.pp.neighbors(adata_mars, n_neighbors=30, use_rep='X')
    sc.tl.umap(adata_mars)
    sc.pl.umap(adata_mars, color=['ground_truth','MARS_labels'],size=15)
    plt.savefig(result_dir+os.sep+prefix+"_prediction_result.png")
    print("=== Finish visualizing..")

def process_data(data_dir, result_dir, args=None):
    '''Load processed datasets from other methods

    I do not want to repeat the configuration generation procedure here for every experiment
    '''
    train_adata, test_adata = None, None
    train_adata_file = result_dir+os.sep+"train_adata.h5ad"
    if os.path.exists(train_adata_file):
        train_adata = anndata.read_h5ad(train_adata_file)

    test_adata_file = result_dir+os.sep+"test_adata.h5ad"
    if os.path.exists(test_adata_file):
        test_adata = anndata.read_h5ad(test_adata_file)

    if train_adata is None or test_adata is None:
        sys.exit("Please run other methods to generate the train/test anndata.")

    ## find common celltypes between train and test
    common_celltypes = set(train_adata.obs["cell.type"]).intersection(set(test_adata.obs["cell.type"]))
    train_cells = train_adata.obs.loc[train_adata.obs["cell.type"].isin(common_celltypes)].index
    test_cells = test_adata.obs.loc[test_adata.obs["cell.type"].isin(common_celltypes)].index

    train_adata = train_adata[train_cells]
    test_adata = test_adata[test_cells]
    return train_adata, test_adata


if __name__ == '__main__':
    data_dir = "/home/wma36/gpu/data"

    ## parse arguments
    import argparse
    parser = argparse.ArgumentParser(description="Celltyping pipeline.")
    parser.add_argument('data_source', help="Load which dataset",
        choices=[
            'PBMC_batch1_ind', 'PBMC_batch1_ABC', 'PBMC_batch2', 'PBMC_batch1_batchtoind',
            'PBMC_protocols_pbmc1', 'PBMC_protocols_batch_smart', 
            "PBMC_Zheng_FACS", "PBMC_Zheng_FACS_curated", "PBMC_cross", 
            'pancreas', 'pancreas_seg_cond', "pancreas_custom", "pancreas_seg_mix", "pancreas_multi_to_multi", 
            "mousebrain_FC", "mousebrain_FC_sub", "mousebrain_HC", "mousebrain_HC_sub", "mousebrain_region",
            "mousebrain_FC_stage", "mousebrain_FC_stage_sub", "mousebrain_FC_datasets", 
            "mousebrain_FC_multiinds", "mousebrain_FC_multiinds_sample", 
            "mousebrain_FC_multiinds_sub", "mousebrain_FC_multiinds_sub_sample",
            "mousebrain_FC_datasets_multiinds", "mousebrain_FC_datasets_multiinds_sample",
            "allenbrain_ss", "allenbrain_10x", "allenbrain_cross"
            ])

    parser.add_argument('--select_on', help="Feature selection on train or test, or None of them",
        choices=['train', 'test'])
    parser.add_argument('--select_method', help="Feature selection method, Seurat/FEAST or None",
            choices=['Seurat', 'FEAST', 'F-test'])
    parser.add_argument('--n_features', help="Number of features selected",
            default=1000, type=int)
    parser.add_argument('--train', help="Specify which as train", required=True)
    parser.add_argument('--test', help="Specify which as test", required=True)

    args = parser.parse_args()
    if args.data_source in ["PBMC_batch1_ind", "PBMC_batch1_ABC", "PBMC_batch2", "PBMC_batch1_batchtoind", 
            "PBMC_protocols_pbmc1", "PBMC_protocols_batch_smart"]:
        pipeline_dir = "/home/wma36/gpu/celltyping_refConstruct/pipelines/result_PBMC_collections"

    if args.data_source in ["PBMC_protocols_pbmc1", "PBMC_protocols_batch_smart"]:
        pipeline_dir = "/home/wma36/gpu/celltyping_refConstruct/pipelines/result_PBMC_protocols_collections"

    if args.data_source in ["PBMC_Zheng_FACS", "PBMC_Zheng_FACS_curated", "PBMC_cross"]:
        pipeline_dir = "/home/wma36/gpu/celltyping_refConstruct/pipelines/result_PBMC_Zheng_collections"

    if args.data_source in ["pancreas", "pancreas_seg_cond", "pancreas_custom", 
            "pancreas_seg_mix", "pancreas_multi_to_multi"]:
        pipeline_dir = "/home/wma36/gpu/celltyping_refConstruct/pipelines/result_Pancreas_collections"

    if args.data_source in ["mousebrain_FC", "mousebrain_FC_sub", "mousebrain_HC", "mousebrain_HC_sub",
            "mousebrain_region", "mousebrain_FC_stage", "mousebrain_FC_stage_sub", "mousebrain_FC_datasets",
            "mousebrain_FC_datasets_multiinds", "mousebrain_FC_datasets_multiinds_sample", 
            "mousebrain_FC_multiinds", "mousebrain_FC_multiinds_sub", "mousebrain_FC_multiinds_sample",
            "mousebrain_FC_multiinds_sub_sample"]:
        pipeline_dir = "/home/wma36/gpu/celltyping_refConstruct/pipelines/result_Mousebrain_collections"

    if args.data_source in ["allenbrain_ss", "allenbrain_10x", "allenbrain_cross"]:
        pipeline_dir = "/home/wma36/gpu/celltyping_refConstruct/pipelines/result_Allenbrain_collections"

    result_prefix = pipeline_dir+os.sep+"result_"+args.data_source+'_'+\
        args.train+'_to_'+args.test
    os.makedirs(result_prefix, exist_ok=True)

    ## create file directory 
    if args.select_on is None and args.select_method is None:
        result_dir = result_prefix+os.sep+"no_feature"
    else:
        result_dir = result_prefix+os.sep+args.select_method+'_'+\
                str(args.n_features)+'_on_'+args.select_on
    os.makedirs(result_dir, exist_ok=True)

    train_adata, test_adata = process_data(data_dir, result_dir, args=args)
    run_MARS(train_adata, test_adata, result_dir)

