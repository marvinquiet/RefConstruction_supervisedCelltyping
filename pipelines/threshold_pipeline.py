'''
Configuration generation for running pipeline on leaving one/two cell type out and tresholding the result
'''
import os, sys, time
import argparse
import numpy as np
import pandas as pd
import scipy
import tensorflow as tf

from sklearn.preprocessing import OneHotEncoder
from sklearn import metrics

import matplotlib.pyplot as plt
import seaborn as sns
import scanpy.api as sc
 
import random

#Set seeds
RANDOM_SEED=0
random.seed(RANDOM_SEED)
np.random.seed(RANDOM_SEED)
tf.set_random_seed(RANDOM_SEED)

## every time reset all graphs and run something new
tf.compat.v1.reset_default_graph()

from pipelines import method_utils, dataloading_utils
from preprocess.process_train_test_data import *

def run_MLP_pipeline(args, train_adata, test_adata, result_dir, prefix=""):
    ''' Run MLP pipeline for left out cell types

    @arg: argparse object
        - args.leaveout: indicate up to how many cell types will be randomly left out
        - args.threshold: use threshold to decide unassigned cells
    '''
    batch_size = 128
    celltype_cols = "cell.type"
    ## Hyperparameters for network
    if train_adata.shape[0] >= 5000:
        ## === parameters for mousebrain (high cell number)
        dims = [128, 32]
        MLP_dims = [128, 64, 32, 16, 8]
    else:
        ## === parameters for PBMC datasets (low cell number)
        dims = [16]
        MLP_dims = [64, 16]

    ## OneHotEncoding the celltypes
    enc = OneHotEncoder(handle_unknown='ignore')
    if scipy.sparse.issparse(train_adata.X):
        x_train = train_adata.X.toarray()
    else:
        x_train = train_adata.X
    y_train = enc.fit_transform(train_adata.obs[[celltype_cols]]).toarray()
    if scipy.sparse.issparse(test_adata.X):
        x_test = test_adata.X.toarray()
    else:
        x_test = test_adata.X

    ### --- run MLP
    print("\n\n=== MLP\n")
    start = time.time()
    y_pred = method_utils.run_MLP(x_train, y_train, x_test, 
            dims=MLP_dims, batch_size=batch_size, seed=RANDOM_SEED) ## run MLP
    end = time.time()
    print("\n\n=== Run time:", end-start)

    ### draw unassigned cells lower than a certain threshold
    thres = args.threshold
    test_adata.obs['pred_celltypes'] = 'unassigned'
    assigned_idx = np.where(np.amax(y_pred, 1) >= thres)[0]
    unassigned_idx = np.where(np.amax(y_pred, 1) < thres)[0]
    ## print out unassigned cells number
    print(test_adata.obs[celltype_cols].value_counts())
    print("Unassigned cells number:", len(unassigned_idx))
    print(test_adata[unassigned_idx].obs[celltype_cols].value_counts())
    ## get assigned cell labels
    assigned_labels = y_pred[assigned_idx].argmax(1)
    n_clusters = len(set(train_adata.obs[celltype_cols]))
    assigned_onehot = np.zeros((assigned_labels.size, n_clusters))
    assigned_onehot[np.arange(assigned_labels.size), assigned_labels] = 1
    assigned_celltypes = enc.inverse_transform(assigned_onehot)
    assigned_cells = test_adata[assigned_idx].obs_names
    test_adata.obs.loc[assigned_cells, 'pred_celltypes'] = assigned_celltypes
    test_adata.obs.to_csv(result_dir+os.sep+prefix+"_predicted_obs.csv")
    ## visualization of unassigned cells
    plt.figure(args.leaveout+1)
    sc.pl.tsne(test_adata, color=[celltype_cols, "pred_celltypes"], size=15)
    plt.savefig(result_dir+os.sep+prefix+'_'+str(args.threshold)+"_prediction_result.png")
    print("=== Finish visualizing..")
    return y_pred

if __name__ == "__main__":
    data_dir = "/home/wma36/gpu/data"

    ## parse arguments
    parser = argparse.ArgumentParser(description="Test for saturation pipeline.")
    parser.add_argument('data_source', help="Load which dataset")
    parser.add_argument('--train', help="Specify which as train", required=True)
    parser.add_argument('--test', help="Specify which as test", required=True)
    parser.add_argument('--sample_seed', help="Downsample seed in combined individual effect", 
            default=0, type=int)
    parser.add_argument('--leaveout', help="Number of cells being left out", 
            default=0, type=int)
    parser.add_argument('--threshold', help="Threshold for setting cells as unassigned",
            default=0.9, type=float)

    args = parser.parse_args()
    ## F-test on train + MLP
    args.method = "MLP"
    args.select_on = "train"
    args.select_method = "F-test"
    args.n_features = 1000

    pipeline_dir = "/home/wma36/gpu/celltyping_refConstruct/pipelines/result_threshold_collections"
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

    load_ind, train_adata, test_adata = load_adata(result_dir)
    if not load_ind:
        train_adata, test_adata = dataloading_utils.load_PBMC_adata(
            data_dir, result_dir, args=args)
        if train_adata is None or test_adata is None: 
            train_adata, test_adata = dataloading_utils.load_Pancreas_adata(
                data_dir, result_dir, args=args)
        if train_adata is None or test_adata is None:
            train_adata, test_adata = dataloading_utils.load_Mousebrain_adata(
                data_dir, result_dir, args=args)
        if train_adata is None or test_adata is None:
            sys.exit("Please check your data source to make sure it matches an input.")

        ## whether to purify reference dataset
        purify_method = ""
        if "purify_dist" in args.data_source:
            purify_method = "distance"
        elif "purify_SVM" in args.data_source:
            purify_method = "SVM"

        train_adata, test_adata = dataloading_utils.process_loaded_data(
                train_adata, test_adata, result_dir, args=args, purify_method=purify_method)
        print("Train anndata: \n", train_adata)
        print("Test anndata: \n", test_adata)

    celltype_cols = "cell.type"
    celltypes = list(set(train_adata.obs[celltype_cols]))
    celltypes.sort()
    if args.leaveout > 0:
        prefix = args.method+'_'+str(args.leaveout)
        plt.figure(0)
        celltypes_palette = sns.color_palette(None, args.leaveout)
        plt.plot([0, 1], [0, 1], color='black', lw=2, linestyle='--')
        plt.xlim([0.0, 1.05])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title(args.data_source.replace('_', ' '))

        ## remove 1,2..args.leaveout cell types
        original_train_adata = train_adata.copy()
        for i in range(args.leaveout):
            leaveout_celltypes = celltypes[:i+1]
            remained_idx = ~original_train_adata.obs[celltype_cols].isin(leaveout_celltypes)
            train_adata = original_train_adata[remained_idx]
            y_pred = run_MLP_pipeline(args, train_adata.copy(), test_adata.copy(), result_dir, prefix=prefix+'_'+str(i+1))

            max_y_pred = y_pred.max(axis=1) ## get max probability of each row
            y_ind = (~test_adata.obs[celltype_cols].isin(leaveout_celltypes)).astype(int).tolist()
            fpr, tpr, thresholds = metrics.roc_curve(y_ind, max_y_pred, pos_label=1)
            auc = metrics.auc(fpr, tpr)
            label = '%s leftout \n (AUC = %0.2f)' % (','.join(leaveout_celltypes), auc)
            plt.figure(0)
            plt.plot(fpr, tpr, color=celltypes_palette[i],
                    lw=2, label=label)
        plt.legend(loc="lower right")
        plt.savefig(result_dir+os.sep+prefix+'_auroc.png')
    else:
        run_MLP_pipeline(args, train_adata.copy(), test_adata.copy(), result_dir, prefix=args.method)
 
