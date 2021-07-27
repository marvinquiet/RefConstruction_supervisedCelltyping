'''
Two-stage pipeline for celltyping, e.g. mousebrain with major cell types and 
minor cell types
'''
import os, sys, math, time
import anndata
import numpy as np
import pandas as pd
import tensorflow as tf

from sklearn.utils import shuffle
from sklearn.preprocessing import OneHotEncoder
from sklearn import metrics

import random
from numpy.random import seed
from tensorflow import set_random_seed

# GPU settings and reproducible
#os.environ["CUDA_VISIBLE_DEVICES"]="0"

#Set seeds
RANDOM_SEED=0
random.seed(RANDOM_SEED)
np.random.seed(RANDOM_SEED)
set_random_seed(RANDOM_SEED)

## every time reset all graphs and run something new
tf.compat.v1.reset_default_graph()

from pipelines.pipeline import *

def prediction(method, train_adata, test_adata):
    print("Train anndata: \n", train_adata)
    print("Test anndata: \n", test_adata)

    ## === first stage in predicting
    ## Hyperparameters for network
    if train_adata.shape[0] >= 5000:
        ## === parameters for mousebrain (high cell number)
        dims = [128, 32]
        MLP_dims = [128, 64, 32, 16, 8]
        batch_size = 128
    else:
        ## === parameters for PBMC datasets (low cell number)
        dims = [16]
        MLP_dims = [64, 16]
        batch_size = 32

    celltype_cols = "cell.type"
    n_clusters = len(set(train_adata.obs[celltype_cols]))
    ## OneHotEncoding the celltypes
    enc = OneHotEncoder(handle_unknown='ignore')
    x_train = np.array(train_adata.X)
    y_train = enc.fit_transform(train_adata.obs[[celltype_cols]]).toarray()

    x_test = np.array(test_adata.X)
    y_test = enc.transform(test_adata.obs[[celltype_cols]]).toarray()
 
    if "MLP" == method:
        graph = tf.get_default_graph()
        with graph.as_default():
            y_pred = run_MLP(x_train, y_train, x_test, 
                    dims=MLP_dims, batch_size=batch_size) ## run MLP
        y_pred = y_pred.argmax(1)

    if "GEDFN" == method or "DFN" == method:
        ## load partition dataframe
        partition_df = pd.read_csv(data_dir+os.sep+"GeneNetworks/HomoSapients_htb_hq_adj.txt", 
                sep=" ", index_col=0)
        partition_df.index = partition_df.index.str.upper()
        partition_df.columns = partition_df.columns.str.upper()

        ## read features
        with open(result_dir+os.sep+"features.txt", 'r') as f:
            features = f.read().splitlines()

        common_genes = set(features).intersection(set(partition_df.index))
        common_genes = list(common_genes)
        common_genes.sort() ## for reproducibility
        partition = partition_df.loc[common_genes, common_genes]
        partition = partition.to_numpy()  ## turn dataframe to numpy

        ## input intersect with partition
        x_train = np.array(train_adata[:, common_genes].X)
        x_test = np.array(test_adata[:, common_genes].X)
        print("After intersecting with partition: ", x_train.shape, x_test.shape)
        y_train = enc.fit_transform(train_adata.obs[[celltype_cols]]).toarray()
        y_test = enc.transform(test_adata.obs[[celltype_cols]]).toarray()

        ### --- run GEDFN
        graph = tf.get_default_graph()
        with graph.as_default():
            if "GE" in method:
                y_pred, var_importance = run_GEDFN(x_train, y_train, x_test, partition, 
                                         dims=dims, batch_size=batch_size, GE=True) ## run GEDFN
            else:
                y_pred, var_importance = run_GEDFN(x_train, y_train, x_test, partition, 
                                         dims=dims, batch_size=batch_size, GE=False) ## run DFN
        tf.keras.backend.clear_session()

        y_pred = y_pred.argmax(1)

    if "ItClust" == method:
        graph = tf.get_default_graph()
        with graph.as_default():
            ## --- run ItClust
            y_pred = run_ItClust(x_train, y_train, x_test, 
                            dims=dims, batch_size=batch_size) ## runItClust
        y_pred = y_pred.argmax(1)
 
    if "SVM_RBF" == method:
        ## --- run SVM
        y_pred = run_SVM(x_train, y_train.argmax(1), x_test)

    if "SVM_linear" == method:
        ## --- run SVM
        y_pred = run_SVM(x_train, y_train.argmax(1), x_test, kernel="linear")

    if "RF" == method:
        ## --- run RF
        y_pred = run_RF(x_train, y_train, x_test)
        y_pred = y_pred.argmax(1)

    return y_pred, enc


def process_mousebrain(args, data_dir, celltype_cols="cell.type"):
    '''Predict mousebrain two-stage

    @args: arguments from argparser
    @data_dir: directory where data is stored, e.g. /home/wma36/gpu/data
    '''
    ## === mouse adult brain
    from preprocess import load_mousebrain_data

    region = args.data_source[-2:]
    train_adata = load_mousebrain_data.load_brain_adata(data_dir, region=region,
            ind=args.train, celltype_gran=0)
    test_adata = load_mousebrain_data.load_brain_adata(data_dir, region=region,
            ind=args.test, celltype_gran=0)
    test_adata.obs["pred_sub_celltypes"] = "unassigned"

    train_adata.obs["subtypes"] = train_adata.obs["cluster"].astype(str)+'.'+train_adata.obs[celltype_cols].astype(str)
    test_adata.obs["subtypes"] = test_adata.obs["cluster"].astype(str)+'.'+test_adata.obs[celltype_cols].astype(str)

    return train_adata, test_adata

def process_PBMC_demultiplex(args, data_dir):
    '''Predict PBMC two-stage

    @args.data_source: PBMC_batch1, PBMC_batch2
    '''

    ### === PBMC demultiplex dataset
    from preprocess import process_PBMC_train_test
    from preprocess import load_PBMC_data
    if "PBMC_batch1" == args.data_source:
        ## extract batch1 info
        train_adata = process_PBMC_train_test.extract_batch1_batch(data_dir, batch = args.train)
        test_adata = process_PBMC_train_test.extract_batch1_batch(data_dir, batch = args.test)
        train_adata.obs.index.name = None
        test_adata.obs.index.name = None

    elif "PBMC_batch2" == args.data_source:
        ## extract batch2 info
        train_adata = process_PBMC_train_test.extract_batch2_condition(data_dir, cond = args.train)
        test_adata = process_PBMC_train_test.extract_batch2_condition(data_dir, cond = args.test)
        train_adata.obs.index.name = None
        test_adata.obs.index.name = None

    ## curate given sub-cell types to major cell types
    train_adata = process_PBMC_train_test.curate_PBMC_demulx_celltypes(train_adata, celltype_gran=0)
    test_adata = process_PBMC_train_test.curate_PBMC_demulx_celltypes(test_adata, celltype_gran=0)

    test_adata.obs["pred_sub_celltypes"] = "unassigned"
    return train_adata, test_adata

if __name__ == "__main__":
    data_dir = "/home/wma36/gpu/data"
    ## parse arguments
    import argparse
    parser = argparse.ArgumentParser(description="Celltyping pipeline.")
    parser.add_argument('data_source', help="Load which dataset",
        choices=["mousebrain_FC", "mousebrain_HC",
            "PBMC_batch1", "PBMC_batch2"])

    parser.add_argument('-m', '--method', help="Run which method",
        choices=['MLP', 'DFN', 'GEDFN', 'ItClust', 'SVM_RBF', 'SVM_linear', 'RF'],
        required=True)
    parser.add_argument('--select_on', help="Feature selection on train or test, or None of them",
        choices=['train', 'test'])
    parser.add_argument('--select_method', help="Feature selection method, Seurat/FEAST or None",
        choices=['Seurat', 'FEAST', "F-test"])
    parser.add_argument('--n_features', help="Number of features selected",
            default=1000, type=int)
    parser.add_argument('--train', help="Specify which as train", required=True)
    parser.add_argument('--test', help="Specify which as test", required=True)

    args = parser.parse_args()

    pipeline_dir = "/home/wma36/gpu/sc_identifier/pipelines"
    result_prefix = pipeline_dir+os.sep+"result_twostage_"+args.data_source+'_'+\
        args.train+'_to_'+args.test
    os.makedirs(result_prefix, exist_ok=True)

    if args.select_on is None and args.select_method is None:
        result_dir = result_prefix+os.sep+"no_feature"
    else:
        result_dir = result_prefix+os.sep+args.select_method+'_'+\
                str(args.n_features)+'_on_'+args.select_on
    os.makedirs(result_dir, exist_ok=True)

    from preprocess import process_train_test_data
    if "mousebrain_FC" == args.data_source or "mousebrain_HC" == args.data_source:
        train_adata, test_adata = process_mousebrain(args, data_dir)
    elif "PBMC_batch1" == args.data_source or "PBMC_batch2" == args.data_source:
        train_adata, test_adata = process_PBMC_demultiplex(args, data_dir)

    tmp_train_adata, tmp_test_adata = train_adata.copy(), test_adata.copy()  ## save a copy of data
    major_train_adata, major_test_adata = \
        process_train_test_data.process_pipeline(tmp_train_adata, tmp_test_adata, result_dir,
                gene_no=args.n_features, select_on=args.select_on, 
                select_method=args.select_method)
    
    ### --- run two-stage classifier
    print("\n\n=== %s \n" % args.method)
    start = time.time()
    # first stage prediction
    y_pred, enc = prediction(args.method, major_train_adata, major_test_adata)
    exec_time = time.time()-start

    celltype_cols = "cell.type"
    n_clusters = len(set(major_train_adata.obs[celltype_cols]))
    pred_onehot = np.zeros((y_pred.size, n_clusters))
    pred_onehot[np.arange(y_pred.size), y_pred] = 1
    pred_celltypes = enc.inverse_transform(pred_onehot)
    major_test_adata.obs['pred_major_celltypes'] = pred_celltypes
    for celltype in set(major_train_adata.obs[celltype_cols]):
        sub_train_cells = major_train_adata.obs.index[major_train_adata.obs["cell.type"] == celltype].tolist()
        sub_test_cells = major_test_adata.obs.index[major_test_adata.obs["pred_major_celltypes"] == celltype].tolist()
        if len(sub_test_cells) == 0: continue  # no need to cluster when the length of sub-cell types is 0

        sub_train_cells = [x.replace("-train", "") for x in sub_train_cells]
        sub_test_cells = [x.replace("-test", "") for x in sub_test_cells]

        sub_train_adata = train_adata[sub_train_cells]
        sub_train_adata.obs.rename(columns={celltype_cols: "majortypes", "subtypes": celltype_cols}, inplace=True)

        n_clusters = len(set(sub_train_adata.obs[celltype_cols]))
        if n_clusters == 1:
            test_adata.obs.loc[sub_test_cells, "pred_sub_celltypes"] = sub_train_adata.obs[celltype_cols][0]
            continue

        sub_test_adata = test_adata[sub_test_cells]
        sub_test_adata.obs.rename(columns={celltype_cols: "majortypes", "subtypes": celltype_cols}, inplace=True)

        tmp_sub_train_adata, tmp_sub_test_adata = sub_train_adata.copy(), sub_test_adata.copy()
        print("******Sub train data:", sub_train_adata.shape, "sub test data:", sub_test_adata.shape, "celltype:", celltype)
        sub_train_adata, sub_test_adata = \
                process_train_test_data.process_pipeline(sub_train_adata, sub_test_adata, 
                        result_dir, gene_no=args.n_features, select_on=args.select_on,
                        select_method=args.select_method, min_genes=0, min_cells=0)

        if sub_test_adata is None:  ## if after filtering, test adata becomes None, randomly assign cell types
            random.seed(RANDOM_SEED)
            sampled_celltypes = random.choices(list(set(tmp_sub_train_adata.obs[celltype_cols])), 
                    k=tmp_sub_test_adata.shape[0])
            test_adata.obs.loc[sub_test_cells, "pred_sub_celltypes"] = sampled_celltypes
            continue

        start = time.time()
        y_pred, enc = prediction(args.method, sub_train_adata, sub_test_adata)
        exec_time += (time.time()-start)

        pred_onehot = np.zeros((y_pred.size, n_clusters))
        pred_onehot[np.arange(y_pred.size), y_pred] = 1
        pred_celltypes = enc.inverse_transform(pred_onehot)
        sub_test_cells = [x.replace("-test", "") for x in sub_test_adata.obs_names.tolist()]
        test_adata.obs.loc[sub_test_cells, "pred_sub_celltypes"] = pred_celltypes

    print("\n\n=== Run time:", exec_time, "\n\n")
    test_adata.obs.to_csv(result_dir+os.sep+args.method+"_predicted_obs.csv")
    result_df = pd.read_csv(result_dir+os.sep+args.method+"_predicted_obs.csv")

    acc = metrics.accuracy_score(result_df["subtypes"], result_df["pred_sub_celltypes"])
    ARI = metrics.cluster.adjusted_rand_score(result_df["subtypes"], result_df["pred_sub_celltypes"])
    macroF1 = metrics.f1_score(result_df["subtypes"], result_df["pred_sub_celltypes"], average="macro")

    print("*****=====", "Testing accuracy: ", acc, " Testing ARI: ", ARI, 
            "Testing MacroF1: ", macroF1, "=====*****")
    ## write metrics result to file
    with open(result_dir+os.sep+args.method+"_metrics.txt", 'w') as f:
        f.write("Acc:%s\n" % str(acc))
        f.write("ARI:%s\n" % str(ARI))
        f.write("macroF1:%s\n" % str(macroF1))
        f.write("runtime:%s\n" % str(exec_time))


