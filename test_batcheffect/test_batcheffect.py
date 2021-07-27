'''
Run MAGIC/scVI to impute gene expression matrix

When running scVI, envs/scVI needs to be activated first
'''
import os, sys, time, argparse

from pipelines import method_utils, dataloading_utils
from preprocess.process_train_test_data import *
import scanpy.api as sc
import numpy as np
import pandas as pd

def run_batcheffectremoval(train_adata, test_adata, result_dir, method="Harmony"):
    '''Run batch effect removal methods to remove batch effect between train and test
    @ method: Harmony or fastMNN in R script
    '''
    # @TODO: may need to remove the absolute path
    batcheffect_dir = "/homelocal/wma36/celltyping_refConstruct/test_batcheffect/test_batcheffect.R"

    adata = train_adata.concatenate(test_adata, join='inner', index_unique=None)
    adata_count = None
    if scipy.sparse.issparse(adata.layers['counts']):
        adata_count = adata.layers['counts'].toarray()
    else:
        adata_count = adata.layers['counts']

    count_df = pd.DataFrame(adata_count,
            index=adata.obs_names.tolist(), columns=adata.var_names.tolist(), dtype="int32")
    count_df.to_csv(result_dir+os.sep+'adata_count.csv')
    adata.obs['batch'].to_csv(result_dir+os.sep+'batchID.csv')

    import subprocess
    subprocess.run('Rscript --vanilla %s %s %s' % (batcheffect_dir, method, result_dir), shell=True)
    subprocess.run("rm %s" % (result_dir+os.sep+'adata_count.csv'), shell=True)
    subprocess.run("rm %s" % (result_dir+os.sep+'batchID.csv'), shell=True)

    corrected_counts = pd.read_csv(result_dir+os.sep+'corrected_counts.csv', 
            index_col=0, header=0)
    corrected_adata = anndata.AnnData(X=corrected_counts)
    corrected_adata.obs = adata.obs
    corrected_adata.obsm = adata.obsm
    corrected_train_adata = corrected_adata[corrected_adata.obs['dataset_batch'] == "train"]
    corrected_test_adata = corrected_adata[corrected_adata.obs['dataset_batch'] == "test"]
 
    sc.pp.scale(corrected_train_adata, zero_center=True, max_value=6)
    sc.pp.scale(corrected_test_adata, zero_center=True, max_value=6)
    subprocess.run("rm %s" % (result_dir+os.sep+'corrected_counts.csv'), shell=True)
    return corrected_train_adata, corrected_test_adata

if __name__ == '__main__':
    data_dir = "/home/wma36/gpu/data"

    ## parse arguments
    parser = argparse.ArgumentParser(description="Batch effect removal pipeline.")
    parser.add_argument('data_source', help="Load which dataset")
    parser.add_argument('--batcheffect', help="Batch effect removal methods", choices=['Harmony', 'fastMNN'])
    parser.add_argument('--train', help="Specify which as train", required=True)
    parser.add_argument('--test', help="Specify which as test", required=True)
    parser.add_argument('--sample_seed', help="Downsample seed in combined individual effect", 
            default=0, type=int)

    args = parser.parse_args()
    ## F-test on train + MLP
    args.method = "MLP"
    args.select_on = "train"
    args.select_method = "F-test"
    args.n_features = 1000

    pipeline_dir = "/home/wma36/gpu/celltyping_refConstruct/pipelines/result_batchremoval_collections"
    result_prefix = pipeline_dir+os.sep+args.batcheffect+os.sep+"result_"+args.data_source+'_'+\
        args.train+'_to_'+args.test
    os.makedirs(result_prefix, exist_ok=True)
    result_dir = result_prefix+os.sep+args.select_method+'_'+str(args.n_features)+'_on_'+args.select_on
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
            train_adata, test_adata, result_dir, args=args, 
            purify_method=purify_method, scale=False, save_raw=True)
    train_adata, test_adata = run_batcheffectremoval(train_adata.copy(), test_adata.copy(), result_dir, args.batcheffect)

    method_utils.run_pipeline(args, train_adata, test_adata, data_dir, result_dir)
