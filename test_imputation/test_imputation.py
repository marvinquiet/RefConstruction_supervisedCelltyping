'''
Run MAGIC/scVI/SAVER to impute gene expression matrix

When running scVI, envs/scvi needs to be activated first
'''
import os, sys, time, argparse

from pipelines import method_utils, dataloading_utils
from preprocess.process_train_test_data import *
import scanpy.api as sc
import pandas as pd

def run_MAGIC(train_adata, test_adata):
    import magic

    train_magic_op = magic.MAGIC()
    if scipy.sparse.issparse(train_adata.X):
        x_train = train_adata.X.toarray()
    else:
        x_train = train_adata.X
    train_emt_magic = train_magic_op.fit_transform(x_train, genes='all_genes')
    train_adata.X = train_emt_magic
    ## standardize the input
    sc.pp.scale(train_adata, zero_center=True, max_value=6)

    test_magic_op = magic.MAGIC()
    if scipy.sparse.issparse(test_adata.X):
        x_test = test_adata.X.toarray()
    else:
        x_test = test_adata.X
    test_emt_magic = test_magic_op.fit_transform(x_test, genes='all_genes')
    test_adata.X = test_emt_magic
    ## standardize the input
    sc.pp.scale(test_adata, zero_center=True, max_value=6)

    return train_adata, test_adata

def run_scVI(result_dir):
    '''Run scVI by activating scVI environment because I do not want to mix up torch and tensorflow

    ## save processed data to file; activate scVI environment; run scVI script; save scVI data to file; deactivate environment; re-activate celltyping environment
    '''

    #TODO: this place is dirty, need to clean up later...
    scVI_envs_dir = "/homelocal/wma36/envs/scvi"
    scVI_script_dir = "/homelocal/wma36/celltyping_refConstruct/test_imputation/test_scVI.py"
    conda_dir = "/home/wma36/miniconda3/etc/profile.d/conda.sh"

    import subprocess
    ## pyton -V is used to check python version
    subprocess.run('bash -c ". %s; conda activate %s; python %s %s"' % 
            (conda_dir, scVI_envs_dir, scVI_script_dir, result_dir), shell=True)

    train_adata = anndata.read_h5ad(result_dir+os.sep+'scVI_train_adata.h5ad')
    test_adata = anndata.read_h5ad(result_dir+os.sep+'scVI_test_adata.h5ad')
    return train_adata, test_adata

def run_SAVER(train_adata, test_adata, result_dir):
    '''Run SAVER through Rscript and store the data
    '''
    SAVER_script_dir = "/homelocal/wma36/celltyping_refConstruct/test_imputation/test_SAVER.R"

    x_train = None
    if scipy.sparse.issparse(train_adata.X):
        x_train = train_adata.X.toarray()
    else:
        x_train = train_adata.X

    x_test = None
    if scipy.sparse.issparse(test_adata.X):
        x_test = test_adata.X.toarray()
    else:
        x_test = test_adata.X

    train_df = pd.DataFrame(data=x_train, 
            index=train_adata.obs_names.tolist(), columns=train_adata.var_names.tolist())
    train_df.to_csv(result_dir+os.sep+'train_data.csv')
    test_df = pd.DataFrame(data=x_test,
            index=test_adata.obs_names.tolist(), columns=test_adata.var_names.tolist())
    test_df.to_csv(result_dir+os.sep+'test_data.csv')

    import subprocess
    subprocess.run('Rscript --vanilla %s %s' % (SAVER_script_dir, result_dir), shell=True)

    x_train = pd.read_csv(result_dir+os.sep+'SAVER_train_data.csv', index_col=0, header=0)
    x_test = pd.read_csv(result_dir+os.sep+'SAVER_test_data.csv', index_col=0, header=0)
    subprocess.run("rm %s" % (result_dir+os.sep+'train_data.csv'), shell=True)
    subprocess.run("rm %s" % (result_dir+os.sep+'test_data.csv'), shell=True)

    ## set data to anndata and then scale
    train_adata.X = np.array(x_train.loc[train_adata.obs_names.tolist(), train_adata.var_names.tolist()])
    test_adata.X = np.array(x_test.loc[test_adata.obs_names.tolist(), test_adata.var_names.tolist()])
    sc.pp.scale(train_adata, zero_center=True, max_value=6)
    sc.pp.scale(test_adata, zero_center=True, max_value=6)

    subprocess.run("rm %s" % (result_dir+os.sep+'SAVER_train_data.csv'), shell=True)
    subprocess.run("rm %s" % (result_dir+os.sep+'SAVER_test_data.csv'), shell=True)
    return train_adata, test_adata


if __name__ == '__main__':
    data_dir = "/home/wma36/gpu/data"

    ## parse arguments
    parser = argparse.ArgumentParser(description="imputation pipeline.")
    parser.add_argument('data_source', help="Load which dataset")
    parser.add_argument('--imputation', help="Imputation method", choices=['MAGIC', 'scVI', 'SAVER'])
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

    pipeline_dir = "/home/wma36/gpu/celltyping_refConstruct/pipelines/result_imputation_collections"
    result_prefix = pipeline_dir+os.sep+args.imputation+os.sep+"result_"+args.data_source+'_'+\
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

        if args.imputation == "MAGIC":
            train_adata, test_adata = dataloading_utils.process_loaded_data(
                train_adata, test_adata, result_dir, args=args, purify_method=purify_method, scale=False)
        if args.imputation == "scVI":
            train_adata, test_adata = dataloading_utils.process_loaded_data(
                train_adata, test_adata, result_dir, args=args, purify_method=purify_method, scale=False, save_raw=True)
        if args.imputation == "SAVER":
            train_adata, test_adata = dataloading_utils.process_loaded_data(
                train_adata, test_adata, result_dir, args=args, purify_method=purify_method, scale=False)
        print("Train anndata: \n", train_adata)
        print("Test anndata: \n", test_adata)

    if args.imputation == "MAGIC":
        train_adata, test_adata = run_MAGIC(train_adata.copy(), test_adata.copy())
    if args.imputation == "scVI":
        train_adata, test_adata = run_scVI(result_dir)
    if args.imputation == "SAVER":
        train_adata, test_adata = run_SAVER(train_adata.copy(), test_adata.copy(), result_dir)

    method_utils.run_pipeline(args, train_adata, test_adata, data_dir, result_dir)
