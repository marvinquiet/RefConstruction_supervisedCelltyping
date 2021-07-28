'''
Configuration generation for running Pancreas datasets
'''

import os, argparse

from pipelines import method_utils, dataloading_utils
from preprocess.process_train_test_data import *

if __name__ == "__main__":
    data_dir = "~/gpu/data"

    ## parse arguments
    import argparse
    parser = argparse.ArgumentParser(description="Celltyping pipeline.")
    parser.add_argument('data_source', help="Load which dataset",
        choices=[
            'pancreas', 'pancreas_seg_cond', 'pancreas_custom', 
            'pancreas_seg_mix', 'pancreas_multi_to_multi'
            ])

    parser.add_argument('-m', '--method', help="Run which method",
        choices=['MLP', 'MLP_GO', 'MLP_CP', 'GEDFN', 'ItClust', 'SVM_RBF', 'SVM_linear', 'RF'], ## remove DFN
        required=True)
    parser.add_argument('--select_on', help="Feature selection on train or test, or None of them",
        choices=['train', 'test'])
    parser.add_argument('--select_method', help="Feature selection method, Seurat/FEAST or None",
            choices=['Seurat', 'FEAST', 'F-test'])
    parser.add_argument('--n_features', help="Number of features selected",
            default=1000, type=int)
    parser.add_argument('--train', help="Specify which as train", required=True)
    parser.add_argument('--test', help="Specify which as test", required=True)
    parser.add_argument('--sample_seed', help="Downsample seed in combined individual effect", 
            default=0, type=int)

    args = parser.parse_args()
    pipeline_dir = "pipelines/result_Pancreas_collections"
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
        train_adata, test_adata = dataloading_utils.load_Pancreas_adata(
            data_dir, result_dir, args=args)

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

    method_utils.run_pipeline(args, train_adata, test_adata, data_dir, result_dir)
