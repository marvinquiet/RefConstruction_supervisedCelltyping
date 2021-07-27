'''
Configuration generation for running PBMC datasets
'''

import os, argparse

from pipelines import method_utils, dataloading_utils
from preprocess.process_train_test_data import *

if __name__ == "__main__":
    data_dir = "/home/wma36/gpu/data"

    ## parse arguments
    parser = argparse.ArgumentParser(description="Celltyping pipeline.")
    parser.add_argument('data_source', help="Load which dataset",
        choices=[
            'PBMC_batch1_ind', 'PBMC_batch1_ind_major',
            'PBMC_batch1_ind_purify_dist', 'PBMC_batch1_ind_purify_SVM',
            'PBMC_batch1_ABC', "PBMC_batch1_ABC_purify_dist", "PBMC_batch1_ABC_purify_SVM",
            'PBMC_batch1_batchtoind', 'PBMC_batch1_multiinds', 'PBMC_batch1_multiinds_sample',
            'PBMC_batch2', "PBMC_batch2_purify_dist", "PBMC_batch2_purify_SVM",
            'PBMC_batch1_batch2_ind', 'PBMC_batch1_batch2_ind_major',
            'PBMC_protocols_pbmc1', 'PBMC_protocols_batch_smart',
            'PBMC_Zheng_FACS', 'PBMC_Zheng_FACS_curated',
            'PBMC_cross'  ## cross PBMC datasets
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
    if args.data_source in ["PBMC_batch1_ind", "PBMC_batch1_ABC", "PBMC_batch2", 
            "PBMC_batch1_batchtoind", "PBMC_batch1_multiinds"
            ]:
        pipeline_dir = "/home/wma36/gpu/celltyping_refConstruct/pipelines/result_PBMC_collections"

    if args.data_source in ["PBMC_protocols_pbmc1", "PBMC_protocols_batch_smart"]:
        pipeline_dir = "/home/wma36/gpu/celltyping_refConstruct/pipelines/result_PBMC_protocols_collections"

    if args.data_source in ["PBMC_Zheng_FACS", "PBMC_Zheng_FACS_curated", "PBMC_cross"]:
        pipeline_dir = "/home/wma36/gpu/celltyping_refConstruct/pipelines/result_PBMC_Zheng_collections"

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
