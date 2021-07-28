'''
Configuration generation for running performance saturation
'''
import os, sys, argparse
import random

from pipelines import method_utils, dataloading_utils
from preprocess.process_train_test_data import *

if __name__ == "__main__":
    data_dir = "~/gpu/data"

    ## parse arguments
    parser = argparse.ArgumentParser(description="Test for saturation pipeline.")
    parser.add_argument('data_source', help="Load which dataset")
    parser.add_argument('--train', help="Specify which as train", required=True)
    parser.add_argument('--test', help="Specify which as test", required=True)
    parser.add_argument('--sample_seed', help="Downsample seed in combined individual effect", 
            default=0, type=int)
    parser.add_argument('--downsample', help="Whether do downsample or not",
            dest='downsample', action='store_true')
    parser.add_argument('--downsample_size', help="Downsample size for testing saturation",
            default=3000, type=int)

    args = parser.parse_args()
    args.method = "MLP"
    args.select_on = "train"
    args.select_method = "F-test"
    args.n_features = 1000

    if args.downsample:
        pipeline_dir = "pipelines/result_saturation_downsample_collections"
    else:
        pipeline_dir = "pipelines/result_saturation_collections"
    result_collections_dir = pipeline_dir+os.sep+"result_"+args.data_source+'_'+args.train.split('_')[0]+'s'
    result_prefix = result_collections_dir+os.sep+str(args.sample_seed) ## to distinguish PBMC tasks
    os.makedirs(result_prefix, exist_ok=True)

    ## create file directory 
    if args.select_on is None and args.select_method is None:
        result_dir = result_prefix+os.sep+"no_feature"
    else:
        result_dir = result_prefix+os.sep+args.select_method+'_'+\
                str(args.n_features)+'_on_'+args.select_on
    os.makedirs(result_dir, exist_ok=True)

    load_ind, train_adata, test_adata = load_adata(result_collections_dir)
    if not load_ind:
        train_adata, test_adata = dataloading_utils.load_PBMC_adata(
            data_dir, result_collections_dir, args=args)
        if train_adata is None or test_adata is None: 
            train_adata, test_adata = dataloading_utils.load_Pancreas_adata(
                data_dir, result_collections_dir, args=args)
        if train_adata is None or test_adata is None:
            train_adata, test_adata = dataloading_utils.load_Mousebrain_adata(
                data_dir, result_collections_dir, args=args)
        if train_adata is None or test_adata is None:
            sys.exit("Please check your data source to make sure it matches an input.")

        ## whether to purify reference dataset
        purify_method = ""
        if "purify_dist" in args.data_source:
            purify_method = "distance"
        elif "purify_SVM" in args.data_source:
            purify_method = "SVM"

        train_adata, test_adata = dataloading_utils.process_loaded_data(
            train_adata, test_adata, result_collections_dir, args=args, purify_method=purify_method)
        print("Train anndata: \n", train_adata)
        print("Test anndata: \n", test_adata)

    if args.downsample: ## add shuffled cells
        train_cells = train_adata.obs_names.tolist()
        random.seed(args.sample_seed)
        random.shuffle(train_cells) ## shuffle original cell list
        original_train_adata = train_adata.copy()[train_cells]
        for i in range(original_train_adata.shape[0]//args.downsample_size+1):
            sampled_number = (i+1)*args.downsample_size if (i+1)*args.downsample_size < original_train_adata.shape[0] else original_train_adata.shape[0] 
            train_adata = original_train_adata[:sampled_number]
            sampled_result_dir = result_dir+os.sep+str(sampled_number)
            os.makedirs(sampled_result_dir, exist_ok=True)
            method_utils.run_pipeline(args, train_adata, test_adata, data_dir, sampled_result_dir)
    else: ## add shuffled individuals
        if args.data_source == "mousebrain_crossdataset_inds": ## a combined version
            pFC_samples = [x for x in train_adata.obs["Sample"].tolist() if x != 'nan']
            allen_samples = [x for x in train_adata.obs["external_donor_name_label"].tolist() if x != 'nan']
            train_adata.obs["ind"] =  pFC_samples + allen_samples

        original_train_adata = train_adata.copy()
        train_inds = list(set(original_train_adata.obs['ind']))
        random.seed(args.sample_seed)
        random.shuffle(train_inds)
        for i in range(len(train_inds)):
            train_idx = original_train_adata.obs['ind'].isin(train_inds[:i+1])
            train_adata = original_train_adata[train_idx]
            ind_result_dir = result_dir+os.sep+str(i+1)
            os.makedirs(ind_result_dir, exist_ok=True)
            method_utils.run_pipeline(args, train_adata, test_adata, data_dir, ind_result_dir)
