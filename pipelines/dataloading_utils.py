'''
Load data for running celltyping experiments
'''
import os, sys
from preprocess.process_train_test_data import *

from preprocess.process_PBMC_train_test import *
from preprocess.process_Pancreas_train_test import *
from preprocess.process_mousebrain_train_test import *

### === process loaded data
def process_loaded_data(train_adata, test_adata, result_dir, 
        args=None, scale=True, plot=True, purify_method="", 
        save_raw=False, save_data=True):

    if args is None:
        sys.exit("Error: Please check your argument parser object!")

    # curate for common cell types
    common_celltypes = set(train_adata.obs["cell.type"]).intersection(set(test_adata.obs["cell.type"]))
    train_cells = train_adata.obs.loc[train_adata.obs["cell.type"].isin(common_celltypes)].index
    test_cells = test_adata.obs.loc[test_adata.obs["cell.type"].isin(common_celltypes)].index
    train_adata = train_adata[train_cells]
    test_adata = test_adata[test_cells]

    if save_raw:
        train_adata.layers["counts"] = train_adata.X.copy() 
        test_adata.layers["counts"] = test_adata.X.copy() 

    ## feature selection
    train_adata, test_adata = feature_selection_train_test(train_adata, test_adata,
            result_dir, args.n_features, args.select_on, args.select_method)
    ## scale and analze
    train_adata, test_adata = scale_and_visualize(train_adata, test_adata,
        result_dir, purify_method=purify_method, scale=scale, plot=plot)
    if save_data:
        save_adata(train_adata, test_adata, result_dir)

    return train_adata, test_adata

### === load PBMC anndata
def load_PBMC_adata(data_dir, result_dir, args=None, scale=True, plot=True):
    train_adata, test_adata = None, None
    ## === PBMC datasets
    if args.data_source == "PBMC_batch1_ind":
        train_adata, test_adata = \
            process_batch1_ind(data_dir, result_dir, 
                ind1=args.train, ind2=args.test)
    if args.data_source == "PBMC_batch1_ind_major":
        train_adata, test_adata = \
            process_batch1_ind(data_dir, result_dir, 
                ind1=args.train, ind2=args.test,
                celltype_gran=0)
    if args.data_source == "PBMC_batch1_ind_purify_dist":
        train_adata, test_adata = \
            process_batch1_ind(data_dir, result_dir, 
                ind1=args.train, ind2=args.test)
    if args.data_source == "PBMC_batch1_ind_purify_SVM":
        train_adata, test_adata = \
            process_batch1_ind(data_dir, result_dir, 
                ind1=args.train, ind2=args.test)
    if args.data_source == "PBMC_batch1_ABC":
        train_adata, test_adata = \
            process_batch1_ABC(data_dir, result_dir, 
                batch1=args.train, batch2=args.test)
    if args.data_source == "PBMC_batch1_ABC_purify_dist":
        train_adata, test_adata = \
            process_batch1_ABC(data_dir, result_dir, 
                batch1=args.train, batch2=args.test)
    if args.data_source == "PBMC_batch1_ABC_purify_SVM":
        train_adata, test_adata = \
            process_batch1_ABC(data_dir, result_dir, 
                batch1=args.train, batch2=args.test)
    if args.data_source == "PBMC_batch1_batchtoind":
        train_adata, test_adata = \
                process_batch1_batchtoind(data_dir, result_dir,
                        batch=args.train, ind=args.test, sample=True, sample_seed=0)
    if args.data_source == "PBMC_batch1_multiinds":
        train_adata, test_adata = \
                process_batch1_multiinds(data_dir, result_dir,
                        sample_ind=args.train, pred_ind=args.test,
                        sample=False)
    if args.data_source == "PBMC_batch1_multiinds_sample":
        train_adata, test_adata = \
                process_batch1_multiinds(data_dir, result_dir,
                        sample_ind=args.train, pred_ind=args.test,
                        sample=True, sample_seed=args.sample_seed)
    if args.data_source == "PBMC_batch2":
        train_adata, test_adata = \
            process_batch2_ctrl_stim(data_dir, result_dir, 
                cond1=args.train, cond2=args.test)
    if args.data_source == "PBMC_batch2_purify_dist":
        train_adata, test_adata = \
            process_batch2_ctrl_stim(data_dir, result_dir, 
                cond1=args.train, cond2=args.test)
    if args.data_source == "PBMC_batch2_purify_SVM":
        train_adata, test_adata = \
            process_batch2_ctrl_stim(data_dir, result_dir, 
                cond1=args.train, cond2=args.test)
    if args.data_source == "PBMC_batch1_batch2_ind":
        train_adata, test_adata = \
            process_batch1_batch2_ind(data_dir, result_dir,
                    input1=args.train, input2=args.test)
    if args.data_source == "PBMC_batch1_batch2_ind_major":
        train_adata, test_adata = \
             process_batch1_batch2_ind(data_dir, result_dir,
                    input1=args.train, input2=args.test,
                    celltype_gran=0)

    ## === PBMC 7 protocols datasets
    if args.data_source == "PBMC_protocols_pbmc1":
        train_adata, test_adata = \
            process_PBMC_protocols_type(data_dir, result_dir, 
                protocols1=args.train, protocols2=args.test, exp="pbmc1")
    if args.data_source == "PBMC_protocols_batch_smart":
        train_adata, test_adata = \
            process_PBMC_protocols_batch(data_dir, result_dir, 
                batch1=args.train, batch2=args.test, protocol="Smart-seq2")

    ## === PBMC Zheng FACS-sorted dataset
    if args.data_source == "PBMC_Zheng_FACS":
        train_adata, test_adata = \
            process_PBMC_Zheng(data_dir, result_dir,
                train_pct=float(args.train), sample_seed=0)
    if args.data_source == "PBMC_Zheng_FACS_curated":
        train_adata, test_adata = \
            process_PBMC_Zheng(data_dir, result_dir,
                train_pct=float(args.train), sample_seed=0,
                curate=True)
    if args.data_source == "PBMC_cross":
        train_adata, test_adata = \
            process_PBMC_comparison(data_dir, result_dir,
                train=args.train, test=args.test, sample_seed=0)

    return train_adata, test_adata

### === load pancreas anndata
def load_Pancreas_adata(data_dir, result_dir, args=None, scale=True, plot=True):
    train_adata, test_adata = None, None
    ## === pancrease datasets
    if args.data_source == "pancreas":
        train_adata, test_adata = \
            process_pancreas(data_dir, result_dir, 
                dataset1=args.train, dataset2=args.test)
    if args.data_source == "pancreas_seg_cond":
        train_adata, test_adata = \
            process_pancreas_seg_cond(data_dir, result_dir, 
                cond1=args.train, cond2=args.test)
    if args.data_source == "pancreas_custom":
        train_adata, test_adata = \
            process_pancreas_custom(data_dir, result_dir, 
                seg_cond=args.train, dataset=args.test)
    if args.data_source == "pancreas_seg_mix":
        train_adata, test_adata = \
            process_pancreas_seg_mix(data_dir, result_dir, 
                main_cond=args.train, pred_cond=args.test)
    if args.data_source == "pancreas_multi_to_multi":
        train_adata, test_adata = \
            process_pancreas_multi_to_multi(data_dir, result_dir, 
                cond1=args.train, cond2=args.test)
    return train_adata, test_adata

### === load mousebrain anndata
def load_Mousebrain_adata(data_dir, result_dir, args=None, scale=True, plot=True):
    train_adata, test_adata = None, None
    ## === mouse adult brain
    if args.data_source == "mousebrain_FC":
        train_adata, test_adata = \
            process_mousebrain_data(data_dir, result_dir, 
                region="FC", ind1=args.train, ind2=args.test)
    if args.data_source == "mousebrain_FC_curate":
        train_adata, test_adata = \
            process_mousebrain_data(data_dir, result_dir, 
                region="FC", ind1=args.train, ind2=args.test,
                curate=True)
    if args.data_source == "mousebrain_FC_sub":
        train_adata, test_adata = \
            process_mousebrain_data(data_dir, result_dir, 
                region="FC", ind1=args.train, ind2=args.test,
                celltype_gran=1)
    if args.data_source == "mousebrain_FC_sub_purify_dist":
        train_adata, test_adata = \
            process_mousebrain_data(data_dir, result_dir, 
                region="FC", ind1=args.train, ind2=args.test,
                celltype_gran=1)
    if args.data_source == "mousebrain_FC_sub_purify_SVM":
        train_adata, test_adata = \
            process_mousebrain_data(data_dir, result_dir, 
                region="FC", ind1=args.train, ind2=args.test,
                celltype_gran=1)
    if args.data_source == "mousebrain_FC_multiinds":
        train_adata, test_adata = \
            process_mousebrain_multiinds_data(
                data_dir, result_dir, region="FC", 
                sample_ind=args.train, pred_ind=args.test)
    if args.data_source == "mousebrain_FC_multiinds_sub":
        train_adata, test_adata = \
            process_mousebrain_multiinds_data(
                data_dir, result_dir, region="FC", 
                sample_ind=args.train, pred_ind=args.test,
                celltype_gran=1)
    if args.data_source == "mousebrain_FC_multiinds_sample":
        train_adata, test_adata = \
            process_mousebrain_multiinds_data(
                data_dir, result_dir, region="FC", sample_seed=args.sample_seed,
                sample_ind=args.train, pred_ind=args.test,
                sample=True)
    if args.data_source == "mousebrain_FC_multiinds_sub_sample":
        train_adata, test_adata = \
            process_mousebrain_multiinds_data(
                data_dir, result_dir, region="FC", sample_seed=args.sample_seed,
                sample_ind=args.train, pred_ind=args.test,
                celltype_gran=1, sample=True)
    if args.data_source == "mousebrain_HC":
        train_adata, test_adata = \
            process_mousebrain_data(data_dir, result_dir, 
                region="HC", ind1=args.train, ind2=args.test)
    if args.data_source == "mousebrain_HC_sub":
        train_adata, test_adata = \
            process_mousebrain_data(data_dir, result_dir, 
                region="HC", ind1=args.train, ind2=args.test,
                celltype_gran=1)
    if args.data_source == "mousebrain_region":
        train_adata, test_adata = \
            process_mousebrain_regions(data_dir, result_dir, 
                region1=args.train, region2=args.test)
    if args.data_source == "mousebrain_region_sub":
        train_adata, test_adata = \
            process_mousebrain_regions(data_dir, result_dir, 
                region1=args.train, region2=args.test,
                celltype_gran=1)

    ## === mousebrain FC
    if args.data_source == "mousebrain_FC_stage":
        train_adata, test_adata = \
            process_mousebrain_FC_stage(data_dir, result_dir, 
                stage1=args.train, stage2=args.test)
    if args.data_source == "mousebrain_FC_stage_sub":
        ## sub-celltypes
        train_adata, test_adata = \
            process_mousebrain_FC_stage(data_dir, result_dir, 
                stage1=args.train, stage2=args.test,
                celltype_gran=1)
    if args.data_source == "mousebrain_FC_datasets":
        ## predict using dataset1 to dataset2
        train_adata, test_adata = \
            process_mousebrain_FC_datasets(data_dir, result_dir, 
                dataset1_input=args.train, dataset2_input=args.test)
    if args.data_source == "mousebrain_FC_datasets_multiinds":
        ## predict using dataset1 to dataset2
        train_adata, test_adata = \
            process_mousebrain_FC_datasets(data_dir, result_dir, 
                dataset1_input=args.train, dataset2_input=args.test)
    if args.data_source == "mousebrain_FC_datasets_multiinds_sample":
        ## predict using dataset1 to dataset2
        train_adata, test_adata = \
            process_mousebrain_FC_datasets(data_dir, result_dir, 
                dataset1_input=args.train, dataset2_input=args.test,
                sample=True, sample_seed=args.sample_seed)

    ## === Mouse cortex data
    if args.data_source == "mousecortex_protocols_mouseFC_ind":
        train_adata, test_adata = \
            process_mousecortex_protocols_mouseFC_data(data_dir, result_dir,
                input1=args.train, input2=args.test, curate=True)

    ## === Combine mouse brain data to predict FC
    if args.data_source == "mousebrain_combined_mouseFC":
        train_adata, test_adata = \
            process_FCdatasets_mouseFC_data(data_dir, result_dir,
                datasets_ref=args.train, mousebrain_tgt=args.test, curate=True)

    ## === allenbrain within dataset
    if args.data_source == "allenbrain_ss":
        train_adata, test_adata = \
            process_allenbrain_SS_data(data_dir, result_dir,
                train_pct=float(args.train), sample_seed=0, curate=True)
    if args.data_source == "allenbrain_10x":
        train_adata, test_adata = \
            process_allenbrain_10x_data(data_dir, result_dir,
                ind1=args.train, ind2=args.test, curate=True)
    ## === allenbrain cross dataset
    if args.data_source == "allenbrain_cross":
        train_adata, test_adata = \
            process_allenbrain_cross_data(data_dir, result_dir,
                dataset1=args.train, dataset2="Allen", sample_seed=0,
                curate=True)

    ## === mousebrain cross dataset -> FC target
    if args.data_source == "mousebrain_crossdataset_inds":
        train_adata, test_adata = \
            process_mousebrain_crossdataset_inds_data(data_dir, result_dir,
                ref_inds=args.train, tgt_inds=args.test, curate=True)

    return train_adata, test_adata
