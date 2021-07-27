import anndata
import random

from preprocess import load_mousebrain_data
from preprocess import load_mouseFC_data
from preprocess import load_mouseprotocol_data
from preprocess import load_mouseAllen_data

# ==== process mouse data
def process_mousebrain_data(data_dir, result_dir, ind1, ind2, region="FC",
        celltype_gran=0, curate=False):
    '''Process Frontal Cortex from one ind to another
    @ind: P60FCAldh1l1Rep1/P60FCCx3cr1Rep1/P60FCRep1/P60FCRep2/P60FCRep3/P60FCRep4/P60FCRep6
    @celltype_gran: granularity of cell types. 0: major cell types, 1: sub-celltypes
    @purify: whether to purify the training dataset
    @purify_method: purify using distance or SVM
    @curate: whether to curate mouse brain data by combining different neurons/interneurons
    '''
    train_adata = load_mousebrain_data.load_brain_adata(data_dir, region=region,
            ind=ind1, celltype_gran=celltype_gran, curate=curate)
    test_adata = load_mousebrain_data.load_brain_adata(data_dir, region=region,
            ind=ind2, celltype_gran=celltype_gran, curate=curate)
    return train_adata, test_adata

def process_mousebrain_multiinds_data(data_dir, result_dir, sample_ind, pred_ind, 
        region="FC", sample=False, sample_seed=0, celltype_gran=0):
    '''Process mousebrain data from multiple individuals to another

    @sample_ind: downsample to a certain individual number
    @pred_ind: the predicted individual
    '''
    ## individuals for brain regions
    FC_inds = ["P60FCAldh1l1Rep1", "P60FCCx3cr1Rep1", "P60FCRep1", "P60FCRep2",
            "P60FCRep3", "P60FCRep4", "P60FCRep6"]
    HC_inds = ["P60HippoRep1", "P60HippoRep2", "P60HippoRep3",
            "P60HippoRep4", "P60HippoRep5", "P60HippoRep6"]

    if region == "FC":
        exclude_list = FC_inds
    elif region == "HC":
        exclude_list = HC_inds
    exclude_list.remove(pred_ind)
    exclude_inds = '_'.join(exclude_list)
    train_adata = load_mousebrain_data.load_brain_adata(data_dir, region=region,
            ind=exclude_inds, celltype_gran=celltype_gran)

    ## sample to number of cells when compared to the sample_id
    if sample:
        avg_number = train_adata.shape[0] // len(exclude_list)

        print("=== Downsample to number:", avg_number)
        random.seed(sample_seed)
        sampled_cells = random.sample(list(train_adata.obs_names), k=avg_number)
        train_adata = train_adata[sampled_cells]

    test_adata = load_mousebrain_data.load_brain_adata(data_dir, region=region,
            ind=pred_ind, celltype_gran=celltype_gran)
    return train_adata, test_adata

def process_mousebrain_regions(data_dir, result_dir, region1, region2,
        celltype_gran=0):
    ''' Use one region to predict another
    can allow region_individual input

    @celltype_gran: with sub-cell types and check ARI
    '''
    ## get region information and individuals
    regions1 = region1.split('_')
    regions2 = region2.split('_')

    region1 = regions1[0]
    if len(regions1) > 1:
        regions1_inds = '_'.join(regions1[1:])
    else:
        regions1_inds = None

    region2 = regions2[0]
    if len(regions2) > 1:
        regions2_inds = '_'.join(regions2[1:])
    else:
        regions2_inds = None

    ## get train/test anndata
    train_adata = load_mousebrain_data.load_brain_adata(data_dir, region=region1,
            ind=regions1_inds, celltype_gran=celltype_gran, curate=True)
    test_adata = load_mousebrain_data.load_brain_adata(data_dir, region=region2,
            ind=regions2_inds, celltype_gran=celltype_gran, curate=True)
    return train_adata, test_adata


## === process mouse FC data
def process_mousebrain_FC_stage(data_dir, result_dir, stage1, stage2,
        celltype_gran=0):
    '''Use one stage to predict another stage
    '''
    train_adata = load_mouseFC_data.load_FC_adata(data_dir, devstage=stage1,
            celltype_gran=celltype_gran)
    test_adata = load_mouseFC_data.load_FC_adata(data_dir, devstage=stage2,
            celltype_gran=celltype_gran)
    return train_adata, test_adata

def process_mousebrain_FC_datasets(data_dir, result_dir, dataset1_input, dataset2_input, 
        sample=False, sample_seed=0, celltype_gran=0):
    '''Use FC from one dataset to predict another

    @dataset1_input: "FC_inds", if sample == True, 
        the inds indicate the number of cells to downsample
    @dataset2_input: "Adult_inds"
    the above two parameters can be interchangable
    '''
    dataset1_split = dataset1_input.split('_')
    dataset2_split = dataset2_input.split('_')

    ## identify where the dataset points to
    regions = ["FC", "HC", "CB", "STR", "SN", "PC"]
    stages = ["Adult", "P21"]

    dataset1_input = dataset1_split[0]
    if len(dataset1_split) > 1:
        dataset1_inds = dataset1_split[1:]
        dataset1_inds = '_'.join(dataset1_inds)
    else:
        dataset1_inds = None

    dataset2_input = dataset2_split[0]
    if len(dataset2_split) > 1:
        dataset2_inds = dataset2_split[1:]
        dataset2_inds = '_'.join(dataset2_inds)
    else:
        dataset2_inds = None

    if dataset1_input in regions:
        train_adata = load_mousebrain_data.load_brain_adata(data_dir, 
                region=dataset1_input, ind=dataset1_inds, curate=True)
    else:
        info = dataset1_input.split('-')
        if info[0] in stages:
            if len(info) > 1:
                train_adata = load_mouseFC_data.load_FC_adata(data_dir, devstage=info[0],
                    treatment=info[1], curate=True)
            else:
                train_adata = load_mouseFC_data.load_FC_adata(data_dir, devstage=info[0],
                    ind=dataset1_inds, curate=True)

    if dataset2_input in regions:
        test_adata = load_mousebrain_data.load_brain_adata(data_dir,
                region=dataset2_input, ind=dataset2_inds, curate=True)
    else:
        info = dataset2_input.split('-')
        if info[0] in stages:
            if len(info) > 1:
                test_adata = load_mouseFC_data.load_FC_adata(data_dir, devstage=info[0],
                    treatment=info[1], curate=True)
            else:
                test_adata = load_mouseFC_data.load_FC_adata(data_dir, devstage=info[0],
                    ind=dataset2_inds, curate=True)
        
    if sample:
        ## sample number to the given individual cell numbers
        avg_number = train_adata.shape[0]
        if dataset1_input == "FC":
            train_adata = load_mousebrain_data.load_brain_adata(data_dir,
                    region=dataset1_input, ind=None, curate=True)
        
        if dataset1_input == "Adult":
            train_adata = load_mouseFC_data.load_FC_adata(data_dir, dataset1_input,
                    ind=None, curate=True)

        random.seed(sample_seed)
        sampled_cells = random.sample(list(train_adata.obs_names), k=avg_number)
        train_adata = train_adata[sampled_cells]
    return train_adata, test_adata

## Process Mouse cortex protocols datasets
def process_mousecortex_protocols_mouseFC_data(data_dir, result_dir, 
        input1, input2, curate=False):
    '''process both mousecortex protocols and mouse FC data to predict

    Note here we can only do major cell types because Cortex data does not 
    provide sub-cell types

    @input1/input2: MouseProtocol+exp+protocols_inds, mouseFC_inds
    @curate: whether to curate the datasets
    '''
    ## load train and test adata
    input1_list = input1.split('_')
    input1_dataset = input1_list[0]
    if len(input1_list) > 1:
        input1_inds = '_'.join(input1_list[1:])
    else:
        input1_inds = None

    input2_list = input2.split('_')
    input2_dataset = input2_list[0]
    if len(input2_list) > 1:
        input2_inds = '_'.join(input2_list[1:])
    else:
        input2_inds = None

    ## extract train and test adata accordingly
    if "mouseFC" in input1_dataset:
        train_adata = load_mousebrain_data.load_brain_adata(data_dir, region="FC",
                ind=input1_inds, celltype_gran=0, curate=True)
    elif "MouseProtocol" in input1_dataset:
        info = input1_dataset.split('+')
        exp = info[1]
        protocol = info[2]

        if exp == "Both":
            train_adata = load_mouseprotocol_data.load_mouseprotocol_adata(data_dir,
                    exp=None, protocol=protocol, curate=True)
        else:
            train_adata = load_mouseprotocol_data.load_mouseprotocol_adata(data_dir,
                    exp=exp, protocol=protocol, curate=True)

    if "mouseFC" in input2_dataset:
        test_adata = load_mousebrain_data.load_brain_adata(data_dir, region="FC",
                ind=input2_inds, celltype_gran=0, curate=True)
    elif "MouseProtocol" in input2_dataset:
        info = input2_dataset.split('+')
        exp = info[1]
        protocol = info[2]

        if exp == "Both":
            test_adata = load_mouseprotocol_data.load_mouseprotocol_adata(data_dir,
                    exp=None, protocol=protocol, curate=True)
        else:
            test_adata = load_mouseprotocol_data.load_mouseprotocol_adata(data_dir,
                    exp=exp, protocol=protocol, curate=True)
    return train_adata, test_adata

def process_FCdatasets_mouseFC_data(data_dir, result_dir, datasets_ref, mousebrain_tgt,
        celltype_gran=0, curate=False):
    '''Combine mouse brain protocol + pFC datasets together to predict the 
    mouse FC target individual

    @datasets_ref: MouseProtocol+Both+DroNc-seq_pFC+Adult+Saline
    @mousebrain_tgt: mousebrain dataset with certain targets
    '''
    datasets_list = datasets_ref.split('_')
    adata_list = []
    for dataset in datasets_list:
        infos = dataset.split('+')
        dataset_name = infos[0]
        if dataset_name == "MouseProtocol":
            exp = infos[1]
            protocol = infos[2]
            if exp == "Both":
                adata = load_mouseprotocol_data.load_mouseprotocol_adata(data_dir,
                        exp=None, protocol=protocol, curate=curate)
            else:
                adata = load_mouseprotocol_data.load_mouseprotocol_adata(data_dur,
                        exp=exp, protocol=protocol, curate=curate)
        elif dataset_name == "pFC":
            devstage = infos[1]
            treatment = infos[2]
            adata = load_mouseFC_data.load_FC_adata(data_dir, devstage=devstage,
                    treatment=treatment, celltype_gran=celltype_gran, curate=curate)
        adata_list.append(adata)

    ## find common genes and concatenate adata
    common_columns = ['barcode', 'nGene', 'nUMI', 'percent.mito', 'cell.type']
    for adata in adata_list:
        adata_obs = adata.obs[common_columns]
        adata.obs = adata_obs
    train_adata = anndata.AnnData.concatenate(*adata_list, join='inner')
    del adata_list  ## release space

    ## get target individuals list
    tgt_list = mousebrain_tgt.split('_')
    region = tgt_list[0]
    if len(tgt_list) > 1:
        tgt_inds = '_'.join(tgt_list[1:])
    else:
        tgt_inds = None

    test_adata = load_mousebrain_data.load_brain_adata(data_dir, region=region,
            ind=tgt_inds, celltype_gran=celltype_gran, curate=curate)
    return train_adata, test_adata


# --- add Allen Brain dataset
def process_allenbrain_SS_data(data_dir, result_dir, train_pct=0.8,
        sample_seed=0, celltype_gran=0, curate=True):
    '''Randomly select train/test from Allen Brain Smart-Seqv4 to do cross-validation

    @ train_pct: percentage selected for training
    '''
    adata = load_mouseAllen_data.load_mouseAllen_SS_data(data_dir, curate=curate)
    random.seed(sample_seed)
    train_cells = random.sample(adata.obs_names.tolist(), round(train_pct*adata.shape[0]))
    train_adata = adata[train_cells]
    test_cells = list(set(adata.obs_names) - set(train_cells))
    test_adata = adata[test_cells]
    return train_adata, test_adata


def process_allenbrain_10x_data(data_dir, result_dir, ind1, ind2, 
        celltype_gran=0, curate=True):
    '''Predict using Allen Brain dataset from one individual to another
    @ind: 371230 (M, P55), 371232 (F, P56), 372312 (M, P52), 372314 (M, P53), 381296 (F, P55)
    @celltype_gran: granularity of cell types. 0: major cell types, 1: sub-celltypes
    @curate: whether to curate mouse brain data by combining different neurons/interneurons
    '''
    ## in case to read data multiple times
    adata = load_mouseAllen_data.load_mouseAllen_10X_data(data_dir, curate=curate)
    train_cells = adata.obs[adata.obs["external_donor_name_label"].isin(ind1.split('_'))].index
    train_adata = adata[train_cells]
    test_cells = adata.obs[adata.obs["external_donor_name_label"].isin(ind2.split('_'))].index
    test_adata = adata[test_cells]
    return train_adata, test_adata

def process_allenbrain_cross_data(data_dir, result_dir, dataset1, dataset2="Allen",
        sample_seed=0, celltype_gran=0, curate=True):
    '''Predict allenbrain FACS-sorted SSv4 data 20% test cells

    @ dataset1: individual FC_P60FCRep1 or pFC_PFCSample1
    @ dataset2: Allen
    @ curate: whether to curate mouse brain cell types
    '''
    ## 20% test cells from allen brain dataset
    adata = load_mouseAllen_data.load_mouseAllen_SS_data(data_dir, curate=curate)
    random.seed(sample_seed)
    train_cells = random.sample(adata.obs_names.tolist(), round(0.8*adata.shape[0]))
    test_cells = list(set(adata.obs_names) - set(train_cells))
    test_adata = adata[test_cells]

    ## construct training dataset
    if "FC" == dataset1:
        train_adata = load_mousebrain_data.load_brain_adata(data_dir, region="FC",
            ind="P60FCRep1", celltype_gran=celltype_gran, curate=curate)
    if "pFC" == dataset1:
        train_adata = load_mouseFC_data.load_FC_adata(data_dir, devstage="Adult",
            ind="PFCSample1", curate=curate)
    return train_adata, test_adata

def process_mousebrain_crossdataset_inds_data(data_dir, result_dir, ref_inds, tgt_inds='mFC_P60FCCx3cr1Rep1',
        celltype_gran=0, curate=True):
    '''Predict mouse FC target individual using pFC individuals and allenbrain

    @ ref_inds: pFC_inds+allen_inds
    @ target_ind: FC_P60FCCx3cr1Rep1
    '''
    adata_list = []
    dataset_inds_list = ref_inds.split('+')
    for dataset_ind in dataset_inds_list:
        if 'mFC' in dataset_ind:
            FC_inds = dataset_ind.replace('mFC_', '')
            FC_adata = load_mousebrain_data.load_brain_adata(data_dir, region="FC",
                    ind=FC_inds, curate=curate)
            adata_list.append(FC_adata)
        if 'pFC' in dataset_ind:
            pFC_inds = dataset_ind.replace('pFC_', '')
            pFC_adata = load_mouseFC_data.load_FC_adata(data_dir, devstage='Adult', ind=pFC_inds, curate=curate)
            adata_list.append(pFC_adata)
        if 'allen' in dataset_ind:
            allen_inds = dataset_ind.replace('allen_', '')
            allen_adata = load_mouseAllen_data.load_mouseAllen_10X_data(data_dir, ind=allen_inds, curate=curate)
            adata_list.append(allen_adata)
    train_adata = anndata.AnnData.concatenate(*adata_list, join='inner')
    del adata_list ## release memory

    ## same as reference, for simplification, I just copied and pasted here
    adata_list = []
    dataset_inds_list = tgt_inds.split('+')
    for dataset_ind in dataset_inds_list:
        if 'mFC' in dataset_ind:
            FC_inds = dataset_ind.replace('mFC_', '')
            FC_adata = load_mousebrain_data.load_brain_adata(data_dir, region="FC",
                    ind=FC_inds, curate=curate)
            adata_list.append(FC_adata)
        if 'pFC' in dataset_ind:
            pFC_inds = dataset_ind.replace('pFC_', '')
            pFC_adata = load_mouseFC_data.load_FC_adata(data_dir, devstage='Adult', ind=pFC_inds, curate=curate)
            adata_list.append(pFC_adata)
        if 'allen' in dataset_ind:
            allen_inds = dataset_ind.replace('allen_', '')
            allen_adata = load_mouseAllen_data.load_mouseAllen_10X_data(data_dir, ind=allen_inds, curate=curate)
            adata_list.append(allen_adata)
    test_adata = anndata.AnnData.concatenate(*adata_list, join='inner')
    del adata_list ## release memory

    return train_adata, test_adata
