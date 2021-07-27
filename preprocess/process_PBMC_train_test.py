import anndata
import random

from preprocess import load_PBMC_data
from preprocess import load_PBMCprotocol_data
from preprocess import load_PBMCZheng_data

def process_batch1_ind(data_dir, result_dir, ind1="1154", ind2="1085", 
        celltype_gran=1):
    ''' Process individual data of PBMC batch1

    @data_dir: where PBMC batch1 data stroes
    @result_dir: where to store PCA/tSNE/UMAP result
    @celltype_gran: granularity of cell types, 0: major cell types, 1:sub-celltypes
        According to the dataset, the give cell types are sub-cell types
    '''
    ## add batch info
    train_adata = load_PBMC_data.load_PBMC_batch1_data(data_dir, ind=ind1)
    test_adata = load_PBMC_data.load_PBMC_batch1_data(data_dir, ind=ind2)

    ## curate given sub-cell types to major cell types
    train_adata = load_PBMC_data.curate_PBMC_demulx_celltypes(train_adata, celltype_gran)
    test_adata = load_PBMC_data.curate_PBMC_demulx_celltypes(test_adata, celltype_gran)
    return train_adata, test_adata

def process_batch1_multiinds(data_dir, result_dir, sample_ind, pred_ind,
        celltype_gran=1, sample=False, sample_seed=0):
    '''Process multi individuals to predict one individual

    @sample_ind: downsample to a certain individual number
    @pred_ind: the predicted individual
    '''
    batch1_inds = [1043, 1079, 1154, 1249, 1493, 1511, 1598, 1085]

    exclude_list = batch1_inds
    exclude_list.remove(int(pred_ind))

    train_adata = load_PBMC_data.load_PBMC_batch1_data(data_dir, ind=exclude_list)
    if sample: ## sample to average number
        avg_number = train_adata.shape[0] // len(exclude_list)

        print("=== Downsample to number:", avg_number)
        random.seed(sample_seed)
        sampled_cells = random.sample(list(train_adata.obs_names), k=avg_number)
        train_adata = train_adata[sampled_cells]

    test_adata = load_PBMC_data.load_PBMC_batch1_data(data_dir, ind=pred_ind)
    train_adata.obs.index.name = None
    test_adata.obs.index.name = None
    return train_adata, test_adata


def process_batch1_ABC(data_dir, result_dir, batch1="A", batch2="B", 
        celltype_gran=1):
    ''' Process two batches for GEDFN

    @data_dir: where PBMC batch1/batch2 data stroes
    @result_dir: where to store PCA/tSNE/UMAP result
    @celltype_gran: default 1, given as sub-cell types
    '''
    ## add batch info
    train_adata = load_PBMC_data.load_PBMC_batch1_data(data_dir, batch = batch1)
    test_adata = load_PBMC_data.load_PBMC_batch1_data(data_dir, batch = batch2)
    train_adata.obs.index.name = None
    test_adata.obs.index.name = None

    ## curate given sub-cell types to major cell types
    train_adata = load_PBMC_data.curate_PBMC_demulx_celltypes(train_adata, celltype_gran)
    test_adata = load_PBMC_data.curate_PBMC_demulx_celltypes(test_adata, celltype_gran)
    return train_adata, test_adata


def process_batch1_batchtoind(data_dir, result_dir, batch="A", ind="1085", 
        sample=True, sample_seed=0, celltype_gran=1):
    '''Process batch1 A -> S5, as opposed to S1 -> S5

    @batch: batch A as training datasets
    @ind: ind S5 as predictor
    @sample: whether to downsample or not
    @sample_seed: random seed for sampling
    '''
    ## get training dataset
    train_adata = load_PBMC_data.load_PBMC_batch1_data(data_dir, batch=batch)
    train_adata.obs.index.name = None

    ## get test dataset
    test_adata = load_PBMC_data.load_PBMC_batch1_data(data_dir, ind=ind)
    test_adata.obs.index.name = None

    ## downsample training dataset to S1 = 1,551 cells
    if sample:
        random.seed(sample_seed)
        sampled_cells = random.sample(list(train_adata.obs_names), k=1551)
        train_adata = train_adata[sampled_cells]

    ## curate given sub-cell types to major cell types
    train_adata = load_PBMC_data.curate_PBMC_demulx_celltypes(train_adata, celltype_gran)
    test_adata = load_PBMC_data.curate_PBMC_demulx_celltypes(test_adata, celltype_gran)
    return train_adata, test_adata


# ---- process batch2 control and stimulated
def process_batch2_ctrl_stim(data_dir, result_dir, cond1="control", cond2="stimulated", 
       celltype_gran=1):
    '''Extract control as train adata, stimulated as test adata

    @cond1/cond2: control or stimultaed
    '''
    train_adata = load_PBMC_data.load_PBMC_batch2_data(data_dir, condition=cond1)
    test_adata = load_PBMC_data.load_PBMC_batch2_data(data_dir, condition=cond2)

    ## curate given sub-cell types to major cell types
    train_adata = load_PBMC_data.curate_PBMC_demulx_celltypes(train_adata, celltype_gran)
    test_adata = load_PBMC_data.curate_PBMC_demulx_celltypes(test_adata, celltype_gran)
    return train_adata, test_adata

# ---- process PBMC batch1/2 individuals
def process_batch1_batch2_ind(data_dir, result_dir, input1, input2, 
        celltype_gran=1):
    '''Use individuals from one batch to predict another batch

    @input1/input2: can be batch1_indID/batc2_indID
    @celltype_gran: 0 major; 1 sub
    '''
    ## split input1
    input1_list = input1.split('_')
    input1_batch = input1_list[0]
    if len(input1_list) > 1:
        input1_inds = '_'.join(input1_list[1:])
    else:
        input1_inds = None

    input2_list = input2.split('_')
    input2_batch = input2_list[0]
    if len(input2_list) > 1:
        input2_inds = '_'.join(input2_list[1:])
    else:
        input2_inds = None

    ## extract train and test adata according to batch information
    if input1_batch == "batch1":
        train_adata = load_PBMC_data.load_PBMC_batch1_data(data_dir, ind=input1_inds)
    elif "batch2" in input1_batch:
        cond = input1_batch.replace("batch2", "")
        if cond == "":
            train_adata = load_PBMC_data.load_PBMC_batch2_data(data_dir, ind=input1_inds)
        else:
            train_adata = load_PBMC_data.load_PBMC_batch2_data(data_dir, ind=input1_inds,
                condition=cond)

    if input2_batch == "batch1":
        test_adata = load_PBMC_data.load_PBMC_batch1_data(data_dir, ind=input2_inds)
    elif "batch2" in input2_batch:
        cond = input2_batch.replace("batch2", "")
        if cond == "":
            test_adata = load_PBMC_data.load_PBMC_batch2_data(data_dir, ind=input2_inds)
        else:
            test_adata = load_PBMC_data.load_PBMC_batch2_data(data_dir, ind=input2_inds,
                condition=cond)

    ## curate given sub-cell types to major cell types
    train_adata = load_PBMC_data.curate_PBMC_demulx_celltypes(train_adata, celltype_gran)
    test_adata = load_PBMC_data.curate_PBMC_demulx_celltypes(test_adata, celltype_gran)
    return train_adata, test_adata


# ---- process PBMC 7 protocols
def process_PBMC_protocols_type(data_dir, result_dir, 
        protocols1="Smart-seq2", protocols2="CEL-Seq2", exp="pbmc1"):
    '''Compare protocols using PBMC1/PBMC2/mixture

    @ protocols1/protocols2: 10x Chromium (v2)/10x Chromium (v2) A/10x Chromium (v2) B/
        10x Chromium (v3)/CEL-Seq2/Drop-seq/Seq-Well/Smart-seq2/inDrops
    @ exp: pbmc1/pbmc2/pbmc1_pbmc2
    '''
    ## load PBMC protocols adata
    adata = load_PBMCprotocol_data.load_PBMC_protocols_data(data_dir)
    exp_cells = adata.obs[adata.obs['Experiment'].isin(exp.split('_'))].index
    exp_adata = adata[exp_cells]

    train_cells = exp_adata.obs[exp_adata.obs['Method'].isin(protocols1.split('_'))].index
    train_adata = exp_adata[train_cells]
    test_cells = exp_adata.obs[exp_adata.obs['Method'].isin(protocols2.split('_'))].index
    test_adata = exp_adata[test_cells]
    return train_adata, test_adata

# ---- process PBMC 7 protocols batch
def process_PBMC_protocols_batch(data_dir, result_dir,
        batch1="pbmc1", batch2="pbmc2", protocol="Smart-seq2"):
    '''Compare same protocol between batch1 and batch2

    @ protocol: 10x Chromium (v2)/10x Chromium (v2) A/10x Chromium (v2) B/
        10x Chromium (v3)/CEL-Seq2/Drop-seq/Seq-Well/Smart-seq2/inDrops
    @ batch1/batch2: pbmc1/pbmc2
    '''
    ## load PBMC protocols adata
    adata = load_PBMCprotocol_data.load_PBMC_protocols_data(data_dir)
    pro_cells = adata.obs[adata.obs['Method'].isin(protocol.split('_'))].index
    pro_adata = adata[pro_cells]

    train_cells = pro_adata.obs[pro_adata.obs['Experiment'].isin([batch1])].index
    train_adata = pro_adata[train_cells]
    test_cells = pro_adata.obs[pro_adata.obs['Experiment'].isin([batch2])].index
    test_adata = pro_adata[test_cells]
    return train_adata, test_adata

# ---- process PBMC Zheng FACS-sorted data
def process_PBMC_Zheng(data_dir, result_dir, train_pct=0.8, sample_seed=0,
        curate=False):
    '''Randomly select train/test dataset to do cross-validation

    @train_pct: percentage of cells selected for training, the rest is remained for testing
    '''
    adata = load_PBMCZheng_data.load_PBMCZheng_data(data_dir, curate=curate) ## no need to curate data for cross-validation
    random.seed(sample_seed)
    train_cells = random.sample(adata.obs_names.tolist(), round(train_pct*adata.shape[0]))
    train_adata = adata[train_cells]
    test_cells = list(set(adata.obs_names) - set(train_cells))
    test_adata = adata[test_cells]
    return train_adata, test_adata

# --- cross dataset prediction using PBMC Zheng as target
def process_PBMC_comparison(data_dir, result_dir, 
        train="Kang", test="Zheng", sample_seed=0):
    '''Use curated PBMC Zheng dataset as target (B cells, CD14+ Monocytes, NK cells, CD4 T cells, CD8 T cells)

    @train: 
        - Kang: select one individual from healthy as reference, 10X
        - Ding: select pbmc fresh sample (pbmc2) sequenced by 10X as reference
    '''
    input = train.split('_')
    dataset = input[0]
    infos = None
    if len(input) > 1:
        infos = input[1:]

    ## build train dataset
    if "Kang" == dataset:
        if infos is None:
            train_adata = load_PBMC_data.load_PBMC_batch1_data(data_dir, ind="1154")
        elif infos[0] == "batch1":
            train_adata = load_PBMC_data.load_PBMC_batch1_data(data_dir)
    if "Ding" == dataset:
        if infos is None:
            train_adata = load_PBMCprotocol_data.load_PBMC_protocols_data(data_dir, 
                exp="pbmc2", protocol="10x-v2", curate=True)
        elif infos[0] == "droplet":
            train_adata = load_PBMCprotocol_data.load_PBMC_protocols_data(data_dir,
                exp="pbmc2", protocol_type="droplet", curate=True)

    ## build test dataset
    adata = load_PBMCZheng_data.load_PBMCZheng_data(data_dir, curate=True)
    random.seed(sample_seed)
    train_cells = random.sample(adata.obs_names.tolist(), round(0.8*adata.shape[0]))
    test_cells = list(set(adata.obs_names) - set(train_cells))
    test_adata = adata[test_cells]
    return train_adata, test_adata
