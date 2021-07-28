source("pipelines/Rcode/run_scmap.R")
source("pipelines/Rcode/run_CHETAH.R")
source("pipelines/Rcode/utils.R")
#source("pipelines/Rcode/run_SC3.R")

data_dir <- "/home/wma36/gpu/data"
partition_features = NA  ## no partition on scmap/CHETAH anymore

## read in arguments
args <- commandArgs(trailingOnly = TRUE)
data_source <- args[1]
method <- args[2]
select_on <- args[3]
select_method <- args[4]
n_features <- args[5]
train <- args[6]
test <- args[7]

if (data_source %in% c("PBMC_batch1_ind", "PBMC_batch1_ABC", "PBMC_batch2", 
                "PBMC_batch1_batchtoind", "PBMC_protocols_pbmc1", "PBMC_protocols_batch_smart"))
    pipeline_dir = "pipelines/result_PBMC_collections"
if (data_source %in% c("PBMC_protocols_pbmc1", "PBMC_protocols_batch_smart"))
    pipeline_dir = "result_PBMC_protocols_collections"
if (data_source %in% c("PBMC_Zheng_FACS", "PBMC_Zheng_FACS_curated", "PBMC_cross"))
    pipeline_dir = "pipelines/result_PBMC_Zheng_collections"
if (data_source %in% c("pancreas", "pancreas_seg_cond", "pancreas_custom",
                "pancreas_seg_mix", "pancreas_multi_to_multi"))
    pipeline_dir = "pipelines/result_Pancreas_collections"
if (data_source %in% c("mousebrain_FC", "mousebrain_FC_sub", "mousebrain_HC", "mousebrain_HC_sub",
                "mousebrain_region", "mousebrain_FC_stage", "mousebrain_FC_stage_sub", "mousebrain_FC_datasets",
                "mousebrain_FC_datasets_multiinds", "mousebrain_FC_datasets_multiinds_sample",
                "mousebrain_FC_multiinds", "mousebrain_FC_multiinds_sub", "mousebrain_FC_multiinds_sample",
                "mousebrain_FC_multiinds_sub_sample"))
    pipeline_dir = "pipelines/result_Mousebrain_collections"
if (data_source %in% c("allenbrain_ss", "allenbrain_10x", "allenbrain_cross"))
    pipeline_dir = "pipelines/result_Allenbrain_collections"

result_prefix <- file.path(pipeline_dir, 
                paste("result", data_source, train, 'to', test, sep='_'))
if (select_on == "NA" && select_method == "NA") {
    result_dir <- file.path(result_prefix, "no_feature")
} else {
    result_dir <- file.path(result_prefix, 
                paste(select_method, n_features, 'on', select_on, sep='_'))
}
dir.create(result_dir)

## make reference
source("preprocess/make_PBMC_reference.R")
if ("PBMC_batch2" == data_source) {
    ref_data <- PBMC_batch2_ref(data_dir, condition=train, method=method)
    target_data <- PBMC_batch2_ref(data_dir, condition=test, method=method)
} 
if ("PBMC_batch1_ABC" == data_source) {
    ref_data <- PBMC_batch1_ref(data_dir, sample=train, method=method)
    target_data <- PBMC_batch1_ref(data_dir, sample=test, method=method)
} 
if ("PBMC_batch1_ind" == data_source) {
    ref_data <- PBMC_batch1_ind_ref(data_dir, ind=train, method=method)
    target_data <- PBMC_batch1_ind_ref(data_dir, ind=test, method=method)
}

if ("PBMC_batch1_batchtoind" == data_source) {
    ref_data <- PBMC_batch1_ref(data_dir, sample=train, method=method)
    ## downsample ref_data to S1=1,551 cells
    ref_data <- ref_downsample(ref_data, size=1551)
    target_data <- PBMC_batch1_ind_ref(data_dir, ind=test, method=method)
}

if ("PBMC_protocols_pbmc1" == data_source) {
    ref <- PBMC_protocols_ref(data_dir, method=method)

    train_idx <- colData(ref)$Method == train & colData(ref)$Experiment == "pbmc1"
    ref_data <- ref[, train_idx]

    test_idx <- colData(ref)$Method == test & colData(ref)$Experiment == "pbmc1"
    target_data <- ref[, test_idx]

    ## curate for same cell types
    res <- ref_target_common_cells(ref_data, target_data, method=method)
    ref_data <- res$ref
    target_data <- res$target
}
if ("PBMC_protocols_batch_smart" == data_source) {
    ref <- PBMC_protocols_ref(data_dir, method=method)

    train_idx <- colData(ref)$Method == "Smart-seq2" & colData(ref)$Experiment == train
    ref_data <- ref[, train_idx]

    test_idx <- colData(ref)$Method == "Smart-seq2" & colData(ref)$Experiment == test
    target_data <- ref[, test_idx]

    ## curate for same cell types
    res <- ref_target_common_cells(ref_data, target_data, method=method)
    ref_data <- res$ref
    target_data <- res$target
}

if ("PBMC_Zheng_FACS" == data_source || "PBMC_Zheng_FACS_curated" == data_source) {
    if (grepl("curated", data_source))
        ref = PBMC_Zheng_ref(data_dir, method=method, curate=TRUE)
    else
        ref = PBMC_Zheng_ref(data_dir, method=method)

    ref_data = ref_downsample(ref, size=round(ncol(ref)*as.numeric(train)))
    target_data = ref[, -which(colnames(ref) %in% colnames(ref_data))]
}

if ("PBMC_cross" == data_source) {
    ref = PBMC_Zheng_ref(data_dir, method=method, curate=TRUE)
    data = ref_downsample(ref, size=round(ncol(ref)*0.8))
    target_data = ref[, -which(colnames(ref) %in% colnames(data))]

    input = unlist(strsplit(train, '_'))
    dataset = input[1]
    infos = input[2]

    ## depends on when train is Kang and Ding
    if ("Kang" == dataset)
        if (is.na(infos)) {
            ref_data = PBMC_batch1_ind_ref(data_dir, ind="1154", method=method)
        } else if (infos == "batch1") {
            ref_data = PBMC_batch1_ind_ref(data_dir, method=method)
        }
    if ("Ding" == dataset)
        if (is.na(infos)) {
            ref_data = PBMC_protocols_ref(data_dir, method=method, 
                            exp="pbmc2", protocol="10x-v2", curate=T)
        } else if (infos == "droplet") {
            ref_data = PBMC_protocols_ref(data_dir, method=method,
                            exp="pbmc2", protocol_type=infos, curate=T)
        }

    ## curate for same cell types
    res <- ref_target_common_cells(ref_data, target_data, method=method)
    ref_data <- res$ref
    target_data <- res$target
}

source("preprocess/make_pancreas_reference.R")
if ("pancreas" == data_source) {
    ref_data <- pancreas_ref(data_dir, sample=train, method=method)
    target_data <- pancreas_ref(data_dir, sample=test, method=method)
}
if ("pancreas_seg_cond" == data_source) {
    ref_data <- pancreas_seg_ref(data_dir, cond=train, method=method)
    target_data <- pancreas_seg_ref(data_dir, cond=test, method=method)
}
if ("pancreas_custom" == data_source) {
    ref_data <- pancreas_seg_ref(data_dir, cond=train, method=method)
    target_data <- pancreas_ref(data_dir, sample=test, method=method)
}
if ("pancreas_seg_mix" == data_source) {
    res <- pancreas_seg_mix_ref(data_dir, main_cond=train, pred_cond=test, method=method)
    ref_data <- res$ref_data
    target_data <- res$target_data
}
if ("pancreas_multi_to_multi" == data_source) {
    cond1 <- unlist(strsplit(train, '_'))
    cond2 <- unlist(strsplit(test, '_'))

    ref_data <- NA; target_data <- NA
    for (i in seq(1, length(cond1))) {
        cond <- cond1[i]
        sce <-  pancreas_ref(data_dir, sample=cond, method=method)

        if (i == 1) {
            ref_data <- sce
        } else {
            common_genes <- intersect(rownames(ref_data), rownames(sce))
            ref_data <- ref_data[common_genes,]
            sce <- sce[common_genes,]

            ref_data <- cbind(ref_data, sce)
        }
    }

    for (i in seq(1, length(cond2))) {
        cond <- cond2[i]
        sce <-  pancreas_ref(data_dir, sample=cond, method=method)

        if (i == 1) {
            target_data <- sce
        } else {
            common_genes <- intersect(rownames(target_data), rownames(sce))
            target_data <- target_data[common_genes,]
            sce <- sce[common_genes,]

            target_data <- cbind(target_data, sce)
        }
    }
}

source("preprocess/make_mousebrain_reference.R")
if ("mousebrain_FC" == data_source) {
    ref_data <- mousebrain_ind_ref(data_dir, region="FC", ind=train, method=method)
    target_data <- mousebrain_ind_ref(data_dir, region="FC", ind=test, method=method)
}
if ("mousebrain_FC_sub" == data_source) {
    ## run using 81 sub-cell types
    ref_data <- mousebrain_ind_ref(data_dir, region="FC", ind=train, 
                                      method=method, celltype_gran=1)
    target_data <- mousebrain_ind_ref(data_dir, region="FC", ind=test, 
                                         method=method, celltype_gran=1)
}
if ("mousebrain_FC_multiinds" == data_source) {
    ## exclude input train
    FC_inds <- c("P60FCAldh1l1Rep1", "P60FCCx3cr1Rep1", "P60FCRep1", "P60FCRep2",
                 "P60FCRep3", "P60FCRep4", "P60FCRep6")
    exclude_list <- FC_inds[FC_inds!=test]
    exclude_inds <- paste(exclude_list, collapse="_")

    ref_data <- mousebrain_ind_ref(data_dir, region="FC", ind=exclude_inds, method=method)
    target_data <- mousebrain_ind_ref(data_dir, region="FC", ind=test, method=method)
}

if ("mousebrain_FC_multiinds_sub" == data_source) {
    ## exclude input train
    FC_inds <- c("P60FCAldh1l1Rep1", "P60FCCx3cr1Rep1", "P60FCRep1", "P60FCRep2",
                 "P60FCRep3", "P60FCRep4", "P60FCRep6")
    exclude_list <- FC_inds[FC_inds!=test]
    exclude_inds <- paste(exclude_list, collapse="_")

    ref_data <- mousebrain_ind_ref(data_dir, region="FC", ind=exclude_inds, 
                                   method=method, celltype_gran=1)
    target_data <- mousebrain_ind_ref(data_dir, region="FC", ind=test, 
                                   method=method, celltype_gran=1)
}

if ("mousebrain_FC_multiinds_sample" == data_source) {
    ## exclude input train
    FC_inds <- c("P60FCAldh1l1Rep1", "P60FCCx3cr1Rep1", "P60FCRep1", "P60FCRep2",
                 "P60FCRep3", "P60FCRep4", "P60FCRep6")
    exclude_list <- FC_inds[FC_inds!=test]
    exclude_inds <- paste(exclude_list, collapse="_")

    ref_data <- mousebrain_ind_ref(data_dir, region="FC", ind=exclude_inds, method=method)
    tmp_data <- mousebrain_ind_ref(data_dir, region="FC", ind=train, method=method)

    ## downsample to 8062 cells which is sample ind
    ref_data <- ref_downsample(ref_data, size=dim(tmp_data)[1])
    target_data <- mousebrain_ind_ref(data_dir, region="FC", ind=test, method=method)
}

if ("mousebrain_FC_multiinds_sub_sample" == data_source) {
    ## exclude input train
    FC_inds <- c("P60FCAldh1l1Rep1", "P60FCCx3cr1Rep1", "P60FCRep1", "P60FCRep2",
                 "P60FCRep3", "P60FCRep4", "P60FCRep6")
    exclude_list <- FC_inds[FC_inds!=train]
    exclude_inds <- paste(exclude_list, collapse="_")

    ref_data <- mousebrain_ind_ref(data_dir, region="FC", ind=exclude_inds, 
                                   method=method, celltype_gran=1)

    tmp_data <- mousebrain_ind_ref(data_dir, region="FC", ind=train, method=method)
    ## downsample to 8062 cells which is sample ind
    ref_data <- ref_downsample(ref_data, size=dim(tmp_data)[1])
    target_data <- mousebrain_ind_ref(data_dir, region="FC", ind=test, 
                                   method=method, celltype_gran=1)
}
if ("mousebrain_HC" == data_source) {
    ref_data <- mousebrain_ind_ref(data_dir, region="HC", ind=train, method=method)
    target_data <- mousebrain_ind_ref(data_dir, region="HC", ind=test, method=method)
}
if ("mousebrain_HC_sub" == data_source) {
    ref_data <- mousebrain_ind_ref(data_dir, region="HC", ind=train, 
                                      method=method, celltype_gran=1)
    target_data <- mousebrain_ind_ref(data_dir, region="HC", ind=test, 
                                         method=method, celltype_gran=1)
}
if ("mousebrain_region" == data_source || "mousebrain_region_sub" == data_source) {
    celltype_gran <- 0
    if ("mousebrain_region_sub" == data_source) celltype_gran <- 1

    celltype_ind <- ifelse("scmap" == method, "cell_type1", "celltypes")

    if ("FC" == train || "FC" == test) {
        FC_data <- mousebrain_ind_ref(data_dir, region="FC", ind=NA, method=method,
                                      celltype_gran=celltype_gran)
        FC_colData <- colData(FC_data)
        idx <- which(FC_colData[, celltype_ind] == "Interneuron_CGE" | FC_colData[, celltype_ind] == "Interneuron_MGE")
        FC_colData[idx, celltype_ind] <- "Interneuron"

        idx <- which(FC_colData[, celltype_ind] == "Neuron_Claustrum" | 
                     FC_colData[, celltype_ind] == "Neuron_L2/3" |
                     FC_colData[, celltype_ind] == "Neuron_L5" | 
                     FC_colData[, celltype_ind] == "Neuron_L5b" |
                     FC_colData[, celltype_ind] == "Neuron_L6")
        FC_colData[idx, celltype_ind] <- "Neuron"
        FC_colData[, celltype_ind] <- as.character(FC_colData[, celltype_ind])
        FC_colData[, celltype_ind] <- as.factor(FC_colData[, celltype_ind])
        colData(FC_data) <- FC_colData

        if ("FC" == train) ref_data <- FC_data
        else if ("FC" == test) target_data <- FC_data
    }

    if ("HC" == train || "HC" == test) {
        HC_data <- mousebrain_ind_ref(data_dir, region="HC", ind=NA, method=method,
                                      celltype_gran=celltype_gran)
        HC_colData <- colData(HC_data)
        ## remove certain rows
        HC_colData <- subset(HC_colData, HC_colData[, celltype_ind] != "Choroid_Plexus" & 
                     HC_colData[, celltype_ind] != "Ependyma" &
                     HC_colData[, celltype_ind] != "Neurongenesis_Mitosis")
        HC_colData[, celltype_ind] <- as.character(HC_colData[, celltype_ind])
        HC_colData[, celltype_ind] <- as.factor(HC_colData[, celltype_ind])
        HC_data <- HC_data[, rownames(HC_colData)]
        colData(HC_data) <- HC_colData

        if ("HC" == train) ref_data <- HC_data
        else if ("HC" == test) target_data <- HC_data
    }
}
if ("mousebrain_FC_stage" == data_source) {
    ref_data <- mousebrain_FC_stage_ref(data_dir, stage=train, method=method)
    target_data <- mousebrain_FC_stage_ref(data_dir, stage=test, method=method)
}
if ("mousebrain_FC_stage_sub" == data_source) {
    ref_data <- mousebrain_FC_stage_ref(data_dir, stage=train, method=method, celltype_gran=1)
    target_data <- mousebrain_FC_stage_ref(data_dir, stage=test, method=method, celltype_gran=1)

    ## curate for same 38 sub-cell types
    res <- ref_target_common_cells(ref_data, target_data, method=method)
    ref_data <- res$ref
    target_data <- res$target
}
if ("mousebrain_FC_datasets" == data_source) {
    FC_data <- mousebrain_ind_ref(data_dir, region="FC", ind=NA, method=method)
    target_data <- mousebrain_FC_stage_ref(data_dir, stage=test, method=method)

    res <- curate_common_FC_celltypes(FC_data, target_data, method=method)
    ref_data <- res$dataset1
    target_data <- res$dataset2

    ## curate for same 7 major cell types
    res <- ref_target_common_cells(ref_data, target_data, method=method)
    ref_data <- res$ref
    target_data <- res$target
}

if ("mousebrain_FC_datasets_multiinds" == data_source | 
    "mousebrain_FC_datasets_multiinds_sample" == data_source) {
    dataset1_input <- unlist(strsplit(train, '_'))
    dataset2_input <- unlist(strsplit(test, '_'))

    dataset1_region <- dataset1_input[1]
    if (length(dataset1_input) > 1) {
        dataset1_inds <- dataset1_input[2:length(dataset1_input)]
        dataset1_inds <- paste(dataset1_inds, sep='_')
    } else {
        dataset1_inds <- NA
    }

    dataset2_stage <- dataset2_input[1]
    if (length(dataset1_input) > 1) {
        dataset2_inds <- dataset2_input[2:length(dataset2_input)]
        dataset2_inds <- paste(dataset2_inds, sep='_')
    } else {
        dataset2_inds <- NA
    }

    ref_data <- mousebrain_ind_ref(data_dir, region=dataset1_region, 
                                   ind=dataset1_inds, method=method)
    target_data <- mousebrain_FC_stage_ref(data_dir, stage=dataset2_stage, 
                                    ind=dataset2_inds, method=method)

    res <- curate_common_FC_celltypes(ref_data, target_data, method=method)
    ref_data <- res$dataset1
    target_data <- res$dataset2

    ## curate for same major cell types
    res <- ref_target_common_cells(ref_data, target_data, method=method)
    ref_data <- res$ref
    target_data <- res$target

    ## downsample to certain amount
    if (grepl("sample", data_source)) {
        ref_data <- ref_downsample(ref_data, size=7972)
    }
}

if ("allenbrain_ss" == data_source) {
    ## Allen brain Smart-Seqv4 80% as training, 20% as test
    ref = allenbrain_ss_ref(data_dir, method=method, curate=TRUE)
    ref_data = ref_downsample(ref, size=round(ncol(ref)*as.numeric(train)))
    target_data = ref[, -which(colnames(ref) %in% colnames(ref_data))]
}

if ("allenbrain_10x" == data_source) {
    ## Allen brain 10X one individual to predict another
    ref_data = allenbrain_10x_ref(data_dir, method=method, ind=train, curate=TRUE) 
    target_data = allenbrain_10x_ref(data_dir, method=method, ind=test, curate=TRUE)
}

if ("allenbrain_cross" == data_source) {
    ## Allen brain Smart-Seqv4 80% as training, 20% as test
    ref = allenbrain_ss_ref(data_dir, method=method, curate=TRUE)
    ref_data = ref_downsample(ref, size=round(ncol(ref)*0.8))
    target_data = ref[, -which(colnames(ref) %in% colnames(ref_data))]

    ## FC/pFC as reference
    if ("FC" == train) {
        ref_data = mousebrain_ind_ref(data_dir, region="FC", ind="P60FCRep1", 
                        method=method, curate=T)
    }
    if ("pFC" == train) {
        ref_data = mousebrain_FC_stage_ref(data_dir, stage="Adult", 
                        ind="PFCSample1", method=method, curate=T)
    }
    res <- ref_target_common_cells(ref_data, target_data, method=method)
    ref_data <- res$ref
    target_data <- res$target
}


## extract number of clusters from reference data
if ("scmap" == method) {
    n_clusters <- length(unique(colData(ref_data)$cell_type1))
} else if ("CHETAH" == method) {
    n_clusters <- length(unique(colData(ref_data)$celltypes))
}
cat("Number of cell types: ", n_clusters, '\n')

## filter data to express in at least 10 cells and 10 genes
ref_data = filter_sce(ref_data)
target_data = filter_sce(target_data)
common_genes = intersect(rownames(ref_data), rownames(target_data))
ref_data = ref_data[common_genes,]
target_data = target_data[common_genes,]

ptm <- proc.time()
if (method == "scmap") {
    cat("\n\n=== scmap:\n")
    res <- run_scmap_cluster(ref_data, target_data, result_dir, partition_features,
                         prefix=method, n_features=n_features, n_clusters=n_clusters, 
                         select_on=select_on, select_method=select_method)
    run_time <- proc.time() - ptm
    cat("\n\n=== Run time:", run_time, "\n")
    analyze_scmap_result(target_data, res, result_dir, 
                         prefix=method, run_time=unname(run_time[1]))
} else if (method == "CHETAH") {
    cat("\n\n=== CHETAH:\n")
    res <- run_CHETAH(ref_data, target_data, result_dir, partition_features,
                         prefix=method, n_features=n_features, n_clusters=n_clusters, 
                         select_on=select_on, select_method=select_method)
    run_time <- proc.time() - ptm
    cat("\n\n=== Run time:", run_time, "\n")
    analyze_CHETAH_result(target_data, res, result_dir, 
                          prefix=method, run_time=unname(run_time[1]))
}


if (FALSE) {
## run SC3 on test datasets and label the individuals
run_SC3_analysis(target_data, k=n_clusters)

common_genes <- intersect(rownames(ref_data), rownames(target_data))
colData(ref_data)$ind <- as.factor(sapply(strsplit(colnames(ref_data), '_'), '[', 1))
colData(target_data)$condition <- "Healthy"
colData(target_data)$ind <- as.factor(sapply(strsplit(colnames(target_data), '\\.'), '[', 1))
sce <- cbind(ref_data[common_genes, ], target_data[common_genes, ])

run_SC3_analysis(common_genes)
}

