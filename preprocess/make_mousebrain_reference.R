suppressMessages(library(Matrix))
suppressMessages(library(SingleCellExperiment))

## make reference for mouse brain
mousebrain_ind_ref <- function(data_dir, region="FC", ind="P60FCRep1", method="scmap",
                                  celltype_gran=0, curate=F) {
    ### === build reference for mousebrain FC data
    # @data_dir: data directory i.e: "/home/wma36/gpu/data"
    # @region: brain regions, e.g. FC/HC
    # @ind: extract individuals; if NA, return all, or using _ as separator to input many
    # @method: scmap/CHETAH
    # @celltype_gran: granularity of cell types. 0 represents major and 1 represents sub-celltypes
    ### === 

    mousebrain_dir <- file.path(data_dir, "/Mousebrain")

    ## load metadata
    mousebrain_metadata <- read.table(file.path(mousebrain_dir, "Mousebrain_metadata.csv"),
                                      sep=",", header=T, row.names=1)
    if (is.na(ind)) {
        mousebrain_metadata <- mousebrain_metadata[mousebrain_metadata$mouse_celltypes != "unknown", ]
    } else {
        inds <- unlist(strsplit(ind, '_'))
        idx <- mousebrain_metadata$mouse_celltypes != "unknown" & mousebrain_metadata$ind %in% inds
        mousebrain_metadata <- mousebrain_metadata[idx, ]
    }

    if ("FC" == region) {
        ## load FC data
        data <- readRDS(file.path(mousebrain_dir, "Cortex_noRep5_FRONTALonly.RDS"))
    } else if ("HC" == region) {
        ## load HC data
        data <- readRDS(file.path(mousebrain_dir, "Hippocampus.RDS"))
    }
    rownames(data) <- toupper(rownames(data))  ## convert genes to upper

    common_cells <- intersect(mousebrain_metadata$barcodes, colnames(data)) 
    data <- data[, common_cells]
    metadata <- mousebrain_metadata[match(common_cells, mousebrain_metadata$barcodes), ]
    
    if (curate) {
        metadata$mouse_celltypes = as.character(metadata$mouse_celltypes)
        idx = which(metadata[, "mouse_celltypes"] == "Interneuron_CGE" | 
                     metadata[, "mouse_celltypes"] == "Interneuron_MGE")
        metadata[idx, "mouse_celltypes"] = "Interneuron"
        idx = which(metadata[, "mouse_celltypes"] %in%
                c("Neuron_Claustrum", "Neuron_L2/3", "Neuron_L5", "Neuron_L5b", "Neuron_L6"))
        metadata[idx, "mouse_celltypes"] <- "Neuron"
        metadata$mouse_celltypes = as.factor(metadata$mouse_celltypes)
    }

    celltypes <- metadata$mouse_celltypes
    if (1 == celltype_gran) {
        subclusters <- metadata$cluster
        celltypes <- paste0(subclusters, '.', celltypes)
    }

    if ("scmap" == method) {
        colData <- DataFrame(cell_type1 = celltypes)
    }
    if ("CHETAH" == method) {
        colData <- DataFrame(celltypes = celltypes)
    }

    rownames(colData) <- colnames(data)
    ref <- SingleCellExperiment(assays = list(counts = as.matrix(data)),
                                colData = colData)
    return(ref)
}

## make reference for mouse FC
mousebrain_FC_stage_ref <- function(data_dir, stage=NA, ind=NA, method="scmap",
                                    celltype_gran=0, curate=F) {
    ### === build reference for mousebrain FC data
    # @data_dir: data directory i.e: "/home/wma36/gpu/data"
    # @stage: development stage, Adult or P21
    # @ind: extract individuals, PFCSample1-12, P21Sample1-3
    # @method: scmap/CHETAH
    # @celltype_gran: granularity of cell types. 0 represents major and 1 represents sub-celltypes
    ### === 
 
    mousebrain_FC_dir <- file.path(data_dir, "MouseFC_GSE124952")

    ## load metadata
    metadata <- read.table(file.path(mousebrain_FC_dir, "GSE124952_meta_data.csv"),
                                      sep=",", header=T, row.names=1)
    if (!is.na(stage)) {
        metadata <- metadata[metadata$DevStage == stage,]
    }

    if (!is.na(ind)) {
        inds <- unlist(strsplit(ind, '_'))
        metadata <- metadata[metadata$Sample %in% inds,]
    }
 
    data <- read.table(file.path(mousebrain_FC_dir, "GSE124952_expression_matrix.csv"),
                       sep=",", header=T, row.names=1)
    rownames(data) <- toupper(rownames(data))  ## convert genes to upper

    common_cells <- intersect(colnames(data), rownames(metadata))
    data <- data[, common_cells]
    metadata <- metadata[common_cells, ]

    if (0 == celltype_gran) {
        if (curate) {
            metadata$CellType = as.character(metadata$CellType)
            idx <- which(metadata[, "CellType"] == 'Oligo'|
                         metadata[, "CellType"] == 'NF Oligo')
            metadata[idx, "CellType"] <- "Oligodendrocytes"
            idx <- which(metadata[, "CellType"] == 'OPC')
            metadata[idx, "CellType"] <- "Polydendrocytes"
            idx <- which(metadata[, "CellType"] == 'Astro')
            metadata[idx, "CellType"] <- "Astrocytes"
            idx <- which(metadata[, "CellType"] == 'Excitatory')
            metadata[idx, "CellType"] <- "Neuron"
            idx <- which(metadata[, "CellType"] == 'Inhibitory')
            metadata[idx, "CellType"] <- "Interneuron"
            idx <- which(metadata[, "CellType"] == 'Endo')
            metadata[idx, "CellType"] <- "Endothelial"
            metadata$CellType = as.factor(metadata$CellType)
        }
        celltypes <- metadata$CellType
    } else if (1 == celltype_gran) {
        celltypes <- metadata$L2_clusters
    }

    if ("scmap" == method) {
        colData <- DataFrame(cell_type1 = celltypes)
    }
    if ("CHETAH" == method) {
        colData <- DataFrame(celltypes = celltypes)
    }

    rownames(colData) <- colnames(data)
    ref <- SingleCellExperiment(assays = list(counts = as.matrix(data)),
                                colData = colData)
    return(ref)
}

## find common cell types between datasets
curate_common_FC_celltypes <- function(mousebrain_dataset1_data, mousebrain_dataset2_data,
                                       method="scmap") {
    ### === curate common FC cell types for two datasets
    # @mousebrain_dataset1_data: sce object
    # @mousebrain_dataset2_data: sce object
    # return: two curated sce objects
    ### ===
    celltype_ind <- ifelse("scmap" == method, "cell_type1", "celltypes")
    dataset1_colData <- colData(mousebrain_dataset1_data)
    dataset1_colData[, celltype_ind] <- as.character(dataset1_colData[, celltype_ind])
    idx <- which(dataset1_colData[, celltype_ind] == "Interneuron_CGE" | 
                 dataset1_colData[, celltype_ind] == "Interneuron_MGE")
    dataset1_colData[idx, celltype_ind] <- "Interneuron"

    idx <- which(dataset1_colData[, celltype_ind] %in%
            c("Neuron_Claustrum", "Neuron_L2/3", "Neuron_L5", "Neuron_L5b", "Neuron_L6"))
    dataset1_colData[idx, celltype_ind] <- "Neuron"
    dataset1_colData[, celltype_ind] <- as.factor(dataset1_colData[, celltype_ind])
    colData(mousebrain_dataset1_data) <- dataset1_colData
    ref_data <- mousebrain_dataset1_data

    dataset2_colData <- colData(mousebrain_dataset2_data)
    dataset2_colData[, celltype_ind] <- as.character(dataset2_colData[, celltype_ind])

    idx <- which(dataset2_colData[, celltype_ind] == 'Oligo'|
                 dataset2_colData[, celltype_ind] == 'NF Oligo')
    dataset2_colData[idx, celltype_ind] <- "Oligodendrocytes"
    idx <- which(dataset2_colData[, celltype_ind] == 'OPC')
    dataset2_colData[idx, celltype_ind] <- "Polydendrocytes"
    idx <- which(dataset2_colData[, celltype_ind] == 'Astro')
    dataset2_colData[idx, celltype_ind] <- "Astrocytes"
    idx <- which(dataset2_colData[, celltype_ind] == 'Excitatory')
    dataset2_colData[idx, celltype_ind] <- "Neuron"
    idx <- which(dataset2_colData[, celltype_ind] == 'Inhibitory')
    dataset2_colData[idx, celltype_ind] <- "Interneuron"
    idx <- which(dataset2_colData[, celltype_ind] == 'Endo')
    dataset2_colData[idx, celltype_ind] <- "Endothelial"

    dataset2_colData[, celltype_ind] <- as.factor(dataset2_colData[, celltype_ind])
    mousebrain_dataset2_data <- mousebrain_dataset2_data[, rownames(dataset2_colData)]
    colData(mousebrain_dataset2_data) <- dataset2_colData

    return(list(dataset1=mousebrain_dataset1_data, 
                dataset2=mousebrain_dataset2_data))
}

## add reference for Allen brain Smart-Seqv4
allenbrain_ss_ref = function(data_dir, method="scmap", ind=NA, curate=FALSE) {
    ### === build reference panel for Allen Brain Smart-SeqV4 data
    # @ data_dir: data directory
    # @ method: scmap or CHETAH
    # @ curate: whether to curate the cell types
    ### === 
    allen_brain_dir = file.path(data_dir, "AllenBrain")
    x = read.table(file.path(allen_brain_dir, "filtered_SS_matrix.csv"), 
                      header=T, sep=',', row.names=1, check.names=F)
    data = t(x)
    rownames(data) <- toupper(rownames(data))  ## convert genes to upper 
    rm(x) ## release memory
    metadata = read.csv(file.path(allen_brain_dir, "filtered_SS_metadata.csv"))
    rownames(metadata) = metadata$barcode

    if (!is.na(ind)) {
        inds <- unlist(strsplit(ind, '_'))
        metadata <- metadata[metadata$donor_label %in% inds,]
    }

    if (curate) {
        metadata$label = as.character(metadata$label)
        idx = which(metadata$label %in% c('L2/3 IT CTX', 'L2/3 IT RSP', 'L2/3 IT RSPv',
                'L2 IT ENTl', 'L4/5 IT CTX', 'L4 RSP-ACA', 'L5/6 IT CTX', 'L5 IT CTX',
                'L5 IT TPE-ENT', 'L5 NP CT CTX ', 'L5 NP CTX', 'L5 PT CTX', 'L6b CTX',
                'L6 CT CTX', 'L6 IT CTX', 'CA1', 'CA1-ProS', 'CA2-IG-FC', 'CA3',
                'CT SUB', 'Car3', 'DG', 'L2 IT ENTm', 'L2 IT PAR', 'L2/3 IT APr',
                'L2/3 IT ENTl', 'L2/3 IT HATA', 'L2/3 IT PPP', 'L3 IT ENTl', 'L3 IT ENTm',
                'L5 PPP', 'L6 IT ENTl', 'L6b/CT ENT', 'NP PPP', 'NP SUB', 'SUB', 'SUB-ProS'))
        metadata[idx, ]$label = "Neuron"
        idx = which(metadata$label %in% c('Lamp5', 'Lamp5 Lhx6', 'Pax6', 'Pvalb',
                'Sncg', 'Sst', 'Sst Chodl', 'Vip'))
        metadata[idx, ]$label = "Interneuron"
        idx = which(metadata$label == "Astro")
        metadata[idx, ]$label = "Astrocytes"
        idx = which(metadata$label == "Endo")
        metadata[idx, ]$label = "Endothelial"
        idx = which(metadata$label == "Micro-PVM")
        metadata[idx, ]$label = "Microglia"
        idx = which(metadata$label == "Oligo")
        metadata[idx, ]$label = "Oligodendrocytes"
        idx = which(metadata$label == "SMC-Peri")
        metadata[idx, ]$label = "Pericytes"

        ## remove certain datasets
        remained_idx = which(!metadata$label %in% c('CR', 'Meis2', 'Meis2 HPF', 'Ndnf HPF', 'VLMC'))
        metadata = metadata[remained_idx, ]
        data = data[, rownames(metadata)]
        metadata$label = as.factor(metadata$label)
    }
 
    if (method == "CHETAH") 
        colData = DataFrame(celltypes=metadata$label)
    if (method == "scmap")
        colData = DataFrame(cell_type1=metadata$label)
    rownames(colData) = rownames(metadata)
    ref = SingleCellExperiment(assays = list(counts = as.matrix(data)),
                               colData = colData)
    return(ref)
}

## add reference for Allen brain 10X Genomics
allenbrain_10x_ref = function(data_dir, method="scmap", ind=NA, curate=FALSE) {
    ### === build reference panel for Allen Brain 10X data
    # @ data_dir: data directory
    # @ method: scmap or CHETAH
    # @ curate: whether to curate the cell types
    ### === 
    allen_brain_dir = file.path(data_dir, "AllenBrain")
    x = read.table(file.path(allen_brain_dir, "filtered_matrix.csv"), 
                      header=T, sep=',', row.names=1, check.names=F)
    data = t(x)
    rownames(data) <- toupper(rownames(data))  ## convert genes to upper
    rm(x) ## release memory
    metadata = read.table(file.path(allen_brain_dir, "filtered_metadata.csv"),
                        header=T, sep=',', row.names=1)
    rownames(metadata) = metadata$sample_name
    if (!is.na(ind)) {
        inds <- unlist(strsplit(ind, '_'))
        metadata <- metadata[metadata$external_donor_name_label %in% inds,]
    }

    if (curate) {
        metadata$cell.type = as.character(metadata$cell.type)
        idx = which(metadata$cell.type %in% c('L2 IT HATA', 'L2 IT RSP-ACA', 'L2/3 IT APr',
                'L2/3 IT CTX', 'L2/3 IT PPP', 'L3 RSP-ACA', 'L4/5 IT CTX', 'L5 IT CTX',
                'L5 NP CTX', 'L5 PT CTX', 'L5 PT RSP-ACA', 'L5/6 IT CTX', 'L6 CT CTX',
                'L6 Car3', 'L6 IT CTX', 'L6 NP CT CTX', 'L6b CTX','NP PPP', 'NP SUB',
                'CA2-IG-FC', 'L6b RHP', 'L6 IT RHP', 'L2/3 IT ENTl', 'L2 IT RSPv',
                'L5 IT TPE-ENT', 'L2/3 IT TPE', 'L2 IT ProS'))
        metadata[idx, "cell.type"] = "Neuron"
        idx = which(metadata$cell.type %in% c('Sncg', 'Pvalb', 'Pvalb Vipr2', 'Sst',
                'Sst Chodl', 'Vip', 'Lamp5', 'Lamp5 Lhx6', 'Pax6', 'Ndnf HPF'))
        metadata[idx, ]$cell.type = "Interneuron"
        idx = which(metadata$cell.type == "Astro")
        metadata[idx, ]$cell.type = "Astrocytes"
        idx = which(metadata$cell.type == "Endo")
        metadata[idx, ]$cell.type = "Endothelial"
        idx = which(metadata$cell.type == "Micro")
        metadata[idx, ]$cell.type = "Microglia"
        idx = which(metadata$cell.type == "OPC")
        metadata[idx, ]$cell.type = "Polydendrocytes"
        idx = which(metadata$cell.type == "Oligo")
        metadata[idx, ]$cell.type = "Oligodendrocytes"
        idx = which(metadata$cell.type == "Peri")
        metadata[idx, ]$cell.type = "Pericytes"

        ## remove certain datasets
        remained_idx = which(!metadata$cell.type %in% c('CR', 'PVM', 'SMC', 'VLMC'))
        metadata = metadata[remained_idx, ]
        data = data[, rownames(metadata)]
        metadata$cell.type = as.factor(metadata$cell.type)
    }

    if (method == "CHETAH") 
        colData = DataFrame(celltypes=metadata$cell.type)
    if (method == "scmap")
        colData = DataFrame(cell_type1=metadata$cell.type)
    ref = SingleCellExperiment(assays = list(counts = as.matrix(data)),
                               colData = colData)
    return(ref)
}
