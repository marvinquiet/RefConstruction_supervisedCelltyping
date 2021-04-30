suppressMessages(library(Matrix))
suppressMessages(library(SingleCellExperiment))

## make reference for mouse brain
mousebrain_ind_ref <- function(data_dir, region="FC", ind="P60FCRep1", method="scmap",
                                  celltype_gran=0) {
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
                                    celltype_gran=0) {
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

    idx <- which(dataset1_colData[, celltype_ind] == "Neuron_Claustrum" | 
                 dataset1_colData[, celltype_ind] == "Neuron_L2/3" |
                 dataset1_colData[, celltype_ind] == "Neuron_L5" | 
                 dataset1_colData[, celltype_ind] == "Neuron_L5b" |
                 dataset1_colData[, celltype_ind] == "Neuron_L6")
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
