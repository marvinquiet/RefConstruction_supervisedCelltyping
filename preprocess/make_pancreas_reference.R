suppressMessages(library(Matrix))
suppressMessages(library(SingleCellExperiment))

## make reference for muraro/xin/seg
pancreas_ref <- function(data_dir, sample="seg", method="scmap") {
    ### ======= build reference for pancreas data
    # @data_dir: data directory i.e: "/home/wma36/gpu/data"
    # @sample: seg/xin/muraro
    # @method: scmap/CHETAH
    ### =======

    pancreas_dir <- file.path(data_dir, paste0("Pancreas6/", sample))
    if (sample == "seg") {
        counts <- read.csv(file.path(pancreas_dir, "Seg_counts.csv"), row.names=1)
        labels <- read.csv(file.path(pancreas_dir, "Seg_label.csv"), row.names=1)
    } else if (sample == "xin") {
        counts <- read.csv(file.path(pancreas_dir, "Xin_counts.csv"), row.names=1)
        labels <- read.csv(file.path(pancreas_dir, "Xin_label.csv"), row.names=1)
    } else if (sample == "muraro") {
        counts <- read.csv(file.path(pancreas_dir, "Muraro_counts.csv"), row.names=1)
        labels <- read.csv(file.path(pancreas_dir, "Muraro_label.csv"), row.names=1)
    }

    if ("scmap" == method) {
        colData <- DataFrame(cell_type1 = labels$x)
    }
    if ("CHETAH" == method) {
        colData <- DataFrame(celltypes = labels$x)
    }

    rownames(colData) <- colnames(counts)
    rownames(counts) <- toupper(rownames(counts)) ## convert to uppercase
    ref <- SingleCellExperiment(assays = list(counts = as.matrix(counts)),
                                colData = colData)
    return(ref)
}


pancreas_seg_ref <- function(data_dir, cond="Healthy", method="scmap") {
    ### ======= build reference for pancreas data
    # @data_dir: data directory i.e: "/home/wma36/gpu/data"
    # @cond: Healthy/T2D
    # @method: scmap/CHETAH
    ### =======

    ref <- pancreas_ref(data_dir, sample="seg", method=method)

    colData(ref)$condition <- ifelse(grepl("T2D", colnames(ref)), "T2D", "Healthy")
    idx <- which(colData(ref)$condition == cond)

    ref <- ref[, idx]
    return(ref)
}

pancreas_seg_mix_ref <- function(data_dir, main_cond="Healthy", 
                                 pred_cond="T2D", method="scmap") {
    ### ======= build reference for pancreas mixed data
    # @data_dir: data directory i.e: "/home/wma36/gpu/data"
    # @main_cond: one main condition mixed with another to predict another
    # @pred_cond: the predicted one
    # @method: scmap/CHETAH
    ### =======

    ref <- pancreas_ref(data_dir, sample="seg", method=method)

    colData(ref)$condition <- ifelse(grepl("T2D", colnames(ref)), "T2D", "Healthy")
    colData(ref)$ind <- sapply(strsplit(colnames(ref), "_"), '[', 1)

    idx <- which(colData(ref)$condition == pred_cond)
    pred_inds <- colData(ref)[idx, ]$ind
    pred_ind <- pred_inds[sort.list(pred_inds)][1]
    cat("===In seg mixed, predicting", pred_ind, "\n")

    main_idx <- which(colData(ref)$ind != pred_ind)
    pred_idx <- which(colData(ref)$ind == pred_ind)

    return(list(ref_data=ref[, main_idx], target_data=ref[, pred_idx]))
}
