suppressMessages(library(scmap))
suppressMessages(library(aricode))

scmap_normdata <- function(sce) {
    ## normalized data
    counts <- assay(sce, "counts")
    lib.sizes <- colSums(counts)
    scale.factor <- 1e4
    # size.factors <- lib.sizes/mean(lib.sizes)
    # normcounts <- t(t(counts)/size.factors)
    # normcounts(sce) <- normcounts
    normcounts(sce) <- scale(counts, center=FALSE, scale=lib.sizes) * scale.factor
    logcounts(sce) <- log2(normcounts(sce) + 1) ## log-normalize the data
    return(sce)
}

scmap_feature_selection <- function(sce, result_dir, partition_features=NA, 
            n_features=1000, n_clusters=NA, select_method="Seurat") {
    ### ---- Feature selection
    # @sce: sce object from scmap
    # @result_dir: load features already selected, or store selected features
    # @partition_features: the features provided by partition in GEDFN
    # @n_features: how many features will be selected
    # @select_method: scmap, Seurat, FEAST_SC3, FEAST_F-test

    # Return:
    # @sce: add scmap_features to the rowData(sce)
    ### ----

    features_file <- file.path(result_dir, "features.txt")
    features <- NA
    ## if features have already been selected
    if (file.exists(features_file)) {
        features <- scan(features_file, what=character())
    } else {
        if (select_method == "Seurat") {
            suppressMessages(library(Seurat))
            seuratObj <- as.Seurat(sce)
            Idents(seuratObj) <- "cell_type1"
            seuratObj <- NormalizeData(seuratObj, 
                normalization.method = "LogNormalize", scale.factor = 10000)
            seuratObj <- FindVariableFeatures(seuratObj, 
                selection.method = "mvp", nfeatures = n_features) # seurat v2
            features <- head(VariableFeatures(seuratObj), n_features)
            rm(seuratObj)  # release memory
        } else if (grepl("FEAST|F-test", select_method)) {
            suppressMessages(library(FEAST))
            counts <- assay(sce, "counts")
            processed_counts <- process_Y(counts, thre=0)
            rm(counts) # release memory

            if ("FEAST" == select_method)
                F_res <- Consensus_F2(processed_counts, top=n_features, k=n_clusters, 
                                split=T, batch_size=1000, num_cores=4)
            else if ("F-test" == select_method)
                F_res <- cal_F2(processed_counts, 
                                as.character(col(sce)$cell_type1))$F_scores

            names(F_res) <- rownames(processed_counts)
            rm(processed_counts) # release memory
            ixs <- order(F_res, decreasing=T)[1:n_features]
            features <- names(F_res[ixs])
        } 
        write(features, features_file)
    }

    ## for intersecting genes
    if (!any(is.na(partition_features))) {
        features <- intersect(partition_features, features)
    }
    cat("Number of features:", length(features), '\n')

    if (select_method == "scmap") {
        ## run scmap selecting features
        png(file.path(result_dir, paste(prefix, "feature_plot.png", sep="_")))
        sce <- selectFeatures(sce, n_features=n_features, suppress_plot = FALSE)
        dev.off()
    }

    if (!any(is.na(features))) {
        rowData(sce)$scmap_features <- rownames(sce) %in% features  ## set features as true
    }
    #cat("Features by scmap:", table(rowData(sce)$scmap_features), '\n')
    return(sce)
}

run_scmap_cluster <- function(ref_data, target_data, result_dir, partition_features=NA, 
                              prefix="scmap_cluster", n_features=1000, n_clusters=8,
                              select_on="train", select_method="Seurat", seed=0) {
    ### ---- Run scmap using ref_data as reference and predict cell types on target
    # @ref_data: reference SCE object
    # @target_data: target SCE object
    # @parrition: the partition from GEDFN
    # @result_dir: directory to store the result
    # @prefix: storing result

    ### ----
    set.seed(seed)

    # set feature symbols
    rowData(ref_data)$feature_symbol <- rownames(ref_data)
    rowData(target_data)$feature_symbol <- rownames(target_data)
    ref_data <- ref_data[!duplicated(rownames(ref_data)), ]  ## deduplicate genes
    target_data <- target_data[!duplicated(rownames(target_data)), ] 

    ## normalize data
    ref_data <- scmap_normdata(ref_data)
    target_data <- scmap_normdata(target_data)

    if (select_on == "train" || select_on == "NA") {
        ref_data <- scmap_feature_selection(ref_data, result_dir, partition_features,
            n_features, n_clusters, select_method)
    } else if (select_on == "test") {
        target_data <- scmap_feature_selection(target_data, result_dir, partition_features,
            n_features, n_clusters, select_method)
        target_features_idx <- which(rowData(target_data)$scmap_features == TRUE)
        target_features <- rownames(target_data)[target_features_idx]
        rowData(ref_data)$scmap_features <- rownames(ref_data) %in% target_features
    }
    cat("Features by scmap:", table(rowData(ref_data)$scmap_features), '\n')

    ## run scmap-cluster
    ref_data <- scmap::indexCluster(ref_data)
    png(file.path(result_dir, paste(prefix, "idxcluster_heatmap.png", sep="_")))
    heatmap(as.matrix(metadata(ref_data)$scmap_cluster_index))
    dev.off()

    ## projection
    scmapCluster_results <- scmapCluster(
            projection = target_data, 
            index_list = list(
                ref = metadata(ref_data)$scmap_cluster_index
            ),
            threshold = 0  # force to assign
    )
    saveRDS(scmapCluster_results, file.path(result_dir, paste(prefix, "result.RDS", sep="_")))
    return(scmapCluster_results)
}

run_scmap_cell <- function(ref_data, target_data, result_dir, partition_features=NA, 
                           prefix="scmap_cell", n_features=1000, n_clusters=8,
                           seed=0, select_on="train", select_method="Seurat") {
    ### ---- Run scmap using re
    set.seed(seed)

    # set feature symbols
    rowData(ref_data)$feature_symbol <- rownames(ref_data)
    rowData(target_data)$feature_symbol <- rownames(target_data)
    ref_data <- ref_data[!duplicated(rownames(ref_data)), ]  ## deduplicate genes
    target_data <- target_data[!duplicated(rownames(target_data)), ] 

    ## normalize data
    ref_data <- scmap_normdata(ref_data)
    target_data <- scmap_normdata(target_data)

    if (select_on == "train") {
        ref_data <- scmap_feature_selection(ref_data, result_dir, partition_features,
            n_features, n_clusters, select_method)
    } else if (select_on == "test") {
        target_data <- scmap_feature_selection(target_data, result_dir, partition_features,
            n_features, n_clusters, select_method)
    }

    ## plot selected features
    png(file.path(result_dir, paste(prefix, "feature_plot.png", sep="_")))
    ref_data <- selectFeatures(ref_data, n_features=n_features, suppress_plot = FALSE)
    dev.off()

    ## run scmap-cluster
    ref_data <- indexCell(ref_data)
    png(file.path(result_dir, paste(prefix, "idxcell_heatmap.png", sep="_")))
    heatmap(as.matrix(metadata(ref_data)$scmap_cell_index))
    dev.off()

    ## projection
    scmapCell_results <- scmapCell(
        projection = target_data, 
        index_list = list(
            ref = metadata(ref_data)$scmap_cell_index
        )
    )
    ## annotate data according to reference
    scmapCell_clusters <- scmapCell2Cluster(
        scmapCell_results, 
        list(
          as.character(colData(ref_data)$cell_type1)
        ),
        threshold = 0.1
    )

    saveRDS(scmapCell_clusters, file.path(result_dir, paste(prefix, "result.RDS", sep='_')))
    return(scmapCell_clusters)
}

analyze_scmap_result <- function(target_data, scmap_result, result_dir, prefix="scmap", 
                                 run_time=run_time) {
    ### === analyze result from scmap
    # @run_time: execution time
    ### ===
    if (FALSE) { ## need browser and shiny app
        ## Sankey plot
        png(file.path(result_dir, paste(prefix, "sankey.png", sep="_")))
        plot(
          getSankey(
            colData(target_data)$cell_type1, 
            scmap_result$scmap_cluster_labs[,'ref'],
            plot_height = 400
          )
        )
        dev.off()
    }

    ref_labels <- colData(target_data)$cell_type1
    pred_labels <- scmap_result$scmap_cluster_labs[,'ref']
    stopifnot(length(ref_labels) == length(pred_labels))

    unassigned_cnt <- sum(pred_labels == "unassigned")
    cat("Number of unassigned cells:", unassigned_cnt, 
            "with total cells as", length(pred_labels), '\n')
    if (unassigned_cnt >= length(pred_labels)*0.1) {
        print("Too many unassigned cells!!")
    }

    ## evaluation
    confusion_mat <- table(ref_labels, pred_labels)
    confusion_mat <- as.data.frame.matrix(confusion_mat)
    if ("unassigned" %in% colnames(confusion_mat)) {
        confusion_mat$unassigned <- NULL
    }

    ## when a certain celltype is missing
    if (dim(confusion_mat)[2] < dim(confusion_mat)[1]){
        missing_celltypes <- setdiff(rownames(confusion_mat), colnames(confusion_mat))
        for (celltype in missing_celltypes) {
            confusion_mat[, celltype] <- 0
        }
    }
    confusion_mat <- confusion_mat[, rownames(confusion_mat)]
    confusion_mat <- as.matrix(confusion_mat)

    accuracy <- sum(diag(confusion_mat))/length(ref_labels)
    precision <- diag(confusion_mat)/colSums(confusion_mat)
    recall <- diag(confusion_mat)/rowSums(confusion_mat)
    f1 <- 2 * precision * recall / (precision + recall) 
    f1[is.na(f1)] <- 0
    macroF1 <- mean(f1)

    ARI <- ARI(ref_labels, pred_labels)
    cat("=== Test Accuracy:", accuracy, "Test ARI:", ARI,
        "Test macroF1:", macroF1, "===\n\n")
    fileConn <- file(file.path(result_dir, paste(prefix, "metrics.txt", sep='_')))
    txt <- c(paste0("Acc:", accuracy), paste0("ARI:", ARI), 
             paste0("macroF1:", macroF1), paste0("runtime:", run_time))
    writeLines(txt, fileConn)
    close(fileConn)
}


