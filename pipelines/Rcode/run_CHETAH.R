suppressMessages(library(CHETAH))
suppressMessages(library(Rtsne))
suppressMessages(library(aricode))

## refer to Kenong's code on normalizing the data
CHETAH_normdata <- function(sce) {
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

CHETAH_feature_selection <- function(sce, result_dir, partition_features=NA, 
            n_features=1000, n_clusters=NA, select_method="Seurat") {
    ### ---- Feature selection
    # @sce: sce object from CHETAH
    # @result_dir: load features already selected, or store selected features
    # @partition_features: the features provided by partition in GEDFN
    # @n_features: how many features will be selected
    # @select_method: Seurat, FEAST_SC3, FEAST_F-test

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
            Idents(seuratObj) <- "celltypes"
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
                                as.character(colData(sce)$celltypes))$F_scores
 
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

    return(features)
}

run_CHETAH <- function(ref_data, target_data, result_dir, partition_features=NA, 
                       prefix="CHETAH", n_features=1000, n_clusters=8,
                       select_on="train", select_method="Seurat", seed=0) {
    ### ---- Run CHETAH using ref_data as reference and predict cell types on target
    # @ref_data: reference SCE object
    # @target_data: target SCE object
    # @parrition: the partition from GEDFN
    # @result_dir: directory to store the result
    # @prefix: storing result

    ### ----
 
    set.seed(seed)
    ref_data <- CHETAH_normdata(ref_data)
    assay(ref_data, "counts") <- logcounts(ref_data)
    target_data <- CHETAH_normdata(target_data)

    if (select_on == "train") {
        features <- CHETAH_feature_selection(ref_data, result_dir, partition_features,
                                             n_features, n_clusters, select_method)
        marker_genes <- Reduce(intersect, list(rownames(ref_data), rownames(target_data), features))
    } else if (select_on == "test") {
        features <- CHETAH_feature_selection(target_data, result_dir, partition_features,
                                             n_features, n_clusters, select_method)
        marker_genes <- Reduce(intersect, list(rownames(ref_data), rownames(target_data), features))
    } else if (select_on == "NA") {
        features_file <- file.path(result_dir, "features.txt")
        features <- NA
        ## if features have already been selected
        if (file.exists(features_file)) {
            features <- scan(features_file, what=character())
        } else {
            features <- intersect(rownames(ref_data), rownames(target_data))
        }
        ## for intersecting genes
        if (!any(is.na(partition_features))) {
            features <- intersect(partition_features, features)
        }
        marker_genes <- intersect(rownames(ref_data), features)
        marker_genes <- intersect(rownames(target_data), marker_genes)
    }

    cat("Features selected:", length(marker_genes), "\n")

    ref_data <- ref_data[marker_genes, ]
    target_data <- target_data[marker_genes, ]
    target_lognorm <- logcounts(target_data)

    ## dimension reduction on the target data, refer to Kenong's to modify a bit
    pca_data <- prcomp(t(target_lognorm), rank=50)
    tsne_data <- Rtsne(pca_data$x[,1:50], pca=FALSE, perplexity=10,
                       check_duplicates=FALSE)
    tsne_mat <- tsne_data$Y
    dimnames(tsne_mat) <- list(colnames(target_data), c("tSNE_1", "tSNE_2"))
    reducedDims(target_data) <- list(TSNE=tsne_mat)

    ## run CHETAH
    ind <- apply(t(counts(ref_data)), 1, var) == 0 # remove cells with sd of zero
    ref_data <- ref_data[, !ind]

    res <- CHETAHclassifier(input = target_data, ref_cells = ref_data)
    res <- Classify(input = res, 0) ## make the threshold
    #saveRDS(res, file.path(result_dir, paste0(prefix, "_result.RDS")))
    return(res)
}

subsample_cells <- function(ref_data, n_cells=200, seed=0) {
    ## --- Subsample cells in reference dataset as taught in the original website
    set.seed(seed)
    cell_selection <- unlist(lapply(unique(ref$celltypes), function(type) {
        type_cells <- which(ref$celltypes == type)
        if (length(type_cells) > n_cells) {
            type_cells[sample(length(type_cells), n_cells)]
        } else type_cells
    }))
    ref_data <- ref_data[, cell_selection]
    return(ref_data)
}


analyze_CHETAH_result <- function(target_data, CHETAH_result, result_dir, 
                                  prefix="CHETAH", run_time=run_time) {
    ref_labels <- colData(target_data)$celltypes
    pred_labels <- unname(CHETAH_result$celltype_CHETAH)
    stopifnot(length(ref_labels) == length(pred_labels))

    png(file.path(result_dir, paste0(prefix, "_tree.png")), height=500, width=1000)
    PlotCHETAH(input = CHETAH_result)
    dev.off()

    ## evaluation
    confusion_mat <- table(ref_labels, pred_labels)
    confusion_mat <- as.data.frame.matrix(confusion_mat)
    ## when a certain celltype is missing
    union_celltypes <- union(rownames(confusion_mat), colnames(confusion_mat))
    if (dim(confusion_mat)[1] < length(union_celltypes) |
        dim(confusion_mat)[2] < length(union_celltypes)) {
        missing_row_celltypes <- setdiff(union_celltypes, rownames(confusion_mat))
        for (celltype in missing_row_celltypes) {
            confusion_mat[celltype, ] <- 0
        }

        missing_col_celltypes <- setdiff(union_celltypes, colnames(confusion_mat))
        for (celltype in missing_col_celltypes) {
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
