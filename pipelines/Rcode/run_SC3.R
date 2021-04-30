suppressMessages(library(SC3))
suppressMessages(library(scater))
suppressMessages(library(aricode))

run_SC3_analysis <- function(sce, n_clusters, result_dir, prefix="sc3",
                             biology=FALSE, n_cores=4, dr_seed=0, 
                             dataset="seg", method="scmap") {
    ### ==== Run SC3 analysis on the dataset
    # @sce: SingleCellExperiment Object
    # @n_clusters: number of celltypes
    # @result_dir: result directory
    # @biology: SC3 biological features based on the identified cell clusters
    # @n_cores: number of cores for running SC3
    # @dr_seed: dimension reduction seed
    ### ====
    ## normalized data
    counts <- assay(sce, "counts")
    lib.sizes <- colSums(counts)
    scale.factor <- 1e4
    normcounts(sce) <- scale(counts, center=FALSE, scale=lib.sizes) * scale.factor
    logcounts(sce) <- log2(normcounts(sce) + 1) ## log-normalize the data

    rowData(sce)$feature_symbol <- rownames(sce)
    sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]

    features_file <- file.path(result_dir, "features.txt")
    if (!file.exists(features_file)) {
        print("Please make sure features.txt exist in the given directory!")
        q()
    }
    features <- scan(features_file, what=character())
    common_features <- intersect(features, rownames(sce))
    sce <- sce[common_features, ] ## feature selection

    ## dimension reduction
    set.seed(dr_seed)
    sce <- runPCA(sce)
    sce <- runTSNE(sce)

    sce <- sc3(sce, ks=n_clusters, biology=biology, n_cores=n_cores)
    colname <- paste('sc3', n_clusters, 'clusters', sep='_')
    cat("===SC3 ARI: ", ARI(colData(sce)[, "cell_type1"], colData(sce)[, colname]), "\n")

    png(file.path(result_dir, paste(prefix, "PCA.png", sep='_')))
    plotPCA(sce, colour_by="cell_type1")
    dev.off()

    png(file.path(result_dir, paste(prefix, "celltype_tSNE.png", sep='_')))
    plotTSNE(sce, colour_by="cell_type1")
    dev.off()

    png(file.path(result_dir, paste(prefix, "sc3cluster_tSNE.png", sep='_')))
    plotTSNE(sce, colour_by=colname)
    dev.off()

    png(file.path(result_dir, paste(prefix, "ind_tSNE.png", sep='_')))
    plotTSNE(sce, colour_by="ind")
    dev.off()

    png(file.path(result_dir, paste(prefix, "condition_tSNE.png", sep='_')))
    plotTSNE(sce, colour_by="condition")
    dev.off()

}
