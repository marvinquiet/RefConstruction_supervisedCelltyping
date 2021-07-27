suppressMessages(library(Matrix))
suppressMessages(library(SingleCellExperiment))

ref_target_common_cells <- function(ref_data, target_data, method="scmap") {
    ### === curate for same cell types
    # @ref_data: sce object
    # @target_data: sce object
    # @method: scmap or CHETAH -> used to indicate colData celltype
    ### ===
    celltype_ind <- ifelse("scmap" == method, "cell_type1", "celltypes")
    common_celltypes <- intersect(colData(ref_data)[, celltype_ind], colData(target_data)[, celltype_ind])
    ref_idx <- colData(ref_data)[, celltype_ind] %in% common_celltypes
    target_idx <- colData(target_data)[, celltype_ind] %in% common_celltypes
    ref_data <- ref_data[, ref_idx]
    target_data <- target_data[, target_idx]
    return(list(ref=ref_data, target=target_data))
}

ref_downsample <- function(ref_data, size, seed=0) {
    ### === downsample reference dataset
    # @ref_data: sce object
    # @size: downsample to the size
    ### ===
    set.seed(seed)
    if (ncol(ref_data) > size) {
        sampled_cells <- sample(colnames(ref_data), size=size)
        ref_data <- ref_data[, sampled_cells]
    }
    return(ref_data)
}

## filter data to express at least in 10 cells with 10 genes
filter_sce = function(sce, min_genes=10, min_cells=10) {
    ### === filter sce object
    # @sce: sce object
    # @min_genes: minimum genes expressed
    # @min_cells: minimum cells expressed
    n_exprsgenes = rowSums(counts(sce) > 0)
    n_exprscells = colSums(counts(sce) > 0)
    selected_genes = which(n_exprsgenes >= 10)
    selected_cells = which(n_exprscells >= 10)
    return (sce[selected_genes, selected_cells])
}
