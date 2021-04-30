##### CHETAH and SCMAP originally written by @KenongSu
library(scmap)
library(CHETAH)

CHETAH_map = function(ref_sce, input_sce, top = 2000){
    # ref_sce includes the counts and celltypes as colData
    ref_col_data = colData(ref_sce)
    if ("celltypes" %in% names(ref_col_data) &  !is.null(assay(ref_sce, "counts"))){
        # input normalization
        input_Y = assay(input_sce, "counts")
        input_L = colSums(input_Y)/median(colSums(input_Y))  ##TODO: using median is more adaptive to datasets
        input_Ynorm = log2(sweep(input_Y, 2, input_L, FUN = "/") + 1)
        # reference normalization
        ref_Y = assay(ref_sce, "counts")
        ref_celltypes = colData(ref_sce)$celltypes
        ref_L = colSums(ref_Y)/median(colSums(ref_Y))
        ref_Ynorm = log2(sweep(ref_Y, 2, ref_L, FUN = "/") + 1)
        # select marker genes from the reference data
        if (top){
            F_res  = cal_F2(ref_Ynorm, classes = ref_celltypes)
            ixs = order(F_res$F_scores, decreasing = T)[1:top]
        }else{
            row_means = rowMeans(Ynorm)
            ixs = order(row_means, decreasing = T)[1:5000]
        }
        ref_genes  = rownames(ref_Y)[ixs]
        marker_genes = intersect(ref_genes, rownames(input_sce))
        ref_Y = ref_Y[marker_genes,]
        input_Y = input_Y[marker_genes, ]
        input_Ynorm = input_Ynorm[marker_genes, ]
        pca_data = prcomp(t(input_Ynorm), rank=50)
        set.seed(123)
        tsne_data = Rtsne(pca_data$x[,1:50], pca = FALSE, perplexity = 10)
        tsne_mat = tsne_data$Y
        dimnames(tsne_mat) = list(colnames(input_sce), c("tSNE_1", "tSNE_2"))
        input_sce = input_sce[marker_genes, ]
        reducedDims(input_sce)  = list(TSNE=tsne_mat)
        # perform the CHETAH mapping
        ref_sce = ref_sce[marker_genes, ]
        res = CHETAHclassifier(input = input_sce, ref_cells = ref_sce)
    }else{
        stop("input the right sce object")
    }
    return(res)
}


scMap_clustermap = function(ref_sce, input_sce, norm_method="norm"){
    ref_col_data = colData(ref_sce)
    ref_row_data = rowData(ref_sce)
    if (!"cell_type1" %in% names(ref_col_data)){
        ref_col_data$cell_type1 = ref_col_data$celltypes
        colData(ref_sce) = ref_col_data
    }
    # process the reference data
    if (is.null(assay(ref_sce, "counts"))){
        stop("input reference data need to have count matrix")
    }else{
        Y = assay(ref_sce, "counts")
        L = colSums(Y) / median(colSums(Y))
        Ynorm = log2(sweep(Y, 2, L, FUN = "/") + 1)
        assay(ref_sce, "logcounts") = Ynorm
        if (! "feature_symbol" %in% names(ref_row_data)){
            ref_row_data$feature_symbol = rownames(ref_sce)
            rowData(ref_sce) = ref_row_data
        }
    }
    # process the input data
    if (is.null(assay(input_sce, "counts"))){
        stop("input input data need to have count matrix")
    }else{
        input_Y = assay(input_sce, "counts")
        input_L = colSums(input_Y) / median(colSums(input_Y))
        input_Ynorm = log2(sweep(input_Y, 2, input_L, FUN = "/") + 1)
        assay(input_sce, "logcounts") = input_Ynorm
        if (is.null(rowData(input_sce)$feature_symbol)){
            input_row_data = rowData(input_sce)
            input_row_data$feature_symbol = rownames(input_sce)
            rowData(input_sce) = input_row_data
        }
    }
    # select features in reference data
    ref_sce = selectFeatures(ref_sce, suppress_plot = T)
    ref_sce = indexCluster(ref_sce)
    ref_mat = metadata(ref_sce)$scmap_cluster_index
    # scmap core
    res = scmapCluster(
        projection = input_sce, 
        index_list = list(ref =ref_mat))
    return(res)
}






