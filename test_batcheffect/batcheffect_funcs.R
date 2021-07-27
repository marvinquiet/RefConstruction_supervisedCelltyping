suppressMessages(library(scRNAseqMethods))  ## import Dr. Wu's package
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(Matrix))
suppressMessages(library(umap))
suppressMessages(library(ggplot2))

## batch effect removal
Harmony_removal <- function(X, batchID, result_dir, plot=T) {
    ### ==== Harmony to remove batch effect
    # Harmony can only return a PCA dimension reduction corrected values
    ### ==== 
    require(harmony)
    Xnorm = normalization(X, method="total")

    library(irlba)
    pca_res = prcomp_irlba(Xnorm, n=20)
    pca_projected = pca_res$rotation %*% diag(pca_res$sdev)
    Xpca_umap = umap(pca_projected)
    Xpca_umap_df = as.data.frame(Xpca_umap$layout)
    colnames(Xpca_umap_df) = c("umap1", "umap2")
    rownames(Xpca_umap_df) = colnames(X)
    Xpca_umap_df$batch = as.factor(batchID$batch)

    if (plot) {
        ggplot(Xpca_umap_df, aes(x=umap1, y=umap2, color=batch)) + 
            geom_point(alpha=0.3)
        ggsave(file.path(result_dir, "UMAP_before_removal.png"))
    }

    Xcorrected = HarmonyMatrix(pca_projected, Xpca_umap_df$batch, do_pca = FALSE)

    Xcorrected_umap = umap(Xcorrected)
    Xcorrected_umap_df = as.data.frame(Xcorrected_umap$layout)
    colnames(Xcorrected_umap_df) = c("umap1", "umap2")
    rownames(Xcorrected_umap_df) = colnames(X)
    Xcorrected_umap_df$batch = Xpca_umap_df$batch
    
    if (plot) {
        ggplot(Xcorrected_umap_df, aes(x=umap1, y=umap2, color=batch)) + 
            geom_point(alpha=0.3)
        ggsave(file.path(result_dir, "UMAP_after_removal.png"))
    }
    rownames(Xcorrected) = rownames(Xpca_umap_df)
    write.csv(Xcorrected, file.path(result_dir, "corrected_counts.csv"), quote=F)
}

MNN_removal = function(X, batchID, result_dir, plot=T) {
    require(batchelor)
    #mnn_res = mnnCorrect(combined_data, batch=batchID)
    fastmnn_res = fastMNN(X, batch=batchID$batch)
    fastmnn_corrected = assays(fastmnn_res)$reconstructed

    library(irlba)
    fastmnn_corrected_pca = prcomp_irlba(as.matrix(t(fastmnn_corrected)), n=20)
    fastmnn_corrected_umap = umap(fastmnn_corrected_pca$x)
    fastmnn_corrected_umap_df = as.data.frame(fastmnn_corrected_umap$layout)
    colnames(fastmnn_corrected_umap_df) = c("umap1", "umap2")
    rownames(fastmnn_corrected_umap_df) = colnames(X)
    fastmnn_corrected_umap_df$batch = as.factor(batchID$batch)

    if (plot) {
        ggplot(fastmnn_corrected_umap_df, aes(x=umap1, y=umap2, color=batch)) +
                    geom_point(alpha=0.3)
        ggsave(file.path(result_dir, "UMAP_after_removal.png"))
    }
    write.csv(t(fastmnn_corrected), file.path(result_dir, "corrected_counts.csv"), quote=F)
}
