## load library for reading Anndata H5AD file and selecting top features
suppressMessages(library(FEAST))

args <- commandArgs(trailingOnly = TRUE)
df_dir <- args[1]
n_celltypes <- args[2]
n_features <- args[3]

counts <- read.csv(df_dir, row.names=1)
processed_counts <- process_Y(counts, thre=0)
rm(counts) # release memory

if (ncol(processed_counts) > 5000) {
    F_res_idx <- FEAST_fast(processed_counts, k=n_celltypes, 
                      split=T, batch_size=2000)
} else if (ncol(processed_counts) > 1000) {
    F_res_idx <- FEAST_fast(processed_counts, k=n_celltypes, 
                      split=T, batch_size=500)
} else if (ncol(processed_counts) > 500){
    F_res_idx <- FEAST_fast(processed_counts, k=n_celltypes, 
                      split=T, batch_size=200, num_cores=4)
} else {
    F_res_idx <- FEAST_fast(processed_counts, k=n_celltypes, 
                      split=F, num_cores=4)
}

#names(F_res) <- rownames(processed_counts)
#ixs <- order(F_res, decreasing=T)[1:n_features]
#features <- names(F_res[ixs])
features = rownames(processed_counts)[F_res_idx[1:n_features]]
write(features, file.path(dirname(df_dir), "FEAST_features.txt"))
