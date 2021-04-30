## load library for reading Anndata H5AD file and selecting top features
suppressMessages(library(FEAST))

args <- commandArgs(trailingOnly = TRUE)
df_dir <- args[1]
annotation_dir <- args[2]
n_features <- args[3]

counts <- read.csv(df_dir, row.names=1)
processed_counts <- process_Y(counts, thre=0)
rm(counts) # release memory

#cell_annotations <- scan(annotation_dir, what=character())
cell_annotations <- readLines(annotation_dir)

F_res <- cal_F2(processed_counts, cell_annotations)$F_scores

names(F_res) <- rownames(processed_counts)
ixs <- order(F_res, decreasing=T)[1:n_features]
features <- names(F_res[ixs])
write(features, file.path(dirname(df_dir), "F-test_features.txt"))
