source("test_batcheffect/batcheffect_funcs.R")

args = commandArgs(trailingOnly = TRUE)
method = args[1]
result_dir = args[2]

X = t(read.table(file.path(result_dir,"adata_count.csv"), sep=',', 
               header=T, row.names=1, check.names=F))
batchID = read.table(file.path(result_dir, "batchID.csv"), sep=',',
               header=T, row.names=1, check.names=F)

if (method == "Harmony") {
    Harmony_removal(X, batchID, result_dir, plot=T)
} else if (method == "fastMNN") {
    MNN_removal(X, batchID, result_dir, plot=T)
}
