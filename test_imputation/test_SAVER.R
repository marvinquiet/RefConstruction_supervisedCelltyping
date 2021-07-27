suppressMessages(library(SAVER))

## read in arguments
args <- commandArgs(trailingOnly = TRUE)
result_dir <- args[1]

if (is.na(result_dir)) {
    print("Error: Please make sure to have the result directory as input")
}

train_data = t(read.table(file.path(result_dir, 'train_data.csv'), 
                    sep=',', row.names=1, header=T, check.names=FALSE))
test_data = t(read.table(file.path(result_dir, 'test_data.csv'),
                    sep=',', row.names=1, header=T, check.names=FALSE))

train_saver = saver(train_data, ncores=10, size.factor=1, estimates.only=TRUE)
test_saver = saver(test_data, ncores=10, size.factor=1, estimates.only=TRUE)
write.csv(t(train_saver), file.path(result_dir, 'SAVER_train_data.csv'), quote=F)
write.csv(t(test_saver), file.path(result_dir, 'SAVER_test_data.csv'), quote=F)

