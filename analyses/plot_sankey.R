library(networkD3)
library(dplyr)
library(aricode)
library(ggplot2)

### ==== Function of this script
# 1. calculate the ARI between predicted label from FC and HC
# 2. generate the sankey plot and check the relationship between sub-cell types from
#    different brain regions
### ====

data_source <- "mousebrain_region_sub_FC_to_HC"
pipeline_dir <- "/home/wma36/gpu/sc_identifier/pipelines"
result_dirs <- list.files(pipeline_dir, pattern=paste0('*', data_source, '*'), 
                          full.names=TRUE)

plot_sankey <- function(original, predicted, 
                        original_prefix="HC-", predicted_prefix="FC-") {
    ### ==== use networkD3 to generate sankey plot
    # @original: original cell type in the data
    # @predicted: predicted cell type in the data
    # @original_prefix: prefix for original data; 
    # @predicted_prefix: prefix for predicted result
    ### ====
    source_nodes <- paste0(original_prefix, 
                           unique(sapply(strsplit(as.character(original), '\\.'), '[', 2)))
    target_nodes <- paste0(predicted_prefix,
                           unique(sapply(strsplit(as.character(predicted), '\\.'), '[', 2)))
    sankey_nodes <- data.frame(node=seq(0, length(source_nodes)+length(target_nodes)-1),
                               name=c(source_nodes, target_nodes))
    
    source_celltypes <- sapply(strsplit(as.character(original), '\\.'), '[', 2)
    source_nodesID <- sankey_nodes[match(paste0(original_prefix, source_celltypes), sankey_nodes$name),]$node 
    target_celltypes <- sapply(strsplit(as.character(predicted), '\\.'), '[', 2)
    target_nodesID <- sankey_nodes[match(paste0(predicted_prefix, target_celltypes), sankey_nodes$name),]$node 
    sankey_links <- data.frame(source=source_nodesID, target=target_nodesID,
                               value=rep(0.1, nrow(data)))
    network <- networkD3::sankeyNetwork(Links=sankey_links, Nodes=sankey_nodes, 
                         Source='source', 
                         Target='target', 
                         Value='value', 
                         NodeID='name')
    return(network)
}

python_suffix <- "_predicted_obs.csv"
R_suffix <- "_result.RDS"
for (result_dir in result_dirs) {
    print(result_dir)

    ## deal with Python result files
    files <- list.files(result_dir, pattern=paste0('*', python_suffix, '$'))
    for (file in files) {
        method <- gsub(python_suffix, '', file)
        file_dir <- file.path(result_dir, file)
        data <- read.csv(file_dir, header=T, row.names=1)
    
        ## calculate the ARI
        ARI <- ARI(data$cell.type, data$pred_celltypes)
        cat(method, ':', ARI, '\n')
    
        ## plot Sankey on major cell types
        network <- plot_sankey(data$cell.type, data$pred_celltypes)
        saveNetwork(network, file.path(result_dir, 
                    paste(method, "sankey.html", sep='_')),
                    selfcontained=TRUE)

        ## plot heatmap for sub-cell types
        df <- data.frame(source=as.character(data$cell.type), 
                         target=as.character(data$pred_celltypes),
                         value=rep(1, nrow(data)))
        aggr_df <- aggregate(df$value, by=list(source=df$source, target=df$target), FUN=sum)
        count_df <- aggregate(df$value, by=list(source=df$source), FUN=sum)
        merged_df <- merge(aggr_df, count_df, by="source")
        merged_df$percent <- merged_df$x.x/merged_df$x.y*100

        g <- ggplot(data=merged_df, aes(x=source, y=target, fill=percent)) +
             geom_tile() +
             scale_fill_gradient(low = "white", high = "red") + 
             theme(axis.text.x=element_text(size=6, angle=90, vjust=1, hjust=1),
                  axis.text.y=element_text(size=6)) +
             labs(x="HC", y="FC")
        ggsave(file.path(result_dir, paste(method, 'heatmap.png', sep='_')))
    }

    ## deal with R result files
    files <- list.files(result_dir, pattern=paste0('*', R_suffix, '$'))
    for (file in files) {
        method <- gsub(R_suffix, '', file)
        file_dir <- file.path(result_dir, file)

        ## calculate the ARI
        if (method == "CHETAH") {
            chetah_res <- read.table(file.path(result_dir, paste0(method, "_metrics.txt")), 
                                    header=F, sep=":")
            ARI <- chetah_res[chetah_res$V1 == "ARI", ]$V2

            #result <- readRDS(file_dir)
            ### calculate the ARI
            #ARI <- ARI(result$celltypes, result$celltype_CHETAH)
            #
            ### plot sankey
            #network <- plot_sankey(result$celltypes, result$celltype_CHETAH)
            #saveNetwork(network, file.path(result_dir, 
            #        paste(method, "sankey.html", sep='_')),
            #        selfcontained=TRUE)
        }

        if (method == "scmap") {
            scmap_res <- read.table(file.path(result_dir, paste0(method, "_metrics.txt")), 
                                    header=F, sep=":")
            ARI <- scmap_res[scmap_res$V1 == "ARI", ]$V2
        }

        cat(method, ':', ARI, '\n')
    }
}
