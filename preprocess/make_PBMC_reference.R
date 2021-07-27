suppressMessages(library(Matrix))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(dplyr))

## make reference for batch2 control/stimulated
PBMC_batch2_ref <- function(data_dir, condition="control", method="scmap") {
    ### ======== build reference panel for batch2 data
    # @ data_dir: data directory i.e: "/home/wma36/gpu/data"
    # @ condition: control or stimulated
    # @ method: scmap or CHETAH
    ### ========

    batch_data_dir <- file.path(data_dir, "PBMC_demuxlet")

    batch2_df <- read.table(file.path(batch_data_dir, "GSE96583_batch2.total.tsne.df.tsv"),
                            header=TRUE, row.names=1, sep="\t")
    idx <- which(batch2_df$multiplets == "singlet" & !is.na(batch2_df$cell))
    batch2_df <- batch2_df[idx, ]

    ## control data
    if (condition == "control") { 
        data <- readMM(file.path(batch_data_dir, "GSM2560248_2.1.mtx"))
        barcodes <- scan(file.path(batch_data_dir, "GSM2560248_barcodes.tsv"), 
                         what=character())
    } else if (condition == "stimulated") {
        data <- readMM(file.path(batch_data_dir, "GSM2560249_2.2.mtx"))
        barcodes <- scan(file.path(batch_data_dir, "GSM2560249_barcodes.tsv"), 
                         what=character())
    }
    genes_df <- read.table(file.path(batch_data_dir, "GSE96583_batch2.genes.tsv"),
                        header=FALSE, sep="\t")
    genes <- genes_df$V2
    rownames(data) <- genes
    colnames(data) <- barcodes

    common_cells <- intersect(rownames(batch2_df), barcodes)
    ## filter matrix by common cells
    data <- data[genes, common_cells]
    rownames(data) <- toupper(rownames(data))  ## convert to uppercase
    if (method == "CHETAH")
        colData <- DataFrame(celltypes = batch2_df[common_cells, ]$cell)
    if (method == "scmap")
        colData <- DataFrame(cell_type1 = batch2_df[common_cells, ]$cell)
    rownames(colData) <- common_cells

    ref <- SingleCellExperiment(assays = list(counts = as.matrix(data)),
                                colData = colData)
    return(ref)
}

## make reference for batch1 A/B/C
PBMC_batch1_ref <- function(data_dir, sample="A", method="scmap") {
    ### ======== build reference panel for batch1 data
    # @ data_dir: data directory
    # @ sample: A/B/C
    # @ method: scmap or CHETAH
    ### ========
    batch_data_dir <- file.path(data_dir, "PBMC_demuxlet")

    batch1_df <- read.table(file.path(batch_data_dir, "GSE96583_batch1.total.tsne.df.tsv"),
                            header=TRUE, row.names=1, sep="\t")
    idx <- which(batch1_df$multiplets == "singlet" & !is.na(batch1_df$cell.type))
    batch1_df <- batch1_df[idx, ]

    ## control data
    if (sample == "A") { 
        data <- readMM(file.path(batch_data_dir, "GSM2560245_A.mat"))
        barcodes <- scan(file.path(batch_data_dir, "GSM2560245_barcodes.tsv"), 
                         what=character())
    } else if (sample == "B") {
        data <- readMM(file.path(batch_data_dir, "GSM2560246_B.mat"))
        barcodes <- scan(file.path(batch_data_dir, "GSM2560246_barcodes.tsv"), 
                         what=character())
    } else if (sample == "C") {
        data <- readMM(file.path(batch_data_dir, "GSM2560247_C.mat"))
        barcodes <- scan(file.path(batch_data_dir, "GSM2560247_barcodes.tsv"), 
                         what=character())
    }
    genes_df <- read.table(file.path(batch_data_dir, "GSE96583_batch1.genes.tsv"),
                        header=FALSE, sep="\t")
    genes <- genes_df$V2
    rownames(data) <- genes
    colnames(data) <- barcodes

    common_cells <- intersect(rownames(batch1_df), barcodes)
    ## filter matrix by common cells
    data <- data[genes, common_cells]
    rownames(data) <- toupper(rownames(data))  ## convert to uppercase
    if (method == "CHETAH")
        colData <- DataFrame(celltypes = batch1_df[common_cells, ]$cell.type)
    if (method == "scmap"){
        colData <- DataFrame(cell_type1 = batch1_df[common_cells, ]$cell.type)
    }
    rownames(colData) <- common_cells
    ref <- SingleCellExperiment(assays = list(counts = as.matrix(data)),
                                colData = colData)
    return(ref) 
}

## make reference for batch1 S1 and S5
PBMC_batch1_ind_ref <- function(data_dir, ind=NA, method="scmap") {
    ### ======== build reference panel for batch1 data
    # @ data_dir: data directory
    # @ ind: individual number; S1: 1154 in batch A, S2: 1085 in batch B
    # @ method: scmap or CHETAH
    ### ========
    batch_data_dir <- file.path(data_dir, "PBMC_demuxlet")

    batch1_df <- read.table(file.path(batch_data_dir, "GSE96583_batch1.total.tsne.df.tsv"),
                            header=TRUE, row.names=1, sep="\t")
    A_ref <- PBMC_batch1_ref(data_dir, sample="A", method=method)
    B_ref <- PBMC_batch1_ref(data_dir, sample="B", method=method)
    C_ref <- PBMC_batch1_ref(data_dir, sample="C", method=method)

    combined_ref <- cbind(A_ref, B_ref, C_ref) ## dim: 32738 11792
    if (is.na(ind)) {
        sample_cells = intersect(rownames(batch1_df), colnames(combined_ref))
    } else {
        idx <- which(batch1_df$ind == ind)
        ## get sample cells
        sample_cells <- intersect(rownames(batch1_df)[idx], colnames(combined_ref))
    }
    ref <- combined_ref[, sample_cells]
    return(ref) 
}

## make reference for PBMC 7 protocols
PBMC_protocols_ref <- function(data_dir, method="scmap", exp=NA, 
                               protocol=NA, protocol_type=NA, curate=F) {
    ### ==== build reference panel for PBMC 7 protocols
    # @ data_dir: data directory
    # @ method: scmap or CHETAH
    ### ====

    plate_protocols <- c('CEL-Seq2', 'Smart-seq2')
    metadata_df <- read.table(file.path(data_dir, "PBMC_protocols/metadata.txt"),
                            header=T, sep="\t")
    ## plate-based data
    plateMM <- readMM(file.path(data_dir, "PBMC_protocols/counts.read.txt"))
    platecells <- scan(file.path(data_dir, "PBMC_protocols/cells.read.new.txt"),
                       what=character())
    colnames(plateMM) <- platecells
    plate_metadata <- metadata_df[metadata_df$Method %in% plate_protocols, ]
    rownames(plate_metadata) <- plate_metadata$NAME
    common_cells <- intersect(rownames(plate_metadata), colnames(plateMM))
    plateMM <- plateMM[, common_cells]
    plategenes <- scan(file.path(data_dir, "PBMC_protocols/genes.read.txt"),
                       what=character())
    rownames(plateMM) <- plategenes
    gene_symbols <- sapply(strsplit(plategenes, '_'), '[', 2)
    uniq_gene_ID <- match(unique(gene_symbols), gene_symbols)
    plateMM <- plateMM[uniq_gene_ID, ]
    rownames(plateMM) <- sapply(strsplit(rownames(plateMM), '_'), '[', 2)
    rownames(plateMM) <- toupper(rownames(plateMM))  ## convert to uppercase
    if (method == "CHETAH")
        colData <- DataFrame(celltypes = plate_metadata[colnames(plateMM), ]$CellType)
    if (method == "scmap"){
        colData <- DataFrame(cell_type1 = plate_metadata[colnames(plateMM), ]$CellType)
    }
    rownames(colData) <- colnames(plateMM)
    plate_ref <- SingleCellExperiment(assays = list(counts = as.matrix(plateMM)),
                                colData = colData)
    merge_df <- merge(
            x = cbind(rownames=rownames(colData(plate_ref)), colData(plate_ref)),
            y = plate_metadata,
            all.x = TRUE,
            by.x="rownames", by.y="NAME")
    rownames(merge_df) <- merge_df$rownames
    merge_df$rownames <- NULL
    colData(plate_ref) <- merge_df[colnames(plate_ref), ]

    ## umi-based data
    umiMM <- readMM(file.path(data_dir, "PBMC_protocols/counts.umi.txt"))
    umicells <- scan(file.path(data_dir, "PBMC_protocols/cells.umi.new.txt"),
                       what=character())
    colnames(umiMM) <- umicells
    umi_metadata <- metadata_df[!metadata_df$Method %in% plate_protocols, ]
    rownames(umi_metadata) <- umi_metadata$NAME
    common_cells <- intersect(rownames(umi_metadata), colnames(umiMM))
    umiMM <- umiMM[, common_cells]
    umigenes <- scan(file.path(data_dir, "PBMC_protocols/genes.read.txt"),
                       what=character())
    rownames(umiMM) <- umigenes
    gene_symbols <- sapply(strsplit(umigenes, '_'), '[', 2)
    uniq_gene_ID <- match(unique(gene_symbols), gene_symbols)
    umiMM <- umiMM[uniq_gene_ID, ]
    rownames(umiMM) <- sapply(strsplit(rownames(umiMM), '_'), '[', 2)
    rownames(umiMM) <- toupper(rownames(umiMM))  ## convert to uppercase
    if (method == "CHETAH")
        colData <- DataFrame(celltypes = umi_metadata[colnames(umiMM), ]$CellType)
    if (method == "scmap"){
        colData <- DataFrame(cell_type1 = umi_metadata[colnames(umiMM), ]$CellType)
    }
    rownames(colData) <- colnames(umiMM)
    umi_ref <- SingleCellExperiment(assays = list(counts = as.matrix(umiMM)),
                                colData = colData)
    merge_df <- merge(
            x = cbind(rownames=rownames(colData(umi_ref)), colData(umi_ref)),
            y = umi_metadata,
            all.x = TRUE,
            by.x="rownames", by.y="NAME")
    rownames(merge_df) <- merge_df$rownames
    merge_df$rownames <- NULL
    colData(umi_ref) <- merge_df[colnames(umi_ref), ]

    ## combine two reference (same genes: all(rownames(plate_ref) == rownames(umi_ref)) == TRUE)
    ref <- cbind(plate_ref, umi_ref)
    ## replace 10X method name
    ref_colData <- colData(ref)
    ref_colData[, 'Method'] <- as.character(ref_colData[, 'Method'])
    idx <- which(ref_colData[, 'Method'] == "10x Chromium (v2)" |
                 ref_colData[, 'Method'] == "10x Chromium (v2) A" |
                 ref_colData[, 'Method'] == "10x Chromium (v2) B")
    ref_colData[idx, 'Method'] <- "10x-v2"

    idx <- which(ref_colData[, 'Method'] == "10x Chromium (v3)")
    ref_colData[idx, 'Method'] <- "10x-v3"
    ref_colData[, 'Method'] <- as.factor(ref_colData[, 'Method'])
    if (curate) { ## curate the cell types
        celltype_ind = ifelse(method == "scmap", "cell_type1", "celltypes")
        ref_colData[, celltype_ind] = as.character(ref_colData[, celltype_ind])
        idx = which(ref_colData[, celltype_ind] == "B cell")
        ref_colData[idx, celltype_ind] = "B cells"
        idx = which(ref_colData[, celltype_ind] == "CD14+ monocyte")
        ref_colData[idx, celltype_ind] = "CD14+ Monocytes"
        idx = which(ref_colData[, celltype_ind] == "CD4+ T cell")
        ref_colData[idx, celltype_ind] = "CD4 T cells"
        idx = which(ref_colData[, celltype_ind] == "Cytotoxic T cell")
        ref_colData[idx, celltype_ind] = "CD8 T cells"
        idx = which(ref_colData[, celltype_ind] == "Natural killer cell")
        ref_colData[idx, celltype_ind] = "NK cells"
        ref_colData[, celltype_ind] = as.factor(ref_colData[, celltype_ind])
    }
    colData(ref) <- ref_colData

    if (!is.na(exp)) {
        exp_idx = colData(ref)$Experiment == exp
        ref = ref[, exp_idx]
    }

    if (!is.na(protocol)) {
        proc_idx = colData(ref)$Method == protocol
        ref = ref[, proc_idx]
    }

    if (!is.na(protocol_type)) {
        if (protocol_type == "plate") {
            proc_type_idx = colData(ref)$Method %in% c("Smart-seq2", "CEL-Seq2")
        } else if (protocol_type == "droplet") {
            proc_type_idx = colData(ref)$Method %in% c("10x-v2", "10x-v3", "Drop-seq", "inDrops", "Seq-Well")
        }
        ref = ref[, proc_type_idx]
    }
    return(ref)
}

## make reference for PBMC Zheng FACS
PBMC_Zheng_ref <- function(data_dir, method="scmap", curate=FALSE) {
    ### ==== build reference panel for PBMC Zheng FACS-sorted
    # @ data_dir: data directory
    # @ method: scmap or CHETAH
    # @ curate: whether to curate the cell types
    ### ====
    Zheng_data_dir <- file.path(data_dir, "PBMC_Zheng_FACS")
    celltypes = c("b_cells", "cd14_monocytes", "cd34", "cd56_nk",
            "cd4_t_helper", "naive_t", "memory_t", "regulatory_t", ## CD4+ T cells
            "cytotoxic_t", "naive_cytotoxic") ## CD8+ T cells
    data = NULL
    celltype_col = c()
    for (celltype in celltypes) {
        celltype_dir = file.path(Zheng_data_dir, celltype)
        if (curate) {
            if (celltype == "cd34") next ## do not bind cd34 cell
            if (celltype %in% c("cd4_t_helper", "naive_t", "memory_t", "regulatory_t")) 
                celltype = "CD4 T cells"
            if (celltype %in% c("cytotoxic_t", "naive_cytotoxic"))
                celltype = "CD8 T cells"
            if (celltype == "b_cells")
                celltype = "B cells"
            if (celltype == "cd14_monocytes")
                celltype = "CD14+ Monocytes"
            if (celltype == "cd56_nk")
                celltype = "NK cells"
        }

        celltype_data = readMM(file.path(celltype_dir, "matrix.mtx"))
        celltype_genes = read.table(file.path(celltype_dir, "genes.tsv"), 
                                    header=FALSE, sep='\t')
        dedup_genes_idx = which(!duplicated(celltype_genes$V2))
        celltype_data = celltype_data[dedup_genes_idx, ]
        rownames(celltype_data) = as.character(celltype_genes[dedup_genes_idx, ]$V2)
        celltype_barcodes = scan(file.path(celltype_dir, "barcodes.tsv"), 
                                 what=character())
        colnames(celltype_data) = paste(celltype_barcodes, celltype, sep='-')
        celltype_col = c(celltype_col, rep(celltype, ncol(celltype_data)))

        if (is.null(data))
            data = celltype_data
        else {
            data = cbind(data, celltype_data[rownames(data), ])
        }
    }

    if (method == "CHETAH")
        colData = DataFrame(celltypes=celltype_col)
    if (method == "scmap")
        colData = DataFrame(cell_type1=celltype_col)
    rownames(colData) = colnames(data)

    ref = SingleCellExperiment(assays = list(counts=data),
                               colData = colData)
    return(ref)
}
