



### deMULTIplex2 ###

# Objectives:
# 1. Include functionality for MULTI-ATAC barcode sequencing
# 2. Improve demultiplexing algorithm
# 3. Improve plotting features
# 4. Include data-scaling options

# full list of 96 barcodes
# barcodes <- readxl::read_excel("/Users/dannyconrad/Box/Data/deMULTIplex2/MULTI-seq_oligos_Nov2019.xlsx")$`Barcode Sequence`



# first 16 barcodes, to move to data
barcodes <- c("GGAGAAGA","CCACAATG","TGAGACCT","GCACACGC", #1-4
              "AGAGAGAG","TCACAGCA","GAAAAGGG","CGAGATTC", #5-8
              "GTAGCACT","CGACCAGC","TTAGCCAG","GGACCCCA", #9-12
              "CCAACCGG","TGACCGAT","GCAACGCC","CAATCGGT") #13-16



### Key functions ###

#################
## revComplement() ##
#################

# Description:
# Takes cell barcodes formatted by Cell Ranger and/or ArchR and returns reverse complement of each cell barcode sequence
#' @export
revComplement <- function(cells,
                      cbLength = 16,
                      keepFormat = F) {
    require(Biostrings)
    require(stringr)

    pattern <- paste("[(A,G,T,C,N)]{", cbLength, "}", sep = "")

    if (keepFormat) {
        cells_tmp <- paste("x",cells,"x", sep = "")
        format <- str_split(cells_tmp, pattern, simplify = T)
        format <- t(apply(format, 1, function(x) {
            tmp <- x
            str_sub(tmp[1], 1, 1) <- ""
            str_sub(tmp[2], -1, -1) <- ""
            tmp
        }))
    }

    cells <- str_extract(cells, pattern)

    res <- as.character(reverseComplement(DNAStringSet(cells)))

    if (keepFormat) {
        res <- paste(format[,1], res, format[,2], sep = "")
    }
    return(res)
}


## readTag()
# previously MULTIseq.preprocess()
# memory intensive for large FASTQ files

# Description:
# Locates and reads into memory FASTQ files produced by sequencing MULTIseq or MULTI-ATAC barcode libraries.
# Read pairs are parsed to extract 10x cell barcodes, MULTI* sample barcodes, and UMIs to count barcodes captured, and these are assembled into a read table
# Read table can be optionally filtered by a set of known cell barcodes (i.e. cells identified by Cell Ranger)
#' @export
readTag <- function(dir,
                      name,
                      barcode.type = c("MULTIseq", "MULTI-ATAC"),
                      assay = c("RNA", "ATAC", "Multiome"),
                      filter.cells = NULL,
                      fastq.A = NA, fastq.B = NA,
                      pos.cbc, pos.sbc, pos.umi) {

    t0 <- Sys.time()
    barcode.type <- match.arg(barcode.type)
    assay <- match.arg(assay)

    if (barcode.type == "MULTIseq" & assay == "RNA") {
        pattern.A <- "R1" # minimum 28bp read that contains cell BC and UMI
        pattern.B <- "R2" # typically 91bp read that contains sample BC
    } else if (barcode.type == "MULTIseq" & assay == "Multiome") {
        pattern.A <- "R1" # minimum 28bp read that contains cell BC and UMI
        pattern.B <- "R3" # typically 91bp read that contains sample BC
    } else if (barcode.type == "MULTI-ATAC") {
        pattern.A <- "R2" # minimum 16bp "index" read that contains cell BC
        pattern.B <- "R3" # typically > 50bp read that contains sample BC and UMI
    } else {
        return(message("Error - barcode should be either MULTIseq or MULTI-ATAC"))
    }

    # set assay = NULL in order to use custom position indices (must manually define all 3)
    if (barcode.type == "MULTIseq" & assay %in% c("RNA","Multiome")) {
        pos.cbc <- c(1,16)
        pos.sbc <- c(1,8)
        pos.umi <- c(17,28)
    } else if (barcode.type == "MULTI-ATAC" & assay == "ATAC") {
        pos.cbc <- c(1,16)
        pos.sbc <- c(1,8)
        pos.umi <- c(9,16)
    } else if (barcode.type == "MULTI-ATAC" & assay == "Multiome") {
        pos.cbc <- c(9,24)
        pos.sbc <- c(1,8)
        pos.umi <- c(17,28)
    } else { cat("Using custom read positions") }

    if (file.exists(as.character(fastq.A))) {
        path.A <- fastq.A
    } else {
        path.A <- list.files(path = dir, pattern = paste(name, pattern.A, ".fastq", sep=".*"), full.names = T)
    }

    if (file.exists(as.character(fastq.B))) {
        path.B <- fastq.B
    } else {
        path.B <- list.files(path = dir, pattern = paste(name, pattern.B, ".fastq", sep=".*"), full.names = T)
    }

    if (length(path.A) < 1 | length(path.B) < 1) {
        return(message("Error - one or more FASTQ files not found"))
    } else if (length(path.A) > 1 | length(path.B) > 1) {
        return(message("Error - too many files with match names found"))
    }


    cat("### Building Read Table ###", fill = T)
    cat("Barcode Type: ", barcode.type, sep = "", fill = T)
    cat("10x Genomics Assay Type: ", assay, sep = "", fill = T)
    cat("Loading first set of reads...", fill = T)
    r <- readFastq(path.A)
    gc(verbose = F)

    if (barcode.type == "MULTIseq") {
        cat("Extracting Cell Barcodes & UMIs...", fill = T)
        read_table <- data.frame(Cell = subseq(sread(r), pos.cbc[1], pos.cbc[2]),
                                UMI = subseq(sread(r), pos.umi[1], pos.umi[2]))
    }
    if (barcode.type == "MULTI-ATAC") {
        cat("Extracting Cell Barcodes...", fill = T)
        read_table <- data.frame(Cell = subseq(sread(r), pos.cbc[1], pos.cbc[2]))
    }
    cat("Finished processing first set of reads; unloading from memory", fill = T)
    r <- NULL
    gc(verbose = F)

    cat("Loading second set of reads...", fill = T)
    r <- readFastq(path.B)
    gc(verbose = F)

    if (barcode.type == "MULTIseq") {
        cat("Extracting Sample Barcodes...", fill = T)
        read_table$Sample <- as.character(subseq(sread(r), pos.sbc[1], pos.sbc[2]))
    }
    if (barcode.type == "MULTI-ATAC") {
        cat("Extracting Sample Barcodes & UMIs...", fill = T)
        read_table$UMI <- as.character(subseq(sread(r), pos.umi[1], pos.umi[2]))
        read_table$Sample <- as.character(subseq(sread(r), pos.sbc[1], pos.sbc[2]))
    }
    cat("Finished processing second set of reads; unloading from memory", fill = T)
    r <- NULL
    gc(verbose = F)

    cat("Finished parsing ", nrow(read_table), " read pairs", sep = "", fill = T)

    if (is.null(filter.cells)) {
        cat("Keeping all reads because filter.cells = NULL", fill = T)
    } else if (!is.null(filter.cells)) {
        cat("Filtering for ", length(filter.cells), " provided cell barcodes...", sep ="", fill = T)
        ind <- which(read_table$Cell %in% filter.cells)
        read_table <- read_table[ind, ]
    }

    read_table <- read_table[,c("Cell","Sample","UMI")]

    cat("Finished building read table", fill = T)
    cat("\n")
    cat("Finished in",
        round(difftime(Sys.time(), t0, units = "mins")[[1]],1),
        "minutes", sep = " ")
    cat("\n")
    return(read_table)
}






##################
## alignMULTI() ##
##################

# previously MULTIseq.align()

# Description:
# Utilizes data.table library to quickly tally total UMI counts of each sample barcode per cell in read table
#' @export
alignMULTI <- function(read_table,
                       barcodes,
                       filter.cells = NULL,
                       string.dist.method = c("hamming", "osa"),
                       min.dist = 1) {
    t0 <- Sys.time()
    string.dist.method <- match.arg(string.dist.method)
    if (is.null(filter.cells)) {
        cells <- unique(read_table$Cell)
        cells <- cells[cells != paste(rep("G",16),collapse = "")]
    } else { cells <- filter.cells }

    if(min.dist > 0 ){
        cat(paste0("Accounting for mismatch with minimal ", string.dist.method, " distance ", min.dist), fill = T)
        unique_bcs <- unique(read_table$Sample)
        corrected_bcs<-sapply(unique_bcs, function(x){
            if(!x %in% barcodes) {
                x2 = barcodes[which(stringdist(x, barcodes, method = string.dist.method)==min.dist)]
                if(length(x2) == 1) x2 else NA
            } else x
        })
        read_table$Sample <- corrected_bcs[read_table$Sample]
    }

    cat("Assembling Barcode Count Table...", fill = T)
    dt <- data.table(read_table[,c("Cell", "Sample")])
    dt <- dt[complete.cases(dt), ]
    assign("dt", dt, env = .GlobalEnv)
    cnt_ind <- dt[, list(Freq =.N), by=list(Cell,Sample)]
    cnt_ind$i <- (1:length(cells))[match(cnt_ind$Cell, cells)] # NA?
    cnt_ind$j <- (1:length(barcodes))[match(cnt_ind$Sample, barcodes)]
    cnt_ind <- cnt_ind[complete.cases(cnt_ind),] # Added to account for provided cell list not in data
    bc_mtx <- sparseMatrix(i = cnt_ind$i, j = cnt_ind$j, x = cnt_ind$Freq, dims = c(length(cells), length(barcodes)))
    rownames(bc_mtx) <- cells
    colnames(bc_mtx) <- barcodes

    cat("Finished in",
        round(difftime(Sys.time(), t0, units = "secs")[[1]]),
        "seconds", sep = " ")

    cat("\n")
    return(bc_mtx)
}

















