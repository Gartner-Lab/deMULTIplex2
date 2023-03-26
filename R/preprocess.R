

################
### readTags ###
################
#' Preprocess Tag Reads
#'
#' Loads the FASTQ files produced by sequencing MULTIseq or MULTI-ATAC barcode libraries. \cr
#' Read pairs are parsed to extract 10X cell barcodes, MULTI* sample tags, and UMIs, and these are assembled into a read table. \cr
#' Read table can be optionally filtered by a set of known cell barcodes (i.e. cells identified by Cell Ranger). \cr
#' Warning: This function is memory intensive for large FASTQ files. \cr
#' Note: This is an updated version of MULTIseq.preprocess() from original deMULTIplex package.
#'
#' @param dir Path to directory containing tag FASTQ files
#' @param name Prefix of FASTQ files for the tag library to be processed
#' @param barcode.type Type of barcoding used: "MULTIseq" or "MULTI-ATAC"
#' @param assay 10X Genomics single-cell assay used: "RNA", "ATAC", or "Multiome"
#' @param filter.cells Optional: a character vector of cell barcodes to filter the read table by (default: NULL)
#' @param fastq.A Optional: provide path to FASTQ containing cell barcodes
#' @param fastq.B Optional: provide path to FASTQ containing sample tags
#' @param pos.cbc Optional: numeric vector specifying start & end indices of cell barcode sequence in FASTQ A
#' @param pos.sbc Optional: numeric vector specifying start & end indices of sample tag sequence in FASTQ B
#' @param pos.umi Optional: numeric vector specifying start & end indices of UMI sequence in FASTQ A or B (depends on barcode & assay)
#'
#' @return Data.frame containing variables "Cell", "UMI", and "Sample"
#'
#' @examples
#' ~/Experiment2
#' | Exp2MULTI_S3_L002_R1_001.fastq.gz
#' | Exp2MULTI_S3_L002_R2_001.fastq.gz
#'
#' read_table <- readTags(dir = "~/Experiment2",
#'                        name = "Exp2MULTI",
#'                        barcode.type = "MULTIseq",
#'                        assay = "RNA",
#'                        filter.cells = exp2_cells)
#'
#' @export
readTags <- function(dir,
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


###############
## alignTags ##
###############
#' Generate cell x tag count matrix from read table
#'
#' Tally total UMI counts of each sample barcode per cell in read table to produce a tag count matrix.
#' Updated version of MULTIseq.align() from original deMULTIplex package.
#'
#' @param read_table Data.frame containing variables "Cell", "UMI", and "Sample" as output by readTags() function
#' @param tag.ref Character vector of sample tag sequences used in experiment
#' @param filter.cells Optional: a character vector of cell barcodes to filter the read table by. Redundant if already filtered with readTags() (default: NULL)
#' @param string.dist.method Specify method for calculating string distance between reference tags and tag reads when correcting for sequencing errors. See ?stringdist for options (default: hamming)
#' @param max.dist Specify maximum string distance allowed for correcting sequencing errors in sample tags (default: 1)
#'
#' @return Sparse matrix of cell x sample tag counts
#'
#' @examples
#' > head(read_table, n = 3)
#'             Cell      UMI   Sample
#' CTTAGGCCTGAACATA TTTCACAC TGAGACCT
#' TCTGAGCCTAAACTGA CAAAGAGG CCACAATG
#' TTCTAGACTGAATTGA GATACGCA TGAGACCT
#'
#' tag_mtx <- alignTags(read_table,
#'                      tag.ref,
#'                      filter.cells = exp2_cells)
#'
#' @export
alignTags <- function(read_table,
                      tag.ref,
                      filter.cells = NULL,
                      string.dist.method = "hamming",
                      max.dist = 1) {
    t0 <- Sys.time()
    string.dist.method <- match.arg(string.dist.method)
    if (is.null(filter.cells)) {
        cells <- unique(read_table$Cell)
        cells <- cells[cells != paste(rep("G",16),collapse = "")]
    } else { cells <- filter.cells }

    if (max.dist > 0 ) {
        cat(paste0("Accounting for mismatch with maximum ", string.dist.method, " distance ", max.dist), fill = T)
        unique_tags <- unique(read_table$Sample)
        corrected_tags<-sapply(unique_tags, function(x){
            if(!x %in% tag.ref) {
                x2 = tag.ref[which(stringdist(x, tag.ref, method = string.dist.method)==max.dist)]
                if(length(x2) == 1) x2 else NA
            } else x
        })
        read_table$Sample <- corrected_tags[read_table$Sample]
    }

    cat("Assembling tag count table...", fill = T)
    dt <- data.table(read_table[,c("Cell", "Sample")])
    dt <- dt[complete.cases(dt), ]
    assign("dt", dt, env = .GlobalEnv)
    cnt_ind <- dt[, list(Freq =.N), by=list(Cell,Sample)]
    cnt_ind$i <- (1:length(cells))[match(cnt_ind$Cell, cells)] # NA?
    cnt_ind$j <- (1:length(tag.ref))[match(cnt_ind$Sample, tag.ref)]
    cnt_ind <- cnt_ind[complete.cases(cnt_ind),] # Added to account for provided cell list not in data
    tag_mtx <- sparseMatrix(i = cnt_ind$i, j = cnt_ind$j, x = cnt_ind$Freq, dims = c(length(cells), length(tag.ref)))
    rownames(tag_mtx) <- cells
    colnames(tag_mtx) <- tag.ref

    cat("Finished in",
        round(difftime(Sys.time(), t0, units = "secs")[[1]]),
        "seconds", sep = " ")

    cat("\n")
    return(tag_mtx)
}


#####################
### revComplement ###
#####################
#' Generate Reverse Complement Sequence
#'
#' Takes cell barcodes formatted by Cell Ranger and/or ArchR and returns reverse complement of each cell barcode sequence.
#' Depending on sequencer used, this may be needed to properly align cell barcodes.
#'
#' @param cells A character vector of cell barcodes (may include non-barcode formatting characters)
#' @param cbLength Length of cell barcode (default: 16)
#' @param keepFormat Whether or not to keep non-barcode formatting characters (e.g. "AGTCAGTCAGTCAGTC-1" or "Sample#AGTCAGTCAGTCAGTC")
#'
#' @return A character vector
#'
#' @export
revComplement <- function(cells,
                          cbLength = 16,
                          keepFormat = FALSE) {
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














