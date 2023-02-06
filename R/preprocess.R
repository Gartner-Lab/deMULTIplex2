



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
## revcompBC() ##
#################

# Description:
# Takes cell barcodes formatted by Cell Ranger and/or ArchR and returns reverse complement of each cell barcode sequence
#' @export
revcompBC <- function(cells,
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

## readMULTI()
# previously MULTIseq.preprocess()
# memory intensive for large FASTQ files

# Description:
# Locates and reads into memory FASTQ files produced by sequencing MULTIseq or MULTI-ATAC barcode libraries.
# Read pairs are parsed to extract 10x cell barcodes, MULTI* sample barcodes, and UMIs to count barcodes captured, and these are assembled into a read table
# Read table can be optionally filtered by a set of known cell barcodes (i.e. cells identified by Cell Ranger)
#' @export
readMULTI <- function(dir,
                      name,
                      barcode, # MULTIseq or MULTI-ATAC
                      assay, # RNA, ATAC, or Multiome
                      filterCells = NULL,
                      A = NA, B = NA,
                      pos.cbc, pos.sbc, pos.umi) {
    require(ShortRead)

    t0 <- Sys.time()

    if (barcode == "MULTIseq" & assay == "RNA") {
        rA <- "R1" # minimum 28bp read that contains cell BC and UMI
        rB <- "R2" # typically 91bp read that contains sample BC
    } else if (barcode == "MULTIseq" & assay == "Multiome") {
        rA <- "R1" # minimum 28bp read that contains cell BC and UMI
        rB <- "R3" # typically 91bp read that contains sample BC
    } else if (barcode == "MULTI-ATAC") {
        rA <- "R2" # minimum 16bp "index" read that contains cell BC
        rB <- "R3" # typically > 50bp read that contains sample BC and UMI
    } else {
        return(message("Error - barcode should be either MULTIseq or MULTI-ATAC"))
    }

    # set assay = NULL in order to use custom position indices (must manually define all 3)
    if (barcode == "MULTIseq" & assay %in% c("RNA","Multiome")) {
        pos.cbc <- c(1,16)
        pos.sbc <- c(1,8)
        pos.umi <- c(17,28)
    } else if (barcode == "MULTI-ATAC" & assay == "ATAC") {
        pos.cbc <- c(1,16)
        pos.sbc <- c(1,8)
        pos.umi <- c(9,16)
    } else if (barcode == "MULTI-ATAC" & assay == "Multiome") {
        pos.cbc <- c(9,24)
        pos.sbc <- c(1,8)
        pos.umi <- c(17,28)
    } else { cat("Using custom read positions") }

    if (file.exists(as.character(A))) {
        rA <- A
    } else {
        rA <- list.files(path = dir, pattern = paste(name, rA, ".fastq", sep=".*"), full.names = T)
    }

    if (file.exists(as.character(B))) {
        rB <- B
    } else {
        rB <- list.files(path = dir, pattern = paste(name, rB, ".fastq", sep=".*"), full.names = T)
    }

    if (length(rA) < 1 | length(rB) < 1) {
        return(message("Error - one or more FASTQ files not found"))
    } else if (length(rA) > 1 | length(rB) > 1) {
        return(message("Error - too many files with match names found"))
    }


    cat("### Building Read Table ###", fill = T)
    cat("Barcode Type: ", barcode, sep = "", fill = T)
    cat("10x Genomics Assay Type: ", assay, sep = "", fill = T)
    cat("Loading first set of reads...", fill = T)
    r <- readFastq(rA)
    gc(verbose = F)

    if (barcode == "MULTIseq") {
        cat("Extracting Cell Barcodes & UMIs...", fill = T)
        readTable <- data.frame(Cell = subseq(sread(r), pos.cbc[1], pos.cbc[2]),
                                UMI = subseq(sread(r), pos.umi[1], pos.umi[2]))
    }
    if (barcode == "MULTI-ATAC") {
        cat("Extracting Cell Barcodes...", fill = T)
        readTable <- data.frame(Cell = subseq(sread(r), pos.cbc[1], pos.cbc[2]))
    }
    cat("Finished processing first set of reads; unloading from memory", fill = T)
    r <- NULL
    gc(verbose = F)

    cat("Loading second set of reads...", fill = T)
    r <- readFastq(rB)
    gc(verbose = F)

    if (barcode == "MULTIseq") {
        cat("Extracting Sample Barcodes...", fill = T)
        readTable$Sample <- as.character(subseq(sread(r), pos.sbc[1], pos.sbc[2]))
    }
    if (barcode == "MULTI-ATAC") {
        cat("Extracting Sample Barcodes & UMIs...", fill = T)
        readTable$UMI <- as.character(subseq(sread(r), pos.umi[1], pos.umi[2]))
        readTable$Sample <- as.character(subseq(sread(r), pos.sbc[1], pos.sbc[2]))
    }
    cat("Finished processing second set of reads; unloading from memory", fill = T)
    r <- NULL
    gc(verbose = F)

    cat("Finished parsing ", nrow(readTable), " read pairs", sep = "", fill = T)

    if (is.null(filterCells)) {
        cat("Keeping all reads because filterCells = NULL", fill = T)
    } else if (!is.null(filterCells)) {
        cat("Filtering for ", length(filterCells), " provided cell barcodes...", sep ="", fill = T)
        ind <- which(readTable$Cell %in% filterCells)
        readTable <- readTable[ind, ]
    }

    readTable <- readTable[,c("Cell","Sample","UMI")]

    cat("Finished building read table", fill = T)
    cat("\n")
    cat("Finished in",
        round(difftime(Sys.time(), t0, units = "mins")[[1]],1),
        "minutes", sep = " ")
    cat("\n")
    cat("Preview of Read Table:", fill = T)
    cat("\n")
    print(head(readTable), quote = F)

    return(readTable)
}






##################
## alignMULTI() ##
##################

# previously MULTIseq.align()

# Description:
# Utilizes data.table library to quickly tally total UMI counts of each sample barcode per cell in read table
#' @export
alignMULTI <- function(readTable,
                       barcodes,
                       filterCells = NULL,
                       names = NULL) {
    require(data.table)
    require(Matrix)
    require(stringdist)

    t0 <- Sys.time()

    if (is.null(filterCells)) {
        cells <- unique(readTable$Cell)
        cells <- cells[cells != paste(rep("G",16),collapse = "")]
    } else { cells <- filterCells }

    if (is.null(names)) {
        names <- paste("Bar", 1:length(barcodes), sep = "")
    }

    cat("Deduplicating & Counting Sample Barcode UMIs...", fill = T)
    dt <- data.table(readTable)
    cnt <- dt[, list(Freq =.N), by=list(Cell,Sample,UMI)] # deduplicate UMIs (Freq not informative here)
    cnt2 <- cnt[, list(Freq =.N), by=list(Cell,Sample)] # tally up UMIs of each Sample per Cell
    cnt_ind <- cnt2

    cnt_ind$i <- (1:length(cells))[match(cnt2$Cell, cells)] # finds index of of each cell barcode in reference list
    cnt_ind$j <- (1:length(barcodes))[match(cnt2$Sample, barcodes)] # finds index of of each sample barcode in reference list

    # Hamming-distance sample barcode correction (for each non-match, see if there is a match with Hamming distance = 1)
    cat("Performing Hamming-Distance Sequencing Error Correction...", fill = T)
    idx <- which(is.na(cnt_ind$j))
    cnt_NA <- cnt_ind[idx,'Sample']
    cnt_NA <- apply(cnt_NA, 1, function(x) {
        tmp <- which(stringdist(x, barcodes, method = "hamming") == 1)
        if (length(tmp) < 1) {
            tmp <- NA
        }
        tmp
    })
    cnt_ind$j[idx] <- cnt_NA

    # Build sparse matrix
    cat("Assembling Barcode Count Table...", fill = T)
    cnt_ind <- cnt_ind[complete.cases(cnt_ind),] # remove rows that didn't match to cell and/or sample reference list (contain NAs)
    cnt_mtx <- sparseMatrix(i = cnt_ind$i,
                            j = cnt_ind$j,
                            x = cnt_ind$Freq,
                            dims = c(length(cells), length(barcodes)))
    colnames(cnt_mtx) <- names
    rownames(cnt_mtx) <- cells
    barTable <- data.frame(as.matrix(cnt_mtx))
    barTable$nUMI <- rowSums(barTable)

    cat("Finished in",
        round(difftime(Sys.time(), t0, units = "secs")[[1]]),
        "seconds", sep = " ")

    cat("\n")
    cat("Preview of Barcode Count Table:", fill = T)
    cat("\n")
    print(head(barTable), quote = F)
    return(barTable)
}

















