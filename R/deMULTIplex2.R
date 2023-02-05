### deMULTIplex2 ###

# Objectives:
# 1. Include functionality for MULTI-ATAC barcode sequencing
# 2. Improve demultiplexing algorithm
# 3. Improve plotting features
# 4. Include data-scaling options

'%ni%' <- Negate('%in%')

# first 16 barcodes
barcodes <- c("GGAGAAGA","CCACAATG","TGAGACCT","GCACACGC", #1-4
              "AGAGAGAG","TCACAGCA","GAAAAGGG","CGAGATTC", #5-8
              "GTAGCACT","CGACCAGC","TTAGCCAG","GGACCCCA", #9-12
              "CCAACCGG","TGACCGAT","GCAACGCC","CAATCGGT") #13-16

# full list of 96 barcodes
# barcodes <- readxl::read_excel("/Users/dannyconrad/Box/Data/deMULTIplex2/MULTI-seq_oligos_Nov2019.xlsx")$`Barcode Sequence`

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



###############
## barHist() ##
###############

# Description:
# Plot log-scaled UMI histograms for each unique barcode (should appear bimodal) and total UMI counts
#' @export
barHist <- function(barTable,
                    minUMI = 10,
                    plotnUMI = T,
                    select = NULL,
                    scale_y = F,
                    bins = 100,
                    colors = NULL) {
    require(tidyr)
    require(ggplot2)
    require(cowplot)

    # make sure there is a nUMI column
    if ("nUMI" %ni% colnames(barTable)) {
        barTable$nUMI <- rowSums(barTable)
    }

    # select specific barcodes if desired
    if (!is.null(select)) {
        barTable <- barTable[,c(colnames(barTable)[select], "nUMI")]
    }

    # remove outliers and low-UMI cells (or empty droplets) to improve visualization
    lo <- max(quantile(barTable$nUMI,0.001), minUMI)
    hi <- quantile(barTable$nUMI,0.999)
    ind <- barTable$nUMI < hi & barTable$nUMI > lo
    barTable <- barTable[ind,]

    if (plotnUMI == FALSE) {
        barTable$nUMI <- NULL
    }

    # reorganize dataframe to facilitate plotting
    x <- pivot_longer(barTable, cols = 1:ncol(barTable), names_to = "Barcode", values_to = "Count")
    x$Barcode <- factor(x$Barcode, levels = colnames(barTable))

    p <- ggplot(x, aes(x = Count, fill = Barcode)) +
        # geom_histogram(bins = bins, color = "white", alpha = 0.6) +
        geom_histogram(bins = bins, alpha = 0.6) +
        scale_x_log10() +
        theme_bw() +
        # ggtitle(paste(nrow(barTable), "cells with minimum", lo, "barcode UMIs", sep=" ")) +
        ggtitle(paste("nCell: ", nrow(barTable), "\n", "minUMI: ", lo, sep="")) +
        theme(plot.title = element_text(size=10))

    # element_text(family = NULL, face = NULL, colour = NULL, size = NULL,
    #              hjust = NULL, vjust = NULL, angle = NULL, lineheight = NULL,
    #              color = NULL)


    if (ncol(barTable) > 10) {
        ncol <- ceiling(sqrt(ncol(barTable)/3))
        p <- p + facet_wrap(facets = "Barcode", ncol = ncol, scales = "free_y")
    } else {
        p <- p + facet_grid(rows = "Barcode", scales = "free_y")
    }

    if (scale_y) {
        p <- p + scale_y_sqrt()
    }

    if (!is.null(colors)) {
        p <- p + scale_fill_manual(values = colors)
    }

    p
}


####################
## classifyLoop() ##
####################

# previously classifyCells() and classloop()

# Description:
# Modified classification algorithm that uses Seurat to normalize UMIs and classify cells based on cluster identity
# Produces plots for instant quality control: barcode UMI histograms showing the thresholds that were used for each barcode as well as UMAPs of barcode space colored by classification
#
# In in each round of classification (max 10):
# 1. perform quantile sweep on all barcodes to find the quantile-threshold that yields the highest proportion of singlets
# 2. classify cells using this quantile-threshold
# 3. remove "Negative"-classified cells
# 4. repeat until no more Negatives produced
#' @export
classifyLoop <- function(barTable,
                         counts = "raw", # raw, sct_count, sct_res
                         doubIdent = F,
                         col = NULL,
                         returnSeu = F) {
    require(Seurat)
    require(deMULTIplex)
    require(dplyr)
    require(cowplot)
    require(lsa)

    # exclude "nUMI" column
    if (tail(colnames(barTable), n = 1) == "nUMI") {
        barTable <- barTable[,1:(ncol(barTable)-1)]
    }

    # set color palette for plotting at the end
    if (is.null(col) & ncol(barTable) <= 8) {
        require(RColorBrewer)
        col <- brewer.pal(ncol(barTable), name = "Set1")[1:ncol(barTable)]
    } else if (is.null(col)) {
        require(viridis)
        col <- viridis(n = ncol(barTable))
    } else if (!is.null(col) & length(col) != ncol(barTable)) {
        return( message("Error: Number of colors must match number of barcodes when setting manual palette") )
    }

    cat("Applying Seurat functions:", fill = T)
    cat("  CreateSeuratObject", fill = T)
    seu <- CreateSeuratObject(data.frame(t(barTable)), verbose = F,)
    cat("  SCTransform", fill = T)
    seu <- SCTransform(seu, verbose = F)
    cat("  RunPCA", fill = T)
    seu <- RunPCA(seu, verbose = F, approx = F)
    cat("  RunUMAP", fill = T)
    seu <- RunUMAP(seu, dims = 1:ncol(seu@reductions$pca@cell.embeddings), verbose = F)
    # cat("  FindNeighbors", fill = T)
    # seu <- FindNeighbors(seu, dims = 1:ncol(seu@reductions$pca@cell.embeddings), verbose = F)
    # cat("  FindClusters", fill = T)
    # seu <- FindClusters(seu, resolution = 2, verbose = F)

    # cat("  FindClusters", fill = T)
    # seu <- FindClusters(seu, resolution = 0.5, verbose = F)
    # g <- list(DimPlot(seu, label = T) + NoLegend(),
    #           FeaturePlot(seu, features = 'nCount_RNA', min.cutoff = 'q05', max.cutoff = 'q95'),
    #           VlnPlot(seu, features = 'nCount_RNA', log = T, pt.size = 0,) + NoLegend())
    # plot_grid(plot_grid(plotlist = g[1:2]), g[[3]], ncol = 1)


    # Use raw counts or transform/normalize?
    if (counts == "raw") {
        barCounts <- barTable # Raw barcode counts
    } else if (counts == "sct_count") {
        barCounts <- seu@assays$SCT@counts %>% as.data.frame() %>% t() %>% as.data.frame() # SCTransform-normalized barcode counts
    } else if (counts == "sct_res") {
        barCounts <- seu@assays$SCT@data %>% as.data.frame() %>% t() %>% as.data.frame() # SCTransform residuals
    } else { message(paste("Error: \"", counts, "\"", " is not a valid datatype", sep = "")) }



    # Calculate cosine similarity against canonical vectors
    cat("Calculating Cosine Similarities...", fill = T)
    mtx_canon <- diag(x = 1, nrow = ncol(barCounts), ncol = ncol(barCounts))
    # mtx_canon <- matrix(0L, nrow = ncol(barCounts), ncol = ncol(barCounts))
    # diag(mtx_canon) <- 1
    cosSim <- t(apply(barCounts, 1, function(x) {
        apply(mtx_canon, 1, function(i) {
            cosine(x, i) })})) %>% as.data.frame()
    colnames(cosSim) <- colnames(barCounts)

    # Start classification loop (max 10 rounds)
    cat("Starting Classification Loop...", fill = T)
    counter <- 0
    neg.cells <- c()

    while (counter < 10) {
        counter <- counter + 1
        cat("Round ", counter, "...", sep="", fill = T)

        # Perform quantile sweep
        quantiles <- seq(0.01, 0.99, by=0.02)
        quantile_sweep <- lapply(quantiles, function(q) {
            classify(barCounts, cosSim, q, neg.cells)[[1]]
        })
        names(quantile_sweep) <- paste("q=",quantiles,sep="")

        # identify optimal threshold from sweep
        res <- findThresh(call.list=quantile_sweep)
        qOpt <- findQ(res$res, res$extrema)

        classification <- classify(barCounts, cosSim, qOpt, neg.cells, doubIdent)
        calls <- classification[[1]]
        thresh <- classification[[2]]

        neg.tmp <- names(calls)[which(calls == "Negative")]
        neg.cells <- c(neg.cells, neg.tmp)

        if (counter == 1) {

            # plot barcode histograms & thresholds determined in first round
            df_bc <- data.frame(pivot_longer(barCounts, cols = 1:ncol(barCounts), names_to = "Barcode", values_to = 'Counts'))
            df_bc$Barcode <- factor(df_bc$Barcode, levels = colnames(barCounts))
            df_thresh <- data.frame(Thresh = do.call(c, thresh),
                                    Barcode = factor(colnames(barCounts), levels = colnames(barCounts)))

            p1 <- barHist(barCounts, plotnUMI = F, scale_y = T, colors = col, minUMI = 0, bins = 50) +
                theme_classic() +
                geom_vline(data = df_thresh, mapping = aes(xintercept = Thresh), size = 1, color = 'grey')
        }

        # print classification results
        cat(paste("\t", "Singlets: ", length(which(calls %in% colnames(barTable))), sep=""),
            paste("\t", "Doublets: ", length(which(calls %ni% c("Negative",colnames(barTable)))), sep=""),
            paste("\t", "Negatives: ", length(which(calls == "Negative")), sep=""),
            sep = "\n")

        # if no more negatives, classification is done!
        if (length(neg.tmp) == 0) { break }

    }
    calls_final <- c(calls, rep("Negative", length(neg.cells)))
    names(calls_final) <- c(names(calls), neg.cells)
    print(table(calls_final))

    seu$Calls <- calls_final[colnames(seu)]

    callTypes <- unique(calls_final)
    doubletTypes <- callTypes[callTypes %ni% c(colnames(barTable), "Negative")]

    seu$Calls <- factor(seu$Calls,
                        levels = c(colnames(barCounts), "Negative", doubletTypes),
                        labels = c(colnames(barCounts), "Negative", rep("Doublet", length(doubletTypes))))

    p2 <- DimPlot(seu, group.by = "Calls", label = T, pt.size = 1, cols = c(col, "#A9A9A9", "#000000"))
    p3 <- FeaturePlot(seu, features = "nCount_RNA", pt.size = 1, max.cutoff = 'q99', min.cutoff = 'q01') + ggtitle("nUMI")
    # p3 <- PCAPlot(seu, group.by = "Calls", pt.size = 1)


    p <- plot_grid(
        p1, p2, p3,
        align = 'hv',
        nrow = 1,
        rel_widths = c(2,3,3)
    )

    print(p)

    if (returnSeu) {
        return(seu)
    } else {
        return(calls_final)
    }
}



################
## classify() ##
################

# Ignore: Function used internally by classifyLoop()

# Description:
# Takes barcode count table and a threshold q-value and classifies cells
#' @export
classify <- function(barCounts,
                     cosSim,
                     q,
                     negCells,
                     doubIdent = F) {
    require(tidyr)

    # Remove negative cells
    barCounts <- barCounts[!rownames(barCounts) %in% negCells,]
    cosSim <- cosSim[!rownames(cosSim) %in% negCells,]

    n_BC <- ncol(barCounts)
    n_cells <- nrow(barCounts)

    # Use upper and lower quantiles of cossine similarity to determine high/low barcode count maxima
    # Calculate count threshold for each barcode from the given quantile between the two maxima
    thresh <- lapply(1:n_BC, function(x) {
        df <- data.frame(Count = barCounts[,x],
                         CosSim = cosSim[,x])
        # high <- median(df$Count[df$CosSim >= quantile(df$CosSim,0.99)])
        # low <- median(df$Count[df$CosSim <= quantile(df$CosSim,0.20)])
        high <- median(df$Count[df$CosSim >= quantile(df$CosSim,0.95)])
        low <- median(df$Count[df$CosSim <= quantile(df$CosSim,0.05)])
        # high <- median(df$Count[df$CosSim > 0.95])
        # low <- median(df$Count[df$CosSim < 0.05])
        quantile(c(high, low),q)
    })

    # Assign classifications based on thresholds
    calls_tmp <- rep("Negative", n_cells)
    names(calls_tmp) <- rownames(barCounts)


    for (bc in 1:n_BC) {
        # Which cells are positive for each barcode
        pos <- which(barCounts[,bc] >= thresh[[bc]])
        if (length(pos) == 0) { return() }

        # Update the calls_tmp vector with this info
        calls_tmp[pos] <- sapply(calls_tmp[pos], function(x) {
            # if a cell was deemed "positive" for a previous barcode, this is now a doublet
            if (x == "Negative") {
                call <- colnames(barCounts)[bc]
            } else if (doubIdent == T) {
                call <- paste(x, colnames(barCounts)[bc], sep="/") # use doubIdent option to independently record each unique doublet/multiplet type (i.e. Bar1/Bar3)
            } else {
                call <- "Doublet"
            }
            return(call)
        })
    }

    return(list(calls_tmp, thresh))

}

##################
## barHeatmap() ##
##################

# Description:
# After classifications are complete, plot barcode barcode x classification heatmap
#' @export
barHeatmap <- function(barTable,
                       calls,
                       log = T,
                       colLow = "white",
                       colHigh = "#c73647") {
    require(pheatmap)
    require(tidyr)

    if (sum(names(calls) %ni% rownames(barTable)) > 0) {
        return(message("Error: cell barcodes do not match"))
    }

    if (log) {
        barTable <- log10(barTable)
    }

    barTable[barTable == -Inf] <- 0

    calltypes <- unique(calls)

    tmp <- sapply(calltypes, function(x) {
        df <- barTable[names(calls[calls == x]),]
        colMeans(df)
    })
    colnames(tmp) <- calltypes
    rownames(tmp) <- colnames(barTable)

    bar <- colnames(barTable)[colnames(barTable) != "nUMI"] %>% as.character()
    doub <- calltypes[calltypes %ni% c(bar, "Negative", NA)] %>% as.character()
    neg <- calltypes[calltypes %ni% c(bar, doub)] %>% as.character()

    order <- c(bar, doub, neg)

    tmp <- tmp[,order]

    tmp_df <- as.data.frame(tmp)
    tmp_df$Barcode <- rownames(tmp_df)
    tmp_df <- pivot_longer(tmp_df, cols = 1:ncol(tmp), names_to = "Call", values_to = "Mean")
    tmp_df$Call <- factor(tmp_df$Call, levels = order)
    tmp_df$Barcode <- factor(tmp_df$Barcode, levels = rev(rownames(tmp)))

    ggplot(tmp_df, aes(x = Call, y = Barcode, fill = Mean, label = round(Mean,1))) +
        geom_tile(color = 'grey') +
        geom_text() +
        theme_classic() +
        scale_fill_gradient(low = colLow, high = colHigh) +
        theme(axis.text.x = element_text(angle = -45, hjust = 0))
}


#####################
## classifyMULTI() ##
#####################

# Description:
# The new demultiplex function with "Correct" Pearson residual and cosine similarity
#' @importFrom Matrix Matrix rowSums
#' @importFrom MASS glm.nb
#' @importFrom ggrastr geom_point_rast
#' @importFrom gridExtra arrangeGrob
#' @importFrom magrittr %>%
#' @export
classifyMULTI <- function(barTable,
                          posThresh = 0.2,
                          cosineThresh = seq(0.5,0.9,0.1),
                          plotUMAP = "cumi.cos",
                          plotDiagnostics = T,
                          plotPath = getwd(),
                          UMAPNeighbors = 30L,
                          seed = 1) {
    set.seed(seed)
    # TODO: Need to handle cells if rowSums = 0
    bc_mtx <- barTable

    if (any(rowSums(bc_mtx)) == 0) {
        stop("Please remove any cells with 0 total barcode counts")
    }

    # Generate canonical vector matrix for calculating cosine similarities
    mtx_canon <- diag(x = 1, nrow = ncol(bc_mtx), ncol = ncol(bc_mtx))
    colnames(mtx_canon) <- colnames(bc_mtx)
    rownames(mtx_canon) <- colnames(bc_mtx)

    # Calculate cosine similarity matrix
    bc_mtx <- Matrix(as.matrix(bc_mtx), sparse = T)
    vec_norm2 <- apply(bc_mtx, 1, function(x) norm(x, type = "2")) # for each cell, calculate its L2 norm (Euclidean distance from origin)
    cos_mtx_raw <- bc_mtx / vec_norm2 # divide matrix by L2 norms to get cosine distances

    # Matrix for storing Pearson residual (from NB fit)
    pr_mtx <- matrix(nrow = nrow(bc_mtx), ncol = ncol(bc_mtx), data = NA)
    rownames(pr_mtx) <- rownames(bc_mtx)
    colnames(pr_mtx) <- colnames(bc_mtx)

    # Matrix for storing corrected UMI counts (inverse NB fit, with library size set to median)
    cumi_mtx <- pr_mtx

    # NB fit
    barcodes <- colnames(bc_mtx)
    model_list <- list()

    for(bc in barcodes) {
        cat("Fitting GLM-NB for ", bc, sep = "", fill = T)
        df <- data.frame(RawUMI = bc_mtx[, bc],
                         CosSim = cos_mtx_raw[,bc])
        df$nUMI <- rowSums(bc_mtx)
        # subset "negative" cells, i.e. below cosine similarity threshold
        neg_cells <- which(df$CosSim < max(df$CosSim, na.rm = T) * (1-posThresh))
        model <- glm.nb(RawUMI ~ log(nUMI), data = df[neg_cells,], link = log)
        model_list[[bc]] <- model

        predicted_res <- predict(model, df, type = "response", se.fit=TRUE)
        df_pred <- cbind(df, predicted_res)
        # Compute Pearson residual using background count: (Count - fit) / sqrt(fit  + fit^2/theta)
        df_pred$pearson_residual <- (df_pred$Count-df_pred$fit)/sqrt(df_pred$fit + df_pred$fit^2/model$theta)

        # Optional: return corrected UMI counts (not necessary for demultiplexing)
        med_scale_df <- data.frame(nUMI = median(df_pred$nUMI))
        expected_mu <- predict(model, med_scale_df, type = "response", se.fit=TRUE)
        mu <- expected_mu$fit
        theta <- model$theta
        variance <- mu + mu^2 / theta
        df_pred$corrected_count <- mu + df_pred$pearson_residual * sqrt(variance)
        df_pred$corrected_umi <- round(df_pred$corrected_count, 0)
        df_pred$corrected_umi[df_pred$corrected_umi < 0] <- 0

        cumi_mtx[,bc] <- df_pred$corrected_umi
        pr_mtx[,bc] <- df_pred$pearson_residual
    }


    # Recompute cosine metric
    pr_mtx <- Matrix(as.matrix(pr_mtx), sparse = T)
    vec_norm2 <- apply(pr_mtx, 1, function(x) norm(x, type = "2"))
    cos_mtx_pr <- pr_mtx / vec_norm2


    # Assign singlets
    assign_tbl <- matrix(nrow = nrow(bc_mtx), ncol = length(cosineThresh), data = NA)
    colnames(assign_tbl) <- paste0("assign.res_", cosineThresh)
    rownames(assign_tbl) <- rownames(bc_mtx)
    for(i in 1:length(cosineThresh)) {
        assign_mtx <- cos_mtx_pr > cosineThresh[i]
        assign_res <- apply(assign_mtx, 1, function(x) {
            if (length(which(x)) == 1) {
                barcodes[which(x)]
            } else {
                NA
            }
        })
        assign_tbl[,i] <- assign_res
    }

    # UMAP Plotting Options
    umap_res <-  NA
    if (plotUMAP == "UMI") { # Raw UMI Counts
        umap_res <- compute_umap(bc_mtx, use_dim = ncol(bc_mtx), n_component=2, n_neighbors = UMAPNeighbors)
    } else if (plotUMAP == "PR") { # Pearson Residuals
        umap_res <- compute_umap(pr_mtx, use_dim = ncol(pr_mtx), n_component=2, n_neighbors = UMAPNeighbors)
    } else if (plotUMAP == "PR_Cos") { # Cosine Similarity calculated from Pearson Residuals
        umap_res <- compute_umap(cos_mtx_pr, use_dim = ncol(cos_mtx_pr), n_component=2, n_neighbors = UMAPNeighbors)
    } else if (plotUMAP == "cUMI") { # Corrected UMI Counts
        umap_res <- compute_umap(cumi_mtx, use_dim = ncol(cumi_mtx), n_component=2, n_neighbors = UMAPNeighbors)
    } else if (plotUMAP == "cUMI_Cos") { # Cosine Similarity calculated from Corrected UMI Counts
        dot_res <- cumi_mtx %*% mtx_canon
        vec_norm2 <- apply(cumi_mtx, 1, function(x) norm(x, type = "2"))
        cos_mtx_cumi <- dot_res / vec_norm2
        umap_res <- compute_umap(cos_mtx_cumi,use_dim = ncol(cos_mtx_cumi), n_component=2, n_neighbors = UMAPNeighbors)
    } else if (plotUMAP == "None") {
    } else { stop("plotUMAP = ", plotUMAP, " is not valid") }

    # Diagnostic Plots
    if (plotDiagnostics) {
        glist <- list()
        for (bc in barcodes) {
            df <- data.frame(
                RawUMI = bc_mtx[, bc],
                CorUMI = cumi_mtx[, bc], # For now always set to TRUE so is computed
                PearsonRes = pr_mtx[,bc],
                CosSim_RawUMI = cos_mtx_raw[,bc],
                CosSim_PearsonRes = cos_mtx_pr[,bc])
            df$nUMI_Raw <- rowSums(bc_mtx)
            df$nUMI_Corrected <- rowSums(cumi_mtx)
            df <- cbind(df, df_pred[,c("fit", "se.fit", "UL", "LL")])
            mappings <- list(
                c("log(nUMI_Raw)", "log(RawUMI)", "CosSim_RawUMI"),
                c("log(nUMI_Raw)", "PearsonRes", "CosSim_RawUMI"),
                c("log(RawUMI)", "CosSim_RawUMI", "CosSim_RawUMI"),
                c("nUMI_Corrected", "CorUMI", "CosSim_PearsonRes"),
                c("log(nUMI_Corrected)", "PearsonRes", "CosSim_PearsonRes"),
                c("log(CorUMI)", "CosSim_PearsonRes", "CosSim_PearsonRes"),
                c("log(PearsonRes)", "CosSim_PearsonRes", "CosSim_PearsonRes")
            )

            if (plotUMAP != "None") {
                df <- cbind(df, umap_res)
                mappings <- append(mappings, list(c("UMAP_1", "UMAP_2", "CosSim_PearsonRes")))
            }

            plot_list <- list()
            for(i in 1:length(mappings)){
                map <- mappings[[i]]
                p <- ggplot(df, aes_string(map[1], map[2])) +
                    geom_point_rast(aes_string(color = map[3]), stroke = 0, size = 1) +
                    scale_color_gradientn(colors = get_numeric_color("BlueGreenRed")) +
                    ggtitle(bc) +
                    theme_bw() +
                    labs(color = gsub("_","\n",map[3]))
                plot_list[[i]] <- p
            }
            # plot_list[[length(mappings) + 1]] <- arrangeGrob(
            #     grobs = list(ggplot(df, aes(x = CosSim_RawUMI)) + geom_histogram(bins = 50, color = "white") + theme_bw() + theme(plot.margin = unit(c(0,1,0,0), "cm")),
            #                  ggplot(df, aes(x = CosSim_PearsonRes)) + geom_histogram(bins = 50, color = "white") + theme_bw() + theme(plot.margin = unit(c(0,1,0,0), "cm"))),
            #     nrow = 2)

            glist[[bc]] <- arrangeGrob(grobs = plot_list, ncol = 3)
        }

        time <- (Sys.time() %>% make.names() %>% strsplit(split = "X") %>% unlist())[2]
        pdf(paste0(plotPath, "/", time, "_diagnostics.pdf"), width = 12, height = 10)
        for(i in 1:length(glist)) {
            grid.newpage()
            grid.draw(glist[[i]])
        }
        dev.off()
    }

    return(
        list(
            res = assign_tbl %>% as.matrix(),
            umap = umap_res %>% as.matrix(),
            pr_mtx = pr_mtx %>% as.matrix(),
            cumi_mtx = cumi_mtx %>% as.matrix(),
            cos_mtx_raw = cos_mtx_raw %>% as.matrix(),
            cos_mtx_pr = cos_mtx_pr %>% as.matrix(),
            models = model_list
        )
    )
}



