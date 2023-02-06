



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
