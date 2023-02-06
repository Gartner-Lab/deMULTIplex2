


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












