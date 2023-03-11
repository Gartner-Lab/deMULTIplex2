

#' @export
confusion_stats <- function(call_label, true_label,
                            call.multiplet = 'NA', call.negative = 'NA',
                            true.multiplet = "doublet", true.negative = "unknown",
                            plot.path =getwd(), plot.name = 'benchmark_', width = 3.5, height = 2.5) {
    call_label <- as.character(call_label[names(true_label)])
    #call_label = sapply(strsplit(call_label, "-|[.]"), function(x)x[1]) # This is to handel multi map to 1 case, name needs to be specified well
    call_label[is.na(call_label) | call_label == 'NA'] = call.negative[1] # Treat NA as negative for now
    acc_mtx <- as.data.frame.matrix(table(call_label, true_label))

    breaksList <- seq(0,1e3,by=1e2)
    pdf(paste0(plot.path, plot.name, "_accmtx", ".pdf"), width = 3.5, height = 2.5)
    pheatmap(acc_mtx, cluster_rows = F, cluster_cols = F, display_numbers =T, number_format ="%.0f", fontsize = 8, fontsize_number=8, color = colorRampPalette(c("white", "#c73647"))(length(breaksList)), breaks = breaksList)
    dev.off()

    use_row <- rownames(acc_mtx)[!rownames(acc_mtx) %in% c(call.multiplet, call.negative)]
    use_col <- colnames(acc_mtx)[!colnames(acc_mtx) %in% c(true.multiplet, true.negative)]
    sub_mtx <- as.matrix(acc_mtx[use_row,use_col])
    #pc <- sum(diag(sub_mtx )) / sum(acc_mtx[use_row,use_col])
    precision.singlet = sum(diag(sub_mtx )) / sum(acc_mtx[use_row,])
    recall.singlet = sum(diag(sub_mtx )) / sum(acc_mtx[,use_col])
    precision.multiplet = sum(acc_mtx[call.multiplet, true.multiplet]) / sum(acc_mtx[call.multiplet,])
    recall.multiplet <- sum(acc_mtx[call.multiplet, true.multiplet]) / sum(acc_mtx[,true.multiplet])

    return(list(
        acc_mtx = acc_mtx,
        stats = c(precision.singlet = precision.singlet,
                  recall.singlet = recall.singlet,
                  precision.multiplet = precision.multiplet,
                  recall.multiplet = recall.multiplet)
    ))
}

#' @export
benchmark_demultiplex2 <- function(tag_mtx, true_label, plot.path =getwd(), plot.name = 'benchmark_demultiplex2', width = 3.5, height = 2.5,
                                   seed = 1,
                                   init.cos.cut = .5,
                                   converge.threshold = 1e-3,
                                   prob.cut = 0.5,
                                   max.cell.fit = 1000,
                                   max.iter = 30) {
    set.seed(seed)
    require(deMULTIplex2)
    start_time <- Sys.time()
    res <- demutiplexTags(tag_mtx,
                          init.cos.cut = init.cos.cut,
                          converge.threshold = converge.threshold,
                          max.iter = max.iter,
                          prob.cut = prob.cut,
                          min.cell.fit = 10,
                          max.cell.fit = max.cell.fit,
                          min.quantile.fit = 0.05, # Remove cells with total umi less than the specified quantile, which could be beads
                          max.quantile.fit = 0.95, # Remove cells with total umi greater than the specified quantile, which could be multiplets
                          residual.type = c("rqr", 'pearson'), # ONLY use RQR for future
                          plot.umap = c("residual", "umi"),
                          plot.diagnostics = T,
                          plot.path = plotPath,
                          plot.name = paste0(test_text),
                          umap.nn = 30L,
                          seed = 1,
                          point.size = 1,
                          label.size = 3,
                          min.bc.show = 50)
    assign("res2", res, env =.GlobalEnv)
    end_time <- Sys.time()
    calls <- res$assign_table$barcode_assign
    calls[res$assign_table$barcode_count == 0] = "Negative"
    calls[res$assign_table$barcode_count > 1] = "Multiplet"
    names(calls) <- rownames(res$assign_table)

    use_time = end_time - start_time
    # Benchmarking
    conf_stats = confusion_stats(calls, true_label,
                                 call.multiplet = 'Multiplet', call.negative = 'Negative',
                                 true.multiplet = "doublet", true.negative = "unknown",
                                 plot.path = plot.path, plot.name =plot.name, width = width, height = height)

    return(
        c(
            list(call = calls),
            conf_stats,
            time = use_time)
    )
}


#' @export
benchmark_demultiplex1 <- function(tag_mtx, true_label, plot.path =getwd(), plot.name = 'benchmark_demultiplex1', width = 3.5, height = 2.5) {
    require(deMULTIplex)

    start_time <- Sys.time()
    ## Round 1 -----------------------------------------------------------------------------------------------------
    ## Perform Quantile Sweep
    bar.table <- tag_mtx
    bar.table_sweep.list <- list()
    n <- 0
    for (q in seq(0.01, 0.99, by=0.02)) {
        print(q)
        n <- n + 1
        bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
        names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
    }
    ## Identify ideal inter-maxima quantile to set barcode-specific thresholds
    threshold.results1 <- findThresh(call.list=bar.table_sweep.list)
    ggplot(data=threshold.results1$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "none") +
        geom_vline(xintercept=threshold.results1$extrema, lty=2) + scale_color_manual(values=c("red","black","blue"))

    round1.calls <- classifyCells(bar.table, q=findQ(threshold.results1$res, threshold.results1$extrema))
    neg.cells <- names(round1.calls)[which(round1.calls == "Negative")]
    bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

    ## Round 2 -----------------------------------------------------------------------------------------------------
    bar.table_sweep.list <- list()
    n <- 0
    for (q in seq(0.01, 0.99, by=0.02)) {
        print(q)
        n <- n + 1
        bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
        names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
    }

    threshold.results2 <- findThresh(call.list=bar.table_sweep.list)
    round2.calls <- classifyCells(bar.table, q=findQ(threshold.results2$res, threshold.results2$extrema))
    neg.cells <- c(neg.cells, names(round2.calls)[which(round2.calls == "Negative")])

    ## Repeat until all no negative cells remain (usually 3 rounds)...
    final.calls <- c(round2.calls, rep("Negative",length(neg.cells)))
    names(final.calls) <- c(names(round2.calls),neg.cells)
    end_time <- Sys.time()

    use_time = end_time - start_time
    # Benchmarking
    conf_stats = confusion_stats(final.calls, true_label,
                                 call.multiplet = 'Doublet', call.negative = 'Negative',
                                 true.multiplet = "doublet", true.negative = "unknown",
                                 plot.path = plot.path, plot.name =plot.name, width = width, height = height)

    return(
        c(list(call = final.calls),
          conf_stats,
          time = use_time)
    )
}


#' @export
benchmark_HTODemux <- function(tag_mtx, rna_mtx, true_label, plot.path =getwd(), plot.name = 'benchmark_HTODemux', width = 3.5, height = 2.5) {
    require(Seurat)
    start_time <- Sys.time()
    seu_obj <- CreateSeuratObject(counts = rna_mtx)
    seu_obj[["HTO"]] <- CreateAssayObject(counts = t(tag_mtx))
    # Normalize HTO data, here we use centered log-ratio (CLR) transformation
    seu_obj <- NormalizeData(seu_obj, assay = "HTO", normalization.method = "CLR")
    seu_obj <- HTODemux(seu_obj, assay = "HTO", positive.quantile = 0.99)

    use_levels = sort(unique(as.character(seu_obj$hash.ID)))
    use_levels <- c(use_levels[!use_levels %in% c("Doublet",  "Negative")], c("Doublet",  "Negative"))
    seu_obj$hash.ID <- factor(as.character(seu_obj$hash.ID), levels = use_levels)
    assign("test1", seu_obj, env =.GlobalEnv)
    seu_call <- seu_obj$hash.ID
    end_time <- Sys.time()
    use_time = end_time - start_time
    conf_stats = confusion_stats(seu_call, true_label,
                                 call.multiplet = 'Doublet', call.negative = 'Negative',
                                 true.multiplet = "doublet", true.negative = "unknown",
                                 plot.path = plot.path, plot.name = plot.name, width = width, height = height)
    return(
        c(list(call=seu_call), conf_stats, time = use_time)
    )
}


#' @export
benchmark_demuxmix_full <- function(tag_mtx, rna_mtx, true_label, plot.path = getwd(), plot.name = 'benchmark_demuxmix_full_', width = 3.5, height = 2.5) {
    start_time <- Sys.time()
    dmm <- demuxmix(as.matrix(t(tag_mtx)), rna = colSums(rna_mtx > 0))
    classLabels <- dmmClassify(dmm)
    classLabels$assign <- classLabels$HTO
    classLabels$assign[classLabels$Type != 'singlet'] <- classLabels$Type[classLabels$Type != 'singlet']
    call_demuxmixfull <- classLabels$assign
    names(call_demuxmixfull) <- rownames(classLabels)
    end_time <- Sys.time()
    use_time = end_time - start_time
    conf_stats = confusion_stats(call_demuxmixfull, true_label,
                                 call.multiplet = 'multiplet', call.negative = c('negative', 'uncertain'),
                                 true.multiplet = "doublet", true.negative = "unknown",
                                 plot.path = plot.path , plot.name = plot.name, width = width, height = height)
    return(
        c(list(call=call_demuxmixfull), conf_stats, time = use_time)
    )
}


#' @export
benchmark_demuxmix_naive <- function(tag_mtx, true_label, plot.path = getwd(), plot.name = 'benchmark_demuxmix_naive', width = 3.5, height = 2.5) {
    start_time <- Sys.time()
    dmmNaive <- demuxmix(as.matrix(t(tag_mtx)), model = "naive")
    classLabelsNaive <- dmmClassify(dmmNaive)
    classLabelsNaive$assign <- classLabelsNaive$HTO
    classLabelsNaive$assign[classLabelsNaive$Type != 'singlet'] <- classLabelsNaive$Type[classLabelsNaive$Type != 'singlet']
    call_demuxmixNaive <- classLabelsNaive$assign
    names(call_demuxmixNaive) <- rownames(classLabelsNaive)
    test_text <- paste0("benchmark_Gaublomme_human_st_", "demuxmixNaive_")
    end_time <- Sys.time()
    use_time = end_time - start_time
    conf_stats = confusion_stats(call_demuxmixNaive, true_label,
                                 call.multiplet = 'multiplet', call.negative = c('negative', 'uncertain'),
                                 true.multiplet = "doublet", true.negative = "unknown",
                                 plot.path = plot.path , plot.name = plot.name, width = width, height = height)
    return(
        c(list(call=call_demuxmixNaive), conf_stats, time = use_time)
    )
}


