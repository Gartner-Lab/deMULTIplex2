

#' @export
confusion_stats <- function(call_label, true_label,
                            tag_mapping,
                            call.multiplet = 'Multiplet',
                            true.multiplet = "doublet",
                            plot.path =getwd(), plot.name = 'benchmark_', width = 3.5, height = 2.5) {
    call_label <- as.character(call_label[names(true_label)])
    call_label[is.na(call_label)] <- 'NA'
    acc_mtx <- as.data.frame.matrix(table(call_label, true_label))
    breaksList <- seq(0,1e3,by=1e2)
    graphics.off() # Shut down open plotting devices if any left over
    pdf(paste0(plot.path, plot.name, "_accmtx", ".pdf"), width = 3.5, height = 2.5)
    pheatmap(acc_mtx, cluster_rows = F, cluster_cols = F, display_numbers =T, number_format ="%.0f", fontsize = 8, fontsize_number=8, color = colorRampPalette(c("white", "#c73647"))(length(breaksList)), breaks = breaksList)
    dev.off()

    unq_tags <- tag_mapping$tag
    tag_stats <- list()
    for(tag in unq_tags) {
        tp <- sum(call_label == tag & true_label %in% tag_mapping$true_label[tag_mapping$tag == tag], na.rm=T)
        fp <- sum(call_label == tag, na.rm=T) - tp
        fn <- sum(call_label != tag & true_label %in% tag_mapping$true_label[tag_mapping$tag == tag], na.rm=T) # Note when multiple-to-one mapping exists this code does not work
        precision = tp/(tp+fp); if(is.na(precision)) precision = 0
        recall = tp/(tp+fn); if(is.na(recall)) recall = 0
        f_score <- tp / (tp + 0.5 * (fp + fn)); if(is.na(f_score)) f_score = 0
        tag_stats[[tag]] <- c(
            tp = tp, fp=fp, fn=fn,
            precision=precision,
            recall=recall,
            f_score=f_score
        )
    }
    assign("tag_stats",tag_stats, env=.GlobalEnv)
    singlet_stats <- c(
        precision = mean(sapply(tag_stats, function(x)x['precision'])), # na.rm?
        recall = mean(sapply(tag_stats, function(x)x['recall'])),
        f_score = mean(sapply(tag_stats, function(x)x['f_score']))
    )

    doublet_called_singlet_rate <- sum(true_label == true.multiplet & call_label %in% tag_mapping$tag, na.rm=T) / sum(true_label == true.multiplet, na.rm=T)
    doublet_called_negative_rate <- sum(true_label == true.multiplet & !call_label %in% c(tag_mapping$tag, call.multiplet), na.rm=T) / sum(true_label == true.multiplet, na.rm=T)
    doublet_called_doublet_rate <- sum(true_label == true.multiplet & call_label == call.multiplet, na.rm=T) / sum(true_label == true.multiplet, na.rm=T)
    doublet_stats <- c(recall = doublet_called_doublet_rate , doublet_called_singlet = doublet_called_singlet_rate, doublet_called_negative = doublet_called_negative_rate)

    return(list(
        acc_mtx = acc_mtx,
        tag_stats = tag_stats,
        singlet_avg_stats = singlet_stats,
        doublet_avg_stats = doublet_stats
    ))
}

#' @export
benchmark_demultiplex2 <- function(tag_mtx, true_label,
                                   tag_mapping,
                                   true.multiplet = "doublet",
                                   plot.path =getwd(), plot.name = 'benchmark_demultiplex2', width = 3.5, height = 2.5,
                                   seed = 1,
                                   init.cos.cut = .5,
                                   converge.threshold = 1e-3,
                                   prob.cut = 0.5,
                                   max.cell.fit = 1000,
                                   max.iter = 30,
                                   min.quantile.fit = 0.05, # Remove cells with total umi less than the specified quantile, which could be beads
                                   max.quantile.fit = 0.95, # Remove cells with total umi greater than the specified quantile, which could be multiplets
                                   residual.type = c("rqr", 'pearson'), # ONLY use RQR for future
                                   plot.umap = c("residual", "umi"),
                                   plot.diagnostics = T) {
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
                          min.quantile.fit = min.quantile.fit, # Remove cells with total umi less than the specified quantile, which could be beads
                          max.quantile.fit = max.quantile.fit, # Remove cells with total umi greater than the specified quantile, which could be multiplets
                          residual.type = residual.type, # ONLY use RQR for future
                          plot.umap = plot.umap,
                          plot.diagnostics = plot.diagnostics,
                          plot.path = plot.path,
                          plot.name = plot.name,
                          umap.nn = 30L,
                          seed = 1,
                          point.size = 1,
                          label.size = 3,
                          min.bc.show = 50)
    end_time <- Sys.time()
    calls <- res$assign_table$barcode_assign
    calls[res$assign_table$barcode_count == 0] = "Negative"
    calls[res$assign_table$barcode_count > 1] = "Multiplet"
    names(calls) <- rownames(res$assign_table)

    use_time = end_time - start_time
    # Benchmarking
    conf_stats <- confusion_stats(calls, true_label,
                             tag_mapping,
                             call.multiplet = 'Multiplet',
                             true.multiplet = true.multiplet,
                             plot.path = plot.path, plot.name = plot.name, width = width, height = height)
    return(
        c(
            list(res = res),
            conf_stats,
            time = use_time)
    )
}


#' @export
benchmark_demultiplex1 <- function(tag_mtx, true_label,
                                   tag_mapping,
                                   true.multiplet = "doublet",
                                   plot.path =getwd(), plot.name = 'benchmark_demultiplex1', width = 3.5, height = 2.5) {
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
                                 tag_mapping,
                                 call.multiplet = 'Doublet',
                                 true.multiplet = true.multiplet,
                                 plot.path = plot.path, plot.name =plot.name, width = width, height = height)

    return(
        c(list(res = final.calls),
          conf_stats,
          time = use_time)
    )
}


#' @export
benchmark_HTODemux <- function(tag_mtx, rna_mtx,
                               true_label,
                               tag_mapping,
                               true.multiplet = "doublet",
                               plot.path =getwd(), plot.name = 'benchmark_HTODemux', width = 3.5, height = 2.5) {
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
    seu_call <- as.character(seu_obj$hash.ID)
    names(seu_call) <- colnames(seu_obj)
    end_time <- Sys.time()
    use_time = end_time - start_time
    conf_stats = confusion_stats(seu_call, true_label,
                                 tag_mapping,
                                 call.multiplet = 'Doublet',
                                 true.multiplet = true.multiplet,
                                 plot.path = plot.path, plot.name = plot.name, width = width, height = height)
    return(
        c(list(res=seu_call), conf_stats, time = use_time)
    )
}


#' @export
benchmark_demuxmix_full <- function(tag_mtx, rna_mtx,
                                    true_label,
                                    tag_mapping,
                                    true.multiplet = "doublet",
                                    plot.path = getwd(), plot.name = 'benchmark_demuxmix_full_', width = 3.5, height = 2.5) {
    start_time <- Sys.time()
    dmm <- demuxmix(as.matrix(t(tag_mtx)), rna = colSums(rna_mtx > 0))
    classLabels <- dmmClassify(dmm)
    classLabels$assign <- classLabels$HTO
    classLabels$assign[classLabels$Type != 'singlet'] <- classLabels$Type[classLabels$Type != 'singlet']
    call_demuxmixfull <- classLabels$assign
    names(call_demuxmixfull) <- rownames(classLabels)
    end_time <- Sys.time()
    use_time = end_time - start_time
    conf_stats <- confusion_stats(call_demuxmixfull, true_label,
                                  tag_mapping = tag_mapping,
                                  call.multiplet = 'multiplet',
                                  true.multiplet = true.multiplet,
                                  plot.path = plot.path , plot.name = plot.name, width = width, height = height)
    return(
        c(list(call=call_demuxmixfull), conf_stats, time = use_time)
    )
}


#' @export
benchmark_demuxmix_naive <- function(tag_mtx,
                                     true_label,
                                     tag_mapping,
                                     true.multiplet = "doublet",
                                     plot.path = getwd(), plot.name = 'benchmark_demuxmix_naive', width = 3.5, height = 2.5) {
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
    conf_stats <- confusion_stats(call_demuxmixNaive, true_label,
                                  tag_mapping = tag_mapping,
                                  call.multiplet = 'multiplet',
                                  true.multiplet = true.multiplet,
                                  plot.path = plot.path , plot.name = plot.name, width = width, height = height)
    return(
        c(list(call=call_demuxmixNaive), conf_stats, time = use_time)
    )
}


