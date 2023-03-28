

#' Demultiplex tag count matrix
#'
#' Demultiplex tag count matrix to assign each cell the most probable sample identity using the deMULTIplex2 algorithm.
#' deMULTIplex2 is built on a statistical model of tag read counts derived from the physical mechanism of tag cross-contamination.
#' It uses generalized linear models and expectation-maximization to probabilistically infer the sample identity of each cell.
#'
#' @param tag_mtx (required) The tag UMI count matrix with each row representing a cell and each column representing a tag. Pre-filtering to remove empty droplets are recommended.
#' @param init.cos.cut Cosine cutoff to initialize the EM algorithm. If input is a single value, same cosine cutoff will be applied to all tags (recommended).
#' If input is a vector, must be the same length as the number of tags. EM will be initialized with the specified cosine cutoff for each tag.
#' @param max.iter Maximum number of EM iterations. For most datasets the default value should suffice.
#' @param converge.threshold If consecutive absolute difference in the expected value of the log likelihood (Q function) is below this value, the algorithm is fully converged.
#' In practice, parameter estimation will stabilize in a few iterations even with down-sampled cells, so Q difference not reaching this threshold is generally fine.
#' @param prob.cut Cutoff of posterior probability to determine if a cell is positive or negative with each tag. Default value (0.5) is recommended.
#' @param min.cell.fit If a tag has too few positive cells, a GLM model cannot be robustly fit.
#' Therefore, tags with positive cells less than this specified amount will not be processed.
#' @param max.cell.fit Robust estimation of GLM generally does not require using all cells. For large datasets,
#' deMULTIplex2 downsample the positive cells and negative cells to the amount specified by this parameter to speed up processing.
#' @param min.quantile.fit If a cell has total tag count below this quantile, it will be excluded from fitting for robustness.
#' @param max.quantile.fit If a cell has total tag count above this quantile, it will be excluded from fitting for robustness.
#' @param residual.type Report which type of residual and use which residual for computing UMAP. "`rqr`": Randomized Quantile Residual. "`pearson`": Pearson residual.
#' @param plot.umap Compute and plot UMAP with residual ("`residual`"), raw umi count ("`umi`") or do not output UMAP summary plots ("`none`").
#' @param plot.diagnostics Whether to plot diagnostic plots.
#' @param plot.path Path to output the plots.
#' @param plot.name Text string to append to the output file name.
#' @param umap.nn corresponds to the `n_neighbors` parameter in the `umap` function in the `uwot` package.
#' @param seed Random seed for reproducibility.
#' @param point.size Size of point in the summary and diagnostic plots.
#' @param label.size Label size in the summary plot.
#' @param min.tag.show Tags with positive cells below this threshold will not be labeled on the summary plot.
#' @return A list containing the assignment results.
#' `final_assign` is the final assignment in a single character vector.
#' `assign_table` stores the detailed assignment result, including the assigned tag for singlets, the total positive tag count of each cell,
#' and whether a cell is a singlet, a multiplet or is not tagged.
#' `umap` is the umap coordinates.
#' `res_mtx` is the residual matrix.
#' `prob_mtx` is the posterior probability of each cell being positively tagged by each tag.
#' `coefs` is the parameter estimates of the GLM models.
#' `df_list` is a list of data frames for each tag containing stats from the final iteration of EM.
#' @references To be added.
#'
#' @examples
#' res <- demultiplexTags(tag_mtx,
#'                        plot.path = "~/",
#'                        plot.name = "test",
#'                        plot.diagnostics = FALSE)
#' table(res$final_assign)
#'
#' @importFrom Matrix Matrix rowSums
#' @importFrom gridExtra arrangeGrob
#' @importFrom magrittr %>%
#' @importFrom grid grid.draw grid.newpage
#' @export
demultiplexTags <- function(tag_mtx,
                            init.cos.cut = 0.5,
                            max.iter = 30,
                            converge.threshold = 1e-3,
                            prob.cut = 0.5,
                            min.cell.fit = 10,
                            max.cell.fit = 1e4,
                            min.quantile.fit = 0.05, # Remove cells with total umi less than the specified quantile, which could be beads
                            max.quantile.fit = 0.95, # Remove cells with total umi greater than the specified quantile, which could be multiplets
                            residual.type = c("rqr", 'pearson'),
                            plot.umap = c("residual", "umi", "none"),
                            plot.diagnostics = F,
                            plot.path = getwd(),
                            plot.name = "",
                            umap.nn = 30L,
                            seed = 1,
                            point.size = 1,
                            label.size = 3,
                            min.tag.show = 50) {
    set.seed(seed)
    # TODO: Need to handle cells if rowSums = 0
    ## evaluate choices
    plot.umap <- match.arg(plot.umap)
    residual.type <- match.arg(residual.type)

    if(any(init.cos.cut < 0.5)) {
        cat("Warning: setting init.cos.cut less than 0.5 is not recommended.", fill=T)
    }

    zero_bc_cells = rowSums(tag_mtx) == 0
    if (sum(zero_bc_cells) > 0) {
        message(paste0("Detected ", sum(zero_bc_cells), " cells with 0 barcode count. These cells will not be classified."))
        tag_mtx = tag_mtx[!zero_bc_cells,]
    }

    barcodes <- colnames(tag_mtx)
    # Calculate cosine similarity matrix
    tag_mtx <- Matrix(as.matrix(tag_mtx), sparse = T)
    vec_norm2 <- apply(tag_mtx, 1, function(x) norm(x, type = "2")) # for each cell, calculate its L2 norm (Euclidean distance from origin)
    cos.umi_mtx <- tag_mtx / vec_norm2 # divide matrix by L2 norms to get cosine distances

    pr_mtx <- matrix(nrow = nrow(tag_mtx), ncol = ncol(tag_mtx), data = NA)
    rownames(pr_mtx) <- rownames(tag_mtx)
    colnames(pr_mtx) <- colnames(tag_mtx)
    rqr_mtx <- pr_mtx
    prob_mtx <- pr_mtx

    coef_list <- list()
    df_list<-list()
    for(bc in barcodes) {
        df <- data.frame(bc.umi = tag_mtx[, bc],
                         cos.umi = cos.umi_mtx[,bc])
        df$tt.umi <- rowSums(tag_mtx)

        cat("Running EM for ", bc, sep = "", fill = T)
        if(length(init.cos.cut) > 1) cos.cut = init.cos.cut[which(barcodes == bc)] else cos.cut = init.cos.cut
        res <- fit.em(df,init.cos.cut = cos.cut, converge.threshold = converge.threshold, max.iter = max.iter, min.cell.fit = min.cell.fit, max.cell.fit = max.cell.fit,
                      min.quantile.fit = min.quantile.fit, max.quantile.fit = max.quantile.fit)

        df <- res$df
        df_list[[bc]] <- df
        pr_mtx[,bc] <- df$pearson_residual
        rqr_mtx[,bc] <- df$rqr
        prob_mtx[,bc] <- df$post1
        coef_list[[bc]] <- res$fit0
    }


    # Recompute cosine metric on residual
    if(residual.type == "pearson") {
        res_mtx = pr_mtx
    } else if(residual.type == "rqr") {
        res_mtx = rqr_mtx
        res_mtx[is.na(res_mtx)] = 0 # Set to 0 for now
        max.rqr = max(res_mtx[is.finite(res_mtx)]) + 1 # Best to cut inf?
        res_mtx[res_mtx > max.rqr] = max.rqr
    }

    res_mtx <- Matrix(as.matrix(res_mtx), sparse = T)

    # Assign barcodes to cells
    call_mtx <- sapply(colnames(prob_mtx), function(bc) {
        prob_mtx[,bc] > prob.cut
    })
    barcode_count = rowSums(call_mtx)
    barcode_assign <- apply(call_mtx, 1, function(x) {
        if (sum(x) == 1) colnames(call_mtx)[which(x)] else NA
    })

    assign_table = data.frame(barcode_assign = barcode_assign,
                              barcode_count = barcode_count,
                              droplet_type = ifelse(barcode_count == 1, 'singlet', ifelse(barcode_count == 0, 'negative', 'multiplet')))
    assign_table$final_assign <- assign_table$barcode_assign
    assign_table$final_assign[is.na(assign_table$final_assign)] <- assign_table$droplet_type[is.na(assign_table$final_assign)]

    final_assign <- setNames(assign_table[, "final_assign"], rownames(assign_table))


    # UMAP Plotting Options
    if (plot.umap == "umi") { # Raw UMI Counts
        umap_res <- compute_umap(tag_mtx, use_dim = ncol(tag_mtx), n_component=2, n_neighbors = umap.nn)
    } else if (plot.umap == "residual") { # Pearson Residuals
        umap_res <- compute_umap(res_mtx, use_dim = ncol(res_mtx), n_component=2, n_neighbors = umap.nn)
    } else {
        umap_res <-  NA
    }

    glist <- list()

    if(plot.umap!="none") { # Always plot
        message("Plotting final umap...")
        umap_df <- cbind(umap_res, assign_table)
        umap_df$barcode_count = as.character(umap_df$barcode_count)
        umap_df$barcode_count[umap_df$barcode_count >= 3] = ">=3"
        umap_df$barcode_count <- factor(umap_df$barcode_count, c("0", "1", "2" ,">=3"))

        glist[["summary"]] <- plotSummary(umap_df, point.size = point.size, label.size = label.size, min.tag.show = min.tag.show)
    }

    # DEBUG
    # if(1) {
    #     bc == barcodes[1]
    #     df <- data.frame(
    #         bc.umi = tag_mtx[, bc],
    #         res = res_mtx[,bc],
    #         cos.umi = cos.umi_mtx[,bc],
    #         cos.res = cos.res_mtx[,bc])
    #     df$tt.umi <- rowSums(tag_mtx)
    #     df <- cbind(df, umap_df)
    #     df$is_negcell = 1:nrow(df) %in% neg_cell_list[[bc]]
    #     assign("df",df, env = .GlobalEnv)
    # }

    # Diagnostic Plots
    if (plot.diagnostics) {
        if(plot.umap=="none") stop("Please specify plot.umap.")
        for (bc in barcodes) {
            df <- data.frame(
                bc.umi = tag_mtx[, bc],
                res = res_mtx[,bc],
                cos.umi = cos.umi_mtx[,bc],
                prob.pos = prob_mtx[,bc])
            df$tt.umi <- rowSums(tag_mtx)
            df <- cbind(df, umap_df)
            mappings <- list(
                c("log(bc.umi)", "cos.umi", "cos.umi"),
                c("log(tt.umi)", "log(bc.umi)", "prob.pos"),
                c("log(tt.umi)", "log(tt.umi-bc.umi)", "prob.pos"),
                c("log(tt.umi)", "res", "prob.pos"),
                #c("qq"),
                #c("log(tt.umi)", "log(bc.umi)", "cos.umi"),
                #c("log(bc.umi)", "cos.umi", "prob.pos"),
                #c("log(tt.umi)", "cos.umi", "cos.umi"),
                #c("log(tt.umi)", "res", "cos.umi"),
                #c("log(res)", "res", "cos.res"),
                #c("log(res)", "cos.res", "cos.res"),
                #c("UMAP_1", "UMAP_2", "res"),
                c("UMAP_1", "UMAP_2", "prob.pos"),
                c("UMAP_1", "UMAP_2", "barcode_count")
            )
            #assign("df", df, env = .GlobalEnv)
            df$barcode_count[call_mtx[,bc] == 0] = NA
            glist[[bc]] <- plot.all.diagnostics(df, mappings, bc = bc, prob.cut = prob.cut, point.size = point.size, ncol=3)
        }
    }

    time <- (Sys.time() %>% make.names() %>% strsplit(split = "X") %>% unlist())[2]
    if(length(glist)) {
        pdf(paste0(plot.path, "/", plot.name, "_", time, "_assignment.pdf"), width = 15, height = 8)
        for(i in 1:length(glist)) {
            grid.newpage()
            message(paste0("Plotting ", names(glist)[i]))
            grid.draw(glist[[i]])
        }
        dev.off()
    }

    return(
        list(
            final_assign = final_assign,
            assign_table = assign_table,
            umap = umap_res %>% as.matrix(),
            res_mtx = res_mtx %>% as.matrix(),
            #rqr_mtx = rqr_mtx,
            #pr_mtx = pr_mtx,
            prob_mtx = prob_mtx,
            coefs = coef_list,
            df_list = df_list
        )
    )
}




# Compute rqr for poisson regression, function adapted from qresiduals {statmod}
rqr.poisson <- function(df, y, fit = "fit") {
    y = df[[y]]
    mu = df[[fit]]
    a <- ppois(y - 1, mu)
    b <- ppois(y, mu)
    u <- runif(n = length(y), min = a, max = b)
    rqr_pred <- qnorm(u)
}


# Compute rqr for negative binomial regression, function adapted from qresiduals {statmod}
rqr.nb <- function (df, y, fit = "fit", model)
{
    y = df[[y]]
    if (is.null(model$theta)) {
        size <- model$call$family[[2]]
    }
    else {
        size <- model$theta
    }
    mu <- df[[fit]]
    p <- size/(mu + size)
    a <- ifelse(y > 0, pbeta(p, size, pmax(y, 1)), 0)
    b <- pbeta(p, size, y + 1)
    u <- runif(n = length(y), min = a, max = b)
    qnorm(u)
}


#' @importFrom magrittr %>%
#' @importFrom dplyr summarize_at group_by_at
#' @importFrom gridExtra arrangeGrob
#' @importFrom ggplot2 ggplot aes aes_string scale_color_manual theme_bw ggtitle guides geom_text
plotSummary <- function(df, point.size = 1, label.size = 3 , min.tag.show = 50) {
    geom_point_plot <- ifelse("ggrastr" %in% rownames(installed.packages()), ggrastr::geom_point_rast, ggplot2::geom_point)

    unq_bcs = unique(df$barcode_assign)
    unq_bcs = unq_bcs[!is.na(unq_bcs)]
    use_color = get_factor_color(unq_bcs, "Set1")
    names(use_color) = unq_bcs

    g1 <- ggplot(df, aes_string("UMAP_1", "UMAP_2")) +
        geom_point_plot(aes_string(color = "barcode_assign"), stroke = 0, size = point.size) +
        #geom_point_plot(data = df[!is.na(df[["barcode_assign"]]), ], aes_string(color = "barcode_assign"), stroke = 0, size = .5) +
        scale_color_manual(values = use_color, na.value='lightgrey') +
        theme_bw() +
        ggtitle("barcode_assign") +
        guides(color = "none")
    show_bc = names(which(table(df[["barcode_assign"]]) > min.tag.show))
    label_data <- df %>% group_by_at("barcode_assign") %>% summarize_at(c("UMAP_1", "UMAP_2"), median)
    label_data <- label_data[label_data[["barcode_assign"]] %in% show_bc,,drop=F]
    g1 <- g1 + geom_text(
        aes_string(
            x="UMAP_1",y="UMAP_2",
            label = "barcode_assign"
        ),
        color = "black",
        size = label.size,
        data = label_data
    )
    count_color = get_gradient_color("BlueGreenRed", cnum = length(levels(df$barcode_count)))
    names(count_color) = levels(df$barcode_count)
    g2 <- ggplot(df, aes_string("UMAP_1", "UMAP_2")) +
        geom_point_plot(aes_string(color = "barcode_count"), stroke = 0, size = point.size) +
        scale_color_manual(values = count_color, na.value='lightgrey') +
        theme_bw() +
        ggtitle("barcode_count")
    return(arrangeGrob(g1,g2, ncol = 2))
}


#' @importFrom gridExtra arrangeGrob
#' @importFrom ggExtra ggMarginal
#' @importFrom ggplot2 ggplot aes aes_string scale_color_manual scale_color_gradientn theme_bw ggtitle guides geom_text stat_qq stat_qq_line xlab ylab geom_abline theme_classic labs
plot.all.diagnostics <- function(df, mappings, bc, prob.cut = 0.5, point.size = 1, ncol = 3) {
    geom_point_plot <- ifelse("ggrastr" %in% rownames(installed.packages()), ggrastr::geom_point_rast, ggplot2::geom_point)

    plot_list <- list()
    for(i in 1:length(mappings)){
        map <- mappings[[i]]
        if(map[1] == "qq") {
            # QQ plot
            p = ggplot(df[df$prob.pos < prob.cut,], aes(sample=res))+stat_qq(size = point.size, stroke = 0) +
                stat_qq_line() +
                xlab("Normal quantiles") +
                ylab("Residual quantiles") +
                theme_bw()
        } else if(map[1] == "qq_init") {
            p = ggplot(dd1, aes(x = theoretical, y = sample)) + geom_point(alpha = 0.3) +
                geom_abline(intercept = 0, slope = 1, color = "blue") +
                theme_classic() +
                ggtitle("QQ plot of standardized residuals")
        } else if (map[1] == "qq_final") {
            p = ggplot(dd2, aes(x = theoretical, y = sample)) + geom_point(alpha = 0.3) +
                geom_abline(intercept = 0, slope = 1, color = "blue") +
                theme_classic() +
                ggtitle("QQ plot of standardized residuals")
        }else {
            p <- ggplot(df, aes_string(map[1], map[2])) +
                geom_point_plot(aes_string(color = map[3]), stroke = 0, size = point.size) +
                ggtitle(bc) +
                theme_bw() +
                labs(color = gsub("_","\n",map[3]))

            count_color = get_gradient_color("BlueGreenRed", cnum = length(levels(df$barcode_count)))
            names(count_color) = levels(df$barcode_count)
            if(map[3] == "barcode_count") {
                p = p + scale_color_manual(values = count_color, na.value='lightgrey')
            } else {
                p = p + scale_color_gradientn(colors = get_gradient_color("BlueGreenRed"))
            }

            if(map[1]!= "UMAP_1"){
                p = ggMarginal(p, type = "histogram")
            }
        }
        plot_list[[i]] <- p
    }
    return(arrangeGrob(grobs = plot_list, ncol = ncol))
}







