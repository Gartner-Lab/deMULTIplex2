

# Description:
# The new demultiplex function
#' @importFrom Matrix Matrix rowSums
#' @importFrom MASS glm.nb
#' @importFrom ggrastr geom_point_rast
#' @importFrom gridExtra arrangeGrob
#' @importFrom magrittr %>%
#' @export
classify.cells <- function(bc_mtx,
                          model = c("nb", "poisson", "poisson.analytical"),
                          neg.thresh = 0.2,
                          residual.type = c("rqr", 'pearson'),
                          residual.thresh = 1e-4, # Re-design for compatibility
                          plot.umap = c("umi", "residual"),
                          plot.diagnostics = T,
                          plot.path = getwd(),
                          umap.nn = 30L,
                          seed = 1) {
    set.seed(seed)
    # TODO: Need to handle cells if rowSums = 0
    ## evaluate choices
    model <- match.arg(model)
    plot.umap <- match.arg(plot.umap)
    residual.type <- match.arg(residual.type)

    if (any(rowSums(bc_mtx)) == 0) {
        stop("Please remove any cells with 0 total barcode counts")
    }


    # Calculate cosine similarity matrix
    bc_mtx <- Matrix(as.matrix(bc_mtx), sparse = T)
    vec_norm2 <- apply(bc_mtx, 1, function(x) norm(x, type = "2")) # for each cell, calculate its L2 norm (Euclidean distance from origin)
    cos.umi_mtx <- bc_mtx / vec_norm2 # divide matrix by L2 norms to get cosine distances

    # Matrix for storing Pearson residual (from NB fit)
    pr_mtx <- matrix(nrow = nrow(bc_mtx), ncol = ncol(bc_mtx), data = NA)
    rownames(pr_mtx) <- rownames(bc_mtx)
    colnames(pr_mtx) <- colnames(bc_mtx)
    rqr_mtx <- pr_mtx

    # NB fit
    barcodes <- colnames(bc_mtx)
    coef_list <- list()
    res_cuts <- list()
    res_fits <- list()

    for(bc in barcodes) {
        df <- data.frame(bc.umi = bc_mtx[, bc],
                         cos.umi = cos.umi_mtx[,bc])
        df$tt.umi <- rowSums(bc_mtx)
        # subset "negative" cells, i.e. below cosine similarity threshold
        neg_cells <- which(df$cos.umi < max(df$cos.umi, na.rm = T) * neg.thresh)
        # if(bc == barcodes[1]) { # For debugging
        #     assign("df", df, env=.GlobalEnv)
        #     assign("neg_cells", neg_cells, env=.GlobalEnv)
        # }

        if(model == "nb") {
            cat("Fitting GLM-NB for ", bc, sep = "", fill = T)
            fit.res <- glm.nb.fit(df, neg_cells = neg_cells, x = "log(tt.umi)", y = "bc.umi")
        } else if(model == "poisson") {
            cat("Fitting GLM-poisson for ", bc, sep = "", fill = T)
            fit.res <- glm.poisson.fit(df, neg_cells = neg_cells, x = "log(tt.umi)", y = "bc.umi")
        } else if(model == "poisson.analytical") {
            cat("Fitting GLM-poisson for ", bc, sep = "", fill = T)
            fit.res <- glm.poisson.fit(df, neg_cells = neg_cells, x = "offset(log(tt.umi))", y = "bc.umi") # Fixing b1 = 1 also fixes b0
        }

        coef_list[[bc]] <- fit.res$model$coefficients
        df_pred <- fit.res$df.pred
        pr_mtx[,bc] <- df_pred$pearson_residual
        rqr_mtx[,bc] <- df_pred$rqr

        x = df_pred$rqr[neg_cells]
        x = x[is.finite(x) & !is.na(x)]
        res_fit <- fitdistr(x, "normal")
        res_cut <- qnorm(residual.thresh, res_fit$estimate[1], res_fit$estimate[2], lower.tail = F)
        max_cut = abs(max(df_pred$rqr[is.finite(df_pred$rqr)], na.rm=T))
        res_cut <- min(abs(res_cut), max_cut)
        res_fits[[bc]] = res_fit
        res_cuts[[bc]] <- res_cut
    }


    # Recompute cosine metric on residual
    if(residual.type == "pr") {
        res_mtx = pr_mtx
    } else if(residual.type == "rqr") {
        res_mtx = rqr_mtx
        res_mtx[is.na(res_mtx)] = 0 # Set to 0 for now
        max.rqr = max(res_mtx[is.finite(res_mtx)]) + 1 # Best to cut inf?
        res_mtx[res_mtx > max.rqr] = max.rqr
    }

    res_mtx <- Matrix(as.matrix(res_mtx), sparse = T)
    vec_norm2 <- apply(res_mtx, 1, function(x) norm(x, type = "2"))
    cos.res_mtx <- res_mtx / vec_norm2

    # Assign barcodes to cells
    res_cuts <- unlist(res_cuts)
    names(res_cuts) <- colnames(res_mtx)

    res_call_mtx <- sapply(colnames(res_mtx), function(bc) {
        res_mtx[,bc] > res_cuts[bc]
    })
    barcode_count = rowSums(res_call_mtx)
    barcode_assign <- apply(res_call_mtx, 1, function(x) {
        if (sum(x) == 1) colnames(res_call_mtx)[which(x)] else NA
    })
    assign_table = data.frame(barcode_assign = barcode_assign, barcode_count = barcode_count)


    # UMAP Plotting Options
    if (plot.umap == "umi") { # Raw UMI Counts
        umap_res <- compute_umap(bc_mtx, use_dim = ncol(bc_mtx), n_component=2, n_neighbors = umap.nn)
    } else if (plot.umap == "residual") { # Pearson Residuals
        umap_res <- compute_umap(pr_mtx, use_dim = ncol(res_mtx), n_component=2, n_neighbors = umap.nn)
    } else {
        umap_res <-  NA
    }

    glist <- list()

    if(1) { # Always plot
        umap_df <- cbind(umap_res, assign_table)
        umap_df$barcode_count = as.character(umap_df$barcode_count)
        umap_df$barcode_count[umap_df$barcode_count >= 3] = ">=3"
        umap_df$barcode_count <- factor(umap_df$barcode_count, c("0", "1", "2" ,">=3"))

        unq_bcs = unique(umap_df$barcode_assign)
        unq_bcs = unq_bcs[!is.na(unq_bcs)]
        use_color = get_factor_color(unq_bcs, "Set1")
        names(use_color) = unq_bcs
        g1 <- ggplot(umap_df, aes_string("UMAP_1", "UMAP_2")) +
            geom_point_rast(data = umap_df[is.na(umap_df[["barcode_assign"]]), ], aes_string(color = "barcode_assign"), stroke = 0, size = 1) +
            geom_point_rast(data = umap_df[!is.na(umap_df[["barcode_assign"]]), ], aes_string(color = "barcode_assign"), stroke = 0, size = 1) +
            scale_color_manual(values = use_color, na.value='lightgrey') +
            theme_bw() +
            ggtitle("barcode_assign") +
            guides(color = "none")
        label_data <- umap_df %>% group_by_at("barcode_assign") %>% summarize_at(c("UMAP_1", "UMAP_2"), median)
        g1 <- g1 + geom_label(
            aes_string(
                x="UMAP_1",y="UMAP_2",
                label = "barcode_assign"
            ),
            color = "black",
            size = 2,
            data = label_data
        )

        g2 <- ggplot(umap_df, aes_string("UMAP_1", "UMAP_2")) +
            geom_point_rast(aes_string(color = "barcode_count"), stroke = 0, size = 1) +
            scale_color_manual(values = get_numeric_color("BlueGreenRed", cnum = length(levels(umap_df$barcode_count))), na.value='lightgrey') +
            theme_bw() +
            ggtitle("barcode_count")
        glist[["summary"]] <- arrangeGrob(g1,g2, ncol = 2)
    }


    # Diagnostic Plots
    if (plot.diagnostics) {
        for (bc in barcodes) {
            df <- data.frame(
                bc.umi = bc_mtx[, bc],
                res = res_mtx[,bc],
                cos.umi = cos.umi_mtx[,bc],
                cos.res = cos.res_mtx[,bc])
            df$tt.umi <- rowSums(bc_mtx)
            df <- cbind(df, df_pred[,c("fit"), drop=FALSE])
            df <- cbind(df, umap_res)
            # if(bc == barcodes[1]) {
            #     assign("df",df, env = .GlobalEnv)
            # }

            mappings <- list(
                c("log(tt.umi)", "log(bc.umi)", "cos.umi"),
                c("log(bc.umi)", "cos.umi", "cos.umi"),
                c("log(tt.umi)", "res", "cos.umi"),
                c("log(res)", "cos.res", "cos.res"),
                c("UMAP_1", "UMAP_2", "res"),
                c("UMAP_1", "UMAP_2", "cos.res")
            )

            plot_list <- list()
            for(i in 1:length(mappings)){
                map <- mappings[[i]]
                p <- ggplot(df, aes_string(map[1], map[2])) +
                    geom_point_rast(aes_string(color = map[3]), stroke = 0, size = 1) +
                    scale_color_gradientn(colors = get_numeric_color("BlueGreenRed")) +
                    ggtitle(bc) +
                    theme_bw() +
                    labs(color = gsub("_","\n",map[3]))

                if(map[2] == "res") {
                    p = p +
                        geom_hline(yintercept = c(res_cuts[bc])) +
                        geom_hline(yintercept = res_fits[[bc]]$estimate[1], color = "red")
                }
                p = ggMarginal(p, type = "histogram")

                plot_list[[i]] <- p
            }

            glist[[bc]] <- arrangeGrob(grobs = plot_list, ncol = 3)
        }
    }

    time <- (Sys.time() %>% make.names() %>% strsplit(split = "X") %>% unlist())[2]
    pdf(paste0(plot.path, "/", time, "_assignment.pdf"), width = 15, height = 10)
    for(i in 1:length(glist)) {
        grid.newpage()
        grid.draw(glist[[i]])
    }
    dev.off()

    return(
        list(
            assign_table = assign_table,
            umap = umap_res %>% as.matrix(),
            res_mtx = res_mtx %>% as.matrix(),
            rqr_mtx = rqr_mtx,
            res_cuts = res_cuts,
            res_fits = res_fits,
            coefs = coef_list
        )
    )
}


#' @export
glm.nb.fit <- function(df, neg_cells, x, y) {
    model <- glm.nb(as.formula(paste0(y,"~",x)), data = df[neg_cells,], link = log)

    predicted_res <- predict(model, df, type = "response", se.fit=TRUE)
    df_pred <- cbind(df, predicted_res)
    # Compute Pearson residual using background count: (Count - fit) / sqrt(fit  + fit^2/theta)
    df_pred$pearson_residual <- (df_pred[[y]]-df_pred$fit)/sqrt(df_pred$fit + df_pred$fit^2/model$theta)
    df_pred$rqr <- rqr.nb(df_pred, y=y,fit="fit", model = model)

    # Optional: return corrected UMI counts (not necessary for demultiplexing)
    med_scale_df <- data.frame(tt.umi = median(df_pred$tt.umi))
    expected_mu <- predict(model, med_scale_df, type = "response", se.fit=TRUE)
    mu <- expected_mu$fit
    theta <- model$theta
    variance <- mu + mu^2 / theta
    df_pred$corrected_count <- mu + df_pred$pearson_residual * sqrt(variance)
    df_pred$corrected_umi <- round(df_pred$corrected_count, 0)
    df_pred$corrected_umi[df_pred$corrected_umi < 0] <- 0
    return(
        list(
            model = model,
            df.pred = df_pred
        )
    )
}


#' @export
glm.poisson.fit <- function(df, neg_cells, x, y) {
    model <- glm(as.formula(paste0(y,"~",x)), fam = poisson(link = log), data = df[neg_cells,])
    predicted_res <- predict(model, df, type = "response", se.fit=TRUE)
    df_pred <- cbind(df, predicted_res)
    # Compute Pearson residual using background count: (Count - fit) / sqrt(fit  + fit^2/theta)
    df_pred$pearson_residual <- (df_pred[[y]]-df_pred$fit)/sqrt(df_pred$fit)
    df_pred$rqr <- rqr.poisson(df_pred, y=y,fit="fit")

    # Optional: return corrected UMI counts (not necessary for demultiplexing)
    med_scale_df <- data.frame(tt.umi = median(df_pred$tt.umi))
    expected_mu <- predict(model, med_scale_df, type = "response", se.fit=TRUE)
    mu <- expected_mu$fit
    variance <- mu
    df_pred$corrected_count <- mu + df_pred$pearson_residual * sqrt(variance)
    df_pred$corrected_umi <- round(df_pred$corrected_count, 0)
    df_pred$corrected_umi[df_pred$corrected_umi < 0] <- 0
    return(
        list(
            model = model,
            df.pred = df_pred
        )
    )
}



# Compute rqr for poisson regression, function adapted from qresiduals {statmod}
#' @export
rqr.poisson <- function(df, y, fit = "fit") {
    y = df[[y]]
    mu = df[[fit]]
    a <- ppois(y - 1, mu)
    b <- ppois(y, mu)
    u <- runif(n = length(y), min = a, max = b)
    rqr_pred <- qnorm(u)
}


# Compute rqr for negative binomial regression, function adapted from qresiduals {statmod}
#' @export
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


