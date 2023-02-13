


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
classifyMULTI <- function(bar.table,
                          model = c("nb", "poisson", "poisson.fixed.b1", "poisson.analytical"),
                          pos.thresh = 0.5,
                          cos.thresh = seq(0.5,0.9,0.1),
                          pr.thresh = seq(2,10,2),
                          plot.umap = "cos.cumi",
                          plot.diagnostics = T,
                          plot.path = getwd(),
                          umap.nn = 30L,
                          seed = 1,
                          gini.cut = 0.4) {
    set.seed(seed)
    # TODO: Need to handle cells if rowSums = 0
    ## evaluate choices
    model <- match.arg(model)
    bc_mtx <- bar.table

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


    cumi_mtx <- pr_mtx # Matrix for storing corrected UMI counts (inverse NB fit, with library size set to median)
    prob_mtx <- pr_mtx # Matrix for storing positive probability based on fit

    # NB fit
    barcodes <- colnames(bc_mtx)
    coef_list <- list()

    for(bc in barcodes) {
        cat("Fitting GLM-NB for ", bc, sep = "", fill = T)
        df <- data.frame(bc.umi = bc_mtx[, bc],
                         cos.umi = cos.umi_mtx[,bc])
        df$tt.umi <- rowSums(bc_mtx)
        # subset "negative" cells, i.e. below cosine similarity threshold
        neg_cells <- which(df$cos.umi < max(df$cos.umi, na.rm = T) * (1-pos.thresh))
        if(bc == barcodes[1]) { # For debugging
            assign("df", df, env=.GlobalEnv)
            assign("neg_cells", neg_cells, env=.GlobalEnv)
        }

        if(model == "nb") {
            fit.res <- glm.nb.fit(df, neg_cells = neg_cells, x = "log(tt.umi)", y = "bc.umi")
        } else if(model == "poisson") {
            fit.res <- glm.poisson.fit(df, neg_cells = neg_cells, x = "log(tt.umi)", y = "bc.umi")
        } else if(model == "poisson.fixed.b1") {
            fit.res <- glm.poisson.fit(df, neg_cells = neg_cells, x = "offset(log(tt.umi))", y = "bc.umi")
        }
        else if(model == "poisson.analytical") { # Same as Poisson fixed! Remove later
            fit.res <- glm.poisson.analytical(df, neg_cells = neg_cells, x = "offset(log(tt.umi))", y = "bc.umi")
        }

        coef_list[[bc]] <- fit.res$model$coefficients
        df_pred <- fit.res$df.pred

        cumi_mtx[,bc] <- df_pred$corrected_umi
        pr_mtx[,bc] <- df_pred$pearson_residual
        prob_mtx[,bc] <- df_pred$prob_pos
    }


    # Recompute cosine metric
    pr_mtx <- Matrix(as.matrix(pr_mtx), sparse = T)
    vec_norm2 <- apply(pr_mtx, 1, function(x) norm(x, type = "2"))
    cos.pr_mtx <- pr_mtx / vec_norm2


    # Assign singlets based on cosine score on Pearson residual
    assign_cos.pr <- matrix(nrow = nrow(bc_mtx), ncol = length(cos.thresh), data = NA)
    colnames(assign_cos.pr) <- paste0("assign_cos.pr_", cos.thresh)
    rownames(assign_cos.pr) <- rownames(bc_mtx)
    for(i in 1:length(cos.thresh)) {
        assign_mtx <- cos.pr_mtx > cos.thresh[i]
        assign_res <- apply(assign_mtx, 1, function(x) {
            if (length(which(x)) == 1) {
                barcodes[which(x)]
            } else {
                NA
            }
        })
        assign_cos.pr[,i] <- assign_res
    }

    # Assign singlets based on Pearson residual
    assign_pr <- matrix(nrow = nrow(bc_mtx), ncol = length(pr.thresh), data = NA)
    colnames(assign_pr) <- paste0("assign_pr_", pr.thresh)
    rownames(assign_pr) <- rownames(bc_mtx)
    for(i in 1:length(pr.thresh)) {
        assign_mtx <- pr_mtx > pr.thresh[i]
        assign_res <- apply(assign_mtx, 1, function(x) {
            if (length(which(x)) == 1) {
                barcodes[which(x)]
            } else {
                NA
            }
        })
        assign_pr[,i] <- assign_res
    }


    # UMAP Plotting Options
    umap_res <-  NA
    if (plot.umap == "UMI") { # Raw UMI Counts
        umap_res <- compute_umap(bc_mtx, use_dim = ncol(bc_mtx), n_component=2, n_neighbors = umap.nn)
    } else if (plot.umap == "pr") { # Pearson Residuals
        umap_res <- compute_umap(pr_mtx, use_dim = ncol(pr_mtx), n_component=2, n_neighbors = umap.nn)
    } else if (plot.umap == "cos.pr") { # Cosine Similarity calculated from Pearson Residuals
        umap_res <- compute_umap(cos.pr_mtx, use_dim = ncol(cos.pr_mtx), n_component=2, n_neighbors = umap.nn)
    } else if (plot.umap == "cumi") { # Corrected UMI Counts
        umap_res <- compute_umap(cumi_mtx, use_dim = ncol(cumi_mtx), n_component=2, n_neighbors = umap.nn)
    } else if (plot.umap == "cos.cumi") { # Cosine Similarity calculated from Corrected UMI Counts
        vec_norm2 <- apply(cumi_mtx, 1, function(x) norm(x, type = "2"))
        cos.cumi_mtx <- cumi_mtx  / vec_norm2
        umap_res <- compute_umap(cos.cumi_mtx,use_dim = ncol(cos.cumi_mtx), n_component=2, n_neighbors = umap.nn)
    } else if (plot.umap == "none") {
    } else { stop("plot.umap = ", plot.umap, " is not valid") }

    # Diagnostic Plots
    if (plot.diagnostics) {
        glist <- list()
        for (bc in barcodes) {
            df <- data.frame(
                bc.umi = bc_mtx[, bc],
                cr.umi = cumi_mtx[, bc], # For now always set to TRUE so is computed
                pr = pr_mtx[,bc],
                cos.umi = cos.umi_mtx[,bc],
                cos.pr = cos.pr_mtx[,bc])
            df$tt.umi <- rowSums(bc_mtx)
            df$tt.cr.umi <- rowSums(cumi_mtx)
            df <- cbind(df, df_pred[,c("fit"), drop=FALSE])
            # if(bc == barcodes[1]) {
            #     assign("df",df, env = .GlobalEnv)
            # }

            mappings <- list(
                c("log(tt.umi)", "log(bc.umi)", "cos.umi"),
                c("log(bc.umi)", "cos.umi", "cos.umi"),
                c("log(tt.umi)", "pr", "cos.umi"),
                c("tt.cr.umi", "cr.umi", "cos.pr"),
                c("log(pr)", "cos.pr", "cos.pr"),
                c("log(tt.cr.umi)", "pr", "cos.pr"),
                c("tt.umi", "bc.umi", "cos.umi"),
                #c("fit", "bc.umi", "cos.pr"),
                #c("fit", "pr", "cos.pr"),
                #c("log(fit)", "log(tt.umi)", "cos.pr"),
                c("log(tt.cr.umi)", "log(cr.umi)", "cos.pr")
            )

            if (plot.umap != "none") {
                df <- cbind(df, umap_res)
                mappings <- append(mappings, list(c("UMAP_1", "UMAP_2", "cos.pr")))
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
            #     grobs = list(ggplot(df, aes(x = cos.umi)) + geom_histogram(bins = 50, color = "white") + theme_bw() + theme(plot.margin = unit(c(0,1,0,0), "cm")),
            #                  ggplot(df, aes(x = cos.pr)) + geom_histogram(bins = 50, color = "white") + theme_bw() + theme(plot.margin = unit(c(0,1,0,0), "cm"))),
            #     nrow = 2)

            glist[[bc]] <- arrangeGrob(grobs = plot_list, ncol = 3)
        }

        time <- (Sys.time() %>% make.names() %>% strsplit(split = "X") %>% unlist())[2]
        pdf(paste0(plot.path, "/", time, "_diagnostics.pdf"), width = 12, height = 10)
        for(i in 1:length(glist)) {
            grid.newpage()
            grid.draw(glist[[i]])
        }
        dev.off()
    }

    # Compute gini as a potential way to call doublets
    cos_for_gini = cos.pr_mtx
    cos_for_gini[cos_for_gini < 0.5] = 0
    gini_res<- apply(cos_for_gini, 1, reldist::gini )

    return(
        list(
            res = cbind(assign_pr, assign_cos.pr) %>% as.matrix(),
            umap = umap_res %>% as.matrix(),
            pr_mtx = pr_mtx %>% as.matrix(),
            prob_mtx = prob_mtx %>% as.matrix(),
            cumi_mtx = cumi_mtx %>% as.matrix(),
            cos.umi_mtx = cos.umi_mtx %>% as.matrix(),
            cos.pr_mtx = cos.pr_mtx %>% as.matrix(),
            gini_res = gini_res,
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
    # Compute truncated probability of being positive for each barcode
    df_pred$prob_neg = dpois(df_pred$bc.umi, lambda = exp(df_pred$fit))
    df_pred$prob_neg[df_pred$pr < 0] = 1
    df_pred$prob_pos = 1- df_pred$prob_neg

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


# Same as Poisson fixed! Remove later
#' @export
glm.poisson.analytical <- function(df, neg_cells, x, y) {
    beta1 = 1
    beta0 = log(sum(df[neg_cells,]$bc.umi) / sum(df[neg_cells,]$tt.umi))
    df_pred <- df
    df_pred$fit <- exp(beta0 + beta1 * log(df$tt.umi))

    # Compute Pearson residual using background count: (Count - fit) / sqrt(fit  + fit^2/theta)
    df_pred$pearson_residual <- (df_pred[[y]]-df_pred$fit)/sqrt(df_pred$fit)
    # df_pred$response_residual <- df_pred[[y]]-df_pred$fit
    # df_pred$working_residual <- (df_pred[[y]]-df_pred$fit)/ df_pred$fit
    # Deviance residual? https://www.datascienceblog.net/post/machine-learning/interpreting_generalized_linear_models/

    # Optional: return corrected UMI counts (not necessary for demultiplexing)
    med_scale_df <- data.frame(tt.umi = median(df_pred$tt.umi))
    mu <- as.numeric(exp(beta0 + beta1 * log(med_scale_df)))
    variance <- mu
    df_pred$corrected_count <- mu + df_pred$pearson_residual * sqrt(variance)
    df_pred$corrected_umi <- round(df_pred$corrected_count, 0)
    df_pred$corrected_umi[df_pred$corrected_umi < 0] <- 0
    return(
        list(
            model = list(
                coefficients = c(beta0 = beta0, beta1 = beta1)
            ),
            df.pred = df_pred
        )
    )
}

# Compute rqr for poisson regression
#' @export
rqr.poisson <- function(df, y, fit = "fit") {
    y = df[[y]]
    mu = df[[fit]]
    a <- ppois(y - 1, mu)
    b <- ppois(y, mu)
    u <- runif(n = length(y), min = a, max = b)
    rqr_pred <- qnorm(u)
}



