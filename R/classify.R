


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
                          model = c("nb", "poisson"),
                          posThresh = 0.2,
                          cosineThresh = seq(0.5,0.9,0.1),
                          plotUMAP = "cumi.cos",
                          plotDiagnostics = T,
                          plotPath = getwd(),
                          UMAPNeighbors = 30L,
                          seed = 1,
                          gini.cut = 0.4) {
    set.seed(seed)
    # TODO: Need to handle cells if rowSums = 0
    ## evaluate choices
    model <- match.arg(model)
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
        #assign("df", df, env=.GlobalEnv)
        #assign("neg_cells", neg_cells, env=.GlobalEnv)
        if(model == "nb") {
            fit.res <- glm.nb.fit(df, use_cells = neg_cells, x = "log(nUMI)", y = "RawUMI")
        } else if(model == "poisson") {
            fit.res <- glm.poisson.fit(df, use_cells = neg_cells, x = "log(nUMI)", y = "RawUMI")
        }

        model_list[[bc]] <- fit.res$model
        df_pred <- fit.res$df.pred

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
            df <- cbind(df, df_pred[,c("fit", "se.fit")])
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

    # Compute gini as a potential way to call doublets
    cos_for_gini = cos_mtx_pr
    cos_for_gini[cos_for_gini < 0.5] = 0
    gini_res<- apply(cos_for_gini, 1, reldist::gini )

    return(
        list(
            res = assign_tbl %>% as.matrix(),
            umap = umap_res %>% as.matrix(),
            pr_mtx = pr_mtx %>% as.matrix(),
            cumi_mtx = cumi_mtx %>% as.matrix(),
            cos_mtx_raw = cos_mtx_raw %>% as.matrix(),
            cos_mtx_pr = cos_mtx_pr %>% as.matrix(),
            gini_res = gini_res,
            models = model_list
        )
    )
}



glm.nb.fit <- function(df, use_cells, x, y) {
    model <- glm.nb(as.formula(paste0(y,"~",x)), data = df[use_cells,], link = log)

    predicted_res <- predict(model, df, type = "response", se.fit=TRUE)
    df_pred <- cbind(df, predicted_res)
    # Compute Pearson residual using background count: (Count - fit) / sqrt(fit  + fit^2/theta)
    df_pred$pearson_residual <- (df_pred[[y]]-df_pred$fit)/sqrt(df_pred$fit + df_pred$fit^2/model$theta)

    # Optional: return corrected UMI counts (not necessary for demultiplexing)
    med_scale_df <- data.frame(nUMI = median(df_pred$nUMI))
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



glm.poisson.fit <- function(df, use_cells, x, y) {
    model <- glm(as.formula(paste0(y,"~",x)), fam = poisson(link = log), data = df[use_cells,])
    predicted_res <- predict(model, df, type = "response", se.fit=TRUE)
    df_pred <- cbind(df, predicted_res)
    # Compute Pearson residual using background count: (Count - fit) / sqrt(fit  + fit^2/theta)
    df_pred$pearson_residual <- (df_pred[[y]]-df_pred$fit)/sqrt(df_pred$fit)

    # Optional: return corrected UMI counts (not necessary for demultiplexing)
    med_scale_df <- data.frame(nUMI = median(df_pred$nUMI))
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





