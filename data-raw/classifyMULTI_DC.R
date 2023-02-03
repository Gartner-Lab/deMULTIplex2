





## This is a test version of the classifyMULTI function where I introduced the Pearson Residual from the positive singlet distribution and introduced some new plots
#' @export
classifyMULTI_DC <- function(barTable,
                             posThresh = 0.2,
                             cosineThresh = seq(0.5,0.9,0.1),
                             plotUMAP = "cumi.cos",
                             plotDiagnostics = T,
                             plotPath = getwd(),
                             UMAPNeighbors = 30L,
                             seed = 1) {
    require(MASS)
    require(ggrastr)
    require(grid)
    require(Matrix)
    require(gridExtra)
    require(ggpubr)
    require(ggExtra)

    set.seed(seed)

    bc_mtx <- barTable[, colnames(barTable) != "nUMI"]

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

    pr_mtx_pos <- pr_mtx

    # Matrix for storing corrected UMI counts (inverse NB fit, with library size set to median)
    cumi_mtx <- pr_mtx

    # NB fit
    barcodes <- colnames(bc_mtx)
    model_list <- list()

    for(bc in barcodes) {
        cat("Fitting GLM-NB for ", bc, sep = "", fill = T)

        df <- data.frame(Count = bc_mtx[, bc],
                         CosSim = cos_mtx_raw[,bc])
        df$nUMI <- rowSums(bc_mtx)

        # subset "negative" cells, i.e. below cosine similarity threshold
        neg_cells <- which(df$CosSim < max(df$CosSim) * (1-posThresh))

        model <- glm.nb(Count ~ log(nUMI), data = df[neg_cells,], link = log)
        model_list[[bc]] <- model

        # define "positive" i.e. Count = nUMI model
        pos_cells <- which(df$CosSim > max(df$CosSim) * (1-posThresh))
        # model_pos <- glm.nb(Count~nUMI, data = data.frame(Count = 1:100, nUMI = 1:100))
        model_pos <- glm.nb(Count ~ nUMI, data = df[pos_cells,], link = identity)

        predicted_res <- predict(model, df, type = "response", se.fit=TRUE)
        df_pred <- cbind(df, predicted_res)
        df_pred <- within(df_pred, {
            pred_mb_count <- fit
            LL <- fit - 1.96 * se.fit
            UL <- fit + 1.96 * se.fit
        })


        predicted_res_pos <- predict(model_pos, df, type = "response", se.fit=TRUE)
        df_pred_pos <- cbind(df, predicted_res_pos)
        df_pred_pos <- within(df_pred_pos, {
            pred_mb_count <- fit
            LL <- fit - 1.96 * se.fit
            UL <- fit + 1.96 * se.fit
        })


        # Compute Pearson residual using background count: (Count - fit) / sqrt(fit  + fit^2/theta)
        df_pred$pearson_residual <- (df_pred$Count-df_pred$fit)/sqrt(df_pred$fit + df_pred$fit^2/model$theta)

        df_pred_pos$pearson_residual <- (df_pred_pos$Count-df_pred_pos$fit)/sqrt(df_pred_pos$fit + df_pred_pos$fit^2/model_pos$theta)


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

        pr_mtx_pos[,bc] <- df_pred_pos$pearson_residual
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
                PearsonRes_Pos = pr_mtx_pos[,bc],
                CosSim_RawUMI = cos_mtx_raw[,bc],
                CosSim_PearsonRes = cos_mtx_pr[,bc])
            df$nUMI_Raw <- rowSums(bc_mtx)
            df$nUMI_Corrected <- rowSums(cumi_mtx)

            mappings <- list(
                c("log(nUMI_Raw)", "log(RawUMI)", "CosSim_RawUMI"),
                c("log(nUMI_Raw)", "PearsonRes", "CosSim_RawUMI"),
                c("log(RawUMI)", "CosSim_RawUMI", "CosSim_RawUMI"),
                # c("log(PearsonRes)", "CosSim_PearsonRes", "CosSim_PearsonRes"),
                c("log(nUMI_Raw)", "log(RawUMI)", "PearsonRes_Pos"),
                c("log(nUMI_Raw)", "PearsonRes_Pos", "PearsonRes_Pos"),
                c("log(RawUMI)", "CosSim_RawUMI", "PearsonRes_Pos")
            )

            if (plotUMAP != "None") {
                df <- cbind(df, umap_res)
                mappings <- append(mappings, list(c("UMAP_1", "UMAP_2", "CosSim_RawUMI"),
                                                  c("UMAP_1", "UMAP_2", "CosSim_PearsonRes")))
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
            #   grobs = list(ggplot(df, aes(x = CosSim_RawUMI)) + geom_histogram(bins = 50, color = "white") + theme_bw() + theme(plot.margin = unit(c(0,1,0,0), "cm")),
            #                ggplot(df, aes(x = CosSim_PearsonRes)) + geom_histogram(bins = 50, color = "white") + theme_bw() + theme(plot.margin = unit(c(0,1,0,0), "cm"))),
            #   nrow = 2)

            plot_list[[length(mappings) + 1]] <- arrangeGrob(
                grobs = list(ggplot(df, aes(x = CosSim_RawUMI)) + geom_histogram(bins = 50, color = "white") + theme_bw() + theme(plot.margin = unit(c(0,1,0,0), "cm")),
                             ggplot(df, aes(x = PearsonRes_Pos)) + geom_histogram(bins = 50, color = "white") + theme_bw() + theme(plot.margin = unit(c(0,1,0,0), "cm"))),
                nrow = 2)


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
            cos_mtx_pr = cos_mtx_pr %>% as.matrix(),
            models = model_list
        )
    )
}






