
require(MASS)
require(ggrastr)
require(grid)
require(gridExtra)
require(ggpubr)
require(Matrix)



# The new demultiplex function with "Correct" Pearson residual and cosine similarity

multi_align <- function(bc_mtx, inter.num = 1, initial.pos = .2, cosine.cut = seq(0.5,0.9, 0.1), plot.umap = "cumi.cos", plot.diagnostics = TRUE, plot.path = "./", umap.neighbors = 30L, seed = 1) {
    set.seed(seed)
    cc_vec <- diag(x=1, nrow=ncol(bc_mtx), ncol = ncol(bc_mtx))
    colnames(cc_vec) <- colnames(bc_mtx); rownames(cc_vec) <- colnames(bc_mtx)

    bc_mtx <- Matrix(as.matrix(bc_mtx), sparse = T)
    dot_res <- bc_mtx %*% cc_vec
    vec_norm2 <- apply(bc_mtx, 1, function(x) norm(x, type = "2"))
    cos_mtx_raw <- dot_res / vec_norm2

    # Matrix for storing Pearson residual (from NB fit)
    pr_mtx = matrix(nrow = nrow(bc_mtx), ncol = ncol(bc_mtx), data = NA)
    rownames(pr_mtx) = rownames(bc_mtx)
    colnames(pr_mtx) = colnames(bc_mtx)

    # Matrix for storing corrected umi counts (inverse NB fit, with library size set to median)
    cumi_mtx = pr_mtx

    # NB fit
    total_bc_count = rowSums(bc_mtx)
    if(any(total_bc_count) == 0) {
        stop("Please remove any cells with 0 barcode count in the bc_mtx.")
    }
    all_bcs = colnames(bc_mtx)
    model_list <- list()
    for(cur_barcode  in all_bcs) {
        message(paste0("Fitting GLM-NB for ", cur_barcode))
        cur_df = data.frame(count = bc_mtx[, cur_barcode], cos_sim = cos_mtx_raw[,cur_barcode])
        cur_df$total_bc_count = total_bc_count
        neg_cells = which(cur_df$cos_sim < max(cur_df$cos_sim) * (1-initial.pos))
        #assign("cur_df", cur_df, env = .GlobalEnv)
        m1 <- glm.nb(count ~ log(total_bc_count), data = cur_df[neg_cells,], link = log)
        model_list[[cur_barcode]] = m1


        predicted_res <- predict(m1, cur_df, type = "response", se.fit=TRUE)
        df_pred <- cbind(cur_df, predicted_res)

        #
        df_pred <- within(df_pred, {
            pred_mb_count <- fit
            LL <- fit - 1.96 * se.fit
            UL <- fit + 1.96 * se.fit
        })
        # Compute Pearson residual using background count
        df_pred$pearson_residual = (df_pred$count-df_pred$fit)/sqrt(df_pred$fit + df_pred$fit^2/m1$theta)

        # If user wants corrected umi to be returned, which is not necessary for demultiplexing purpose, for now always compute
        med_scale_df = data.frame(total_bc_count = median(df_pred$total_bc_count))
        expected_mu <- predict(m1, med_scale_df, type = "response", se.fit=TRUE)
        mu = expected_mu$fit
        theta = m1$theta
        variance <- mu + mu^2 / theta
        df_pred$corrected_count <- mu + df_pred$pearson_residual * sqrt(variance)
        df_pred$corrected_umi <- round(df_pred$corrected_count, 0)
        df_pred$corrected_umi[df_pred$corrected_umi <0] = 0
        cumi_mtx[,cur_barcode] = df_pred$corrected_umi
        pr_mtx[,cur_barcode] = df_pred$pearson_residual
    }

    # Recompute cosine metric
    pr_mtx <- Matrix(as.matrix(pr_mtx), sparse = T)
    dot_res <- pr_mtx %*% cc_vec
    vec_norm2 <- apply(pr_mtx, 1, function(x) norm(x, type = "2"))
    cos_mtx_pr <- dot_res / vec_norm2


    # Assign singlets
    assign_tbl = matrix(nrow = nrow(bc_mtx), ncol = length(cosine.cut), data = NA)
    colnames(assign_tbl) <- paste0("assign.res_", cosine.cut)
    rownames(assign_tbl) <- rownames(bc_mtx)
    for(i in 1:length(cosine.cut)) {
        assign_mtx = cos_mtx_pr > cosine.cut[i]
        assign_res = apply(assign_mtx, 1, function(x) {if(length(which(x)) == 1) all_bcs[which(x)] else NA})
        assign_tbl[,i] = assign_res
    }

    umap_res = NA
    if(plot.umap!="none") {
        if(plot.umap == "pr") {
            umap_res = compute_umap(pr_mtx,use_dim = ncol(pr_mtx), n_component=2, n_neighbors = umap.neighbors)
        }
        else if(plot.umap == "pr.cos") {
            umap_res = compute_umap(cos_mtx_pr,use_dim = ncol(cos_mtx_pr), n_component=2, n_neighbors = umap.neighbors)
        }
        else if(plot.umap == "cumi") {
            umap_res = compute_umap(cumi_mtx,use_dim = ncol(cumi_mtx), n_component=2, n_neighbors = umap.neighbors)
        } else if (plot.umap == "cumi.cos") {
            dot_res <- cumi_mtx %*% cc_vec
            vec_norm2 <- apply(cumi_mtx, 1, function(x) norm(x, type = "2"))
            cos_mtx_cumi <- dot_res / vec_norm2
            umap_res = compute_umap(cos_mtx_cumi,use_dim = ncol(cos_mtx_cumi), n_component=2, n_neighbors = umap.neighbors)
        }
    }

    if(plot.diagnostics) {
        glist <- list()
        for(cur_barcode  in all_bcs) {
            cur_df = data.frame(
                raw_umi = bc_mtx[, cur_barcode],
                corrected_umi = cumi_mtx[, cur_barcode], # For now always set to TRUE so is computed
                pearson_residual = pr_mtx[,cur_barcode],
                cos_sim_raw = cos_mtx_raw[,cur_barcode],
                cos_sim_pr = cos_mtx_pr[,cur_barcode])
            cur_df$raw_total= rowSums(bc_mtx)
            cur_df$corrected_total = rowSums(cumi_mtx)
            cur_df <- cbind(cur_df, umap_res)
            show_pairs <- list(
                c("log(raw_total)", "log(raw_umi)", "cos_sim_raw"),
                c("log(raw_total)", "pearson_residual", "cos_sim_raw"),
                c("log(raw_umi)", "cos_sim_raw", "cos_sim_raw"),
                #c("log(corrected_total)", "log(corrected_umi)", "cos_sim_pr"),
                c("corrected_total", "corrected_umi", "cos_sim_pr"),
                c("log(corrected_total)", "pearson_residual", "cos_sim_pr"),
                c("log(corrected_umi)", "cos_sim_pr", "cos_sim_pr"),
                c("UMAP_1", "UMAP_2", "cos_sim_pr"),
                c("log(pearson_residual)", "cos_sim_pr", "cos_sim_pr")
            )

            plot_list <- list()
            for(i in 1:length(show_pairs)){
                pp = show_pairs[[i]]
                plot_list[[i]] = ggplot(cur_df, aes_string(pp[1], pp[2])) +
                    geom_point_rast(aes_string(color = pp[3]), stroke = 0, size = 1) +
                    scale_color_gradientn(colors = get_numeric_color("BlueGreenRed")) +
                    ggtitle(cur_barcode) +
                    theme_bw()
            }
            glist[[cur_barcode]] = arrangeGrob(grobs = plot_list, ncol = 3)
        }
        pdf(paste0(plot.path, "/", make.names(Sys.time()), "diagnostics.pdf"), width = 10, height =10)
        for(i in 1:length(glist)) {
            grid.newpage()
            grid.draw(glist[[i]])
        }
        dev.off()
    }

    return(
        list(
            res = assign_tbl,
            umap = umap_res,
            pr_mtx = pr_mtx,
            cumi_mtx = cumi_mtx,
            cos_mtx_pr = cos_mtx_pr,
            models = model_list
        )
    )
}

