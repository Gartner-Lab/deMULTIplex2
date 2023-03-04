

# Description:
# The new demultiplex function
#' @importFrom Matrix Matrix rowSums
#' @importFrom MASS glm.nb
#' @importFrom ggrastr geom_point_rast
#' @importFrom gridExtra arrangeGrob
#' @importFrom magrittr %>%
#' @export
classifyCells <- function(bc_mtx,
                          init.cos.cut = 0.9,
                          converge.threshold = 1e-3,
                          iter.max = 1e2,
                          prob.cut = 0.5,
                          residual.type = c("rqr", 'pearson'), # ONLY use RQR for future
                          plot.umap = c("residual", "umi"),
                          plot.diagnostics = F,
                          plot.path = getwd(),
                          plot.name = "",
                          umap.nn = 30L,
                          seed = 1,
                          point.size = 1,
                          label.size = 3,
                          min.bc.show = 50) {
    set.seed(seed)
    # TODO: Need to handle cells if rowSums = 0
    ## evaluate choices
    plot.umap <- match.arg(plot.umap)
    residual.type <- match.arg(residual.type)

    if(init.cos.cut < 0.5) {
        cat("Warning: setting init.cos.cut less than 0.5 is not recommended.", fill=T)
    }

    zero_bc_cells = rowSums(bc_mtx) == 0
    if (sum(zero_bc_cells) > 0) {
        message(paste0("Detected ", sum(zero_bc_cells), " cells with 0 barcode count. These cells will not be classified."))
        bc_mtx = bc_mtx[!zero_bc_cells,]
    }

    barcodes <- colnames(bc_mtx)
    # Calculate cosine similarity matrix
    bc_mtx <- Matrix(as.matrix(bc_mtx), sparse = T)
    vec_norm2 <- apply(bc_mtx, 1, function(x) norm(x, type = "2")) # for each cell, calculate its L2 norm (Euclidean distance from origin)
    cos.umi_mtx <- bc_mtx / vec_norm2 # divide matrix by L2 norms to get cosine distances

    pr_mtx <- matrix(nrow = nrow(bc_mtx), ncol = ncol(bc_mtx), data = NA)
    rownames(pr_mtx) <- rownames(bc_mtx)
    colnames(pr_mtx) <- colnames(bc_mtx)
    rqr_mtx <- pr_mtx
    prob_mtx <- pr_mtx

    coef_list <- list()

    for(bc in barcodes) {
        df <- data.frame(bc.umi = bc_mtx[, bc],
                         cos.umi = cos.umi_mtx[,bc])
        df$tt.umi <- rowSums(bc_mtx)

        cat("Running EM for ", bc, sep = "", fill = T)
        res <- fit.em(df,init.cos.cut = init.cos.cut, converge.threshold = converge.threshold, iter.max = iter.max)

        df <- res$df
        pr_mtx[,bc] <- df$pearson_residual
        rqr_mtx[,bc] <- df$rqr
        prob_mtx[,bc] <- df$p1
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
    vec_norm2 <- apply(res_mtx, 1, function(x) norm(x, type = "2"))
    cos.res_mtx <- res_mtx / vec_norm2

    # Assign barcodes to cells
    call_mtx <- sapply(colnames(prob_mtx), function(bc) {
        prob_mtx[,bc] > prob.cut
    })
    barcode_count = rowSums(call_mtx)
    barcode_assign <- apply(call_mtx, 1, function(x) {
        if (sum(x) == 1) colnames(call_mtx)[which(x)] else NA
    })
    assign_table = data.frame(barcode_assign = barcode_assign, barcode_count = barcode_count)

    # UMAP Plotting Options
    if (plot.umap == "umi") { # Raw UMI Counts
        umap_res <- compute_umap(bc_mtx, use_dim = ncol(bc_mtx), n_component=2, n_neighbors = umap.nn)
    } else if (plot.umap == "residual") { # Pearson Residuals
        umap_res <- compute_umap(res_mtx, use_dim = ncol(res_mtx), n_component=2, n_neighbors = umap.nn)
    } else {
        umap_res <-  NA
    }

    glist <- list()

    if(1) { # Always plot
        message("Plotting final umap...")
        umap_df <- cbind(umap_res, assign_table)
        umap_df$barcode_count = as.character(umap_df$barcode_count)
        umap_df$barcode_count[umap_df$barcode_count >= 3] = ">=3"
        umap_df$barcode_count <- factor(umap_df$barcode_count, c("0", "1", "2" ,">=3"))

        glist[["summary"]] <- plot_summary(umap_df, point.size = point.size, label.size = label.size, min.bc.show = min.bc.show)
    }

    # DEBUG
    # if(1) {
    #     bc == barcodes[1]
    #     df <- data.frame(
    #         bc.umi = bc_mtx[, bc],
    #         res = res_mtx[,bc],
    #         cos.umi = cos.umi_mtx[,bc],
    #         cos.res = cos.res_mtx[,bc])
    #     df$tt.umi <- rowSums(bc_mtx)
    #     df <- cbind(df, umap_df)
    #     df$is_negcell = 1:nrow(df) %in% neg_cell_list[[bc]]
    #     assign("df",df, env = .GlobalEnv)
    # }

    # Diagnostic Plots
    if (plot.diagnostics) {
        for (bc in barcodes) {
            df <- data.frame(
                bc.umi = bc_mtx[, bc],
                res = res_mtx[,bc],
                cos.umi = cos.umi_mtx[,bc],
                cos.res = cos.res_mtx[,bc],
                prob.pos = prob_mtx[,bc])
            df$tt.umi <- rowSums(bc_mtx)
            df <- cbind(df, umap_df)
            mappings <- list(
                c("log(tt.umi)", "log(bc.umi)", "cos.umi"),
                c("log(tt.umi)", "log(bc.umi)", "prob.pos"),
                c("log(tt.umi)", "log(tt.umi-bc.umi)", "prob.pos"),
                #c("log(bc.umi)", "cos.umi", "cos.umi"),
                c("log(bc.umi)", "cos.umi", "prob.pos"),
                c("log(tt.umi)", "res", "prob.pos"),
                #c("log(tt.umi)", "cos.umi", "cos.umi"),
                #c("log(tt.umi)", "res", "cos.umi"),
                #c("qq"),
                #c("log(res)", "res", "cos.res"),
                #c("log(res)", "cos.res", "cos.res"),
                #c("UMAP_1", "UMAP_2", "res"),
                c("UMAP_1", "UMAP_2", "prob.pos"),
                c("UMAP_1", "UMAP_2", "barcode_count")
            )
            assign("df", df, env = .GlobalEnv)
            df$barcode_count[call_mtx[,bc] == 0] = NA
            glist[[bc]] <- plot.all.diagnostics(df, mappings, bc = bc, point.size = point.size, ncol=3)
        }
    }

    time <- (Sys.time() %>% make.names() %>% strsplit(split = "X") %>% unlist())[2]
    pdf(paste0(plot.path, "/", plot.name, "_", time, "_assignment.pdf"), width = 20, height = 15)
    for(i in 1:length(glist)) {
        grid.newpage()
        message(paste0("Plotting ", names(glist)[i]))
        grid.draw(glist[[i]])
    }
    dev.off()

    return(
        list(
            assign_table = assign_table,
            umap = umap_res %>% as.matrix(),
            res_mtx = res_mtx %>% as.matrix(),
            rqr_mtx = rqr_mtx,
            pr_mtx = pr_mtx,
            prob_mtx = prob_mtx,
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
    df_pred$prob_neg <- dnbinom(df_pred$bc.umi, mu=df_pred$fit, size=model$theta)

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
    df_pred$prob_neg = dpois(df_pred$bc.umi, lambda = df_pred$fit)

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


plot_summary <- function(df, point.size = 1, label.size = 3 , min.bc.show = 50) {
    unq_bcs = unique(df$barcode_assign)
    unq_bcs = unq_bcs[!is.na(unq_bcs)]
    use_color = get_factor_color(unq_bcs, "Set1")
    names(use_color) = unq_bcs

    g1 <- ggplot(df, aes_string("UMAP_1", "UMAP_2")) +
        geom_point_rast(aes_string(color = "barcode_assign"), stroke = 0, size = point.size) +
        #geom_point_rast(data = df[!is.na(df[["barcode_assign"]]), ], aes_string(color = "barcode_assign"), stroke = 0, size = .5) +
        scale_color_manual(values = use_color, na.value='lightgrey') +
        theme_bw() +
        ggtitle("barcode_assign") +
        guides(color = "none")
    show_bc = names(which(table(df[["barcode_assign"]]) > min.bc.show))
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
        geom_point_rast(aes_string(color = "barcode_count"), stroke = 0, size = point.size) +
        scale_color_manual(values = count_color, na.value='lightgrey') +
        theme_bw() +
        ggtitle("barcode_count")
    return(arrangeGrob(g1,g2, ncol = 2))
}


plot.all.diagnostics <- function(df, mappings, bc, point.size = 1, ncol = 3) {
    plot_list <- list()
    for(i in 1:length(mappings)){
        map <- mappings[[i]]
        if(map == "qq") {
            # QQ plot
            p = ggplot(df, aes(sample=res))+stat_qq(size = point.size, stroke = 0) +
                stat_qq_line() +
                xlab("Normal quantiles") +
                ylab("Residual quantiles") +
                theme_bw()
        } else if(map == "qq_init") {
            p = ggplot(dd1, aes(x = theoretical, y = sample)) + geom_point(alpha = 0.3) +
                geom_abline(intercept = 0, slope = 1, color = "blue") +
                theme_classic() +
                ggtitle("QQ plot of standardized residuals")
        } else if (map == "qq_final") {
            p = ggplot(dd2, aes(x = theoretical, y = sample)) + geom_point(alpha = 0.3) +
                geom_abline(intercept = 0, slope = 1, color = "blue") +
                theme_classic() +
                ggtitle("QQ plot of standardized residuals")
        }else {
            p <- ggplot(df, aes_string(map[1], map[2])) +
                geom_point_rast(aes_string(color = map[3]), stroke = 0, size = point.size) +
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



fit.em <- function(df, init.cos.cut = .9, converge.threshold = 1e-3, iter.max = 1e2) {
    # Initialization
    mem.iter <- as.numeric(df$cos.umi > init.cos.cut)
    pi0 <- sum(mem.iter==0)/length(mem.iter)
    pi1 <- sum(mem.iter==1)/length(mem.iter)

    fail.fit.I = 0

    if(sum(mem.iter==1) < 10) {
        cat("Less than 10 poistive cells detected in initialization. Consider changing the initial cosine cutoff to lower, and check if the well contain stained cells.", fill=T)
        fail.fit.I = 1
    } else if(sum(mem.iter==0) < 10) {
        cat("Less than 10 negative cells detected in initialization. Consider changing the initial cosine cutoff to higher.", fill=T)
        fail.fit.I = 1
    } else {
        tryCatch({
            # Fit 0 part with glm.NB
            fit0 = glm.nb("bc.umi~log(tt.umi)", df[mem.iter ==0, ], link = log)
            fit1 = glm.nb("(tt.umi - bc.umi)~log(tt.umi)", df[mem.iter ==1, ], link = log)
        }, error = function(e) {
            cat("Fitting failed.", fill=T)
            fail.fit.I <<- 1
        })
    }

    glist <- list()

    if(!fail.fit.I) {
        pred0 <- predict(fit0, df, type = "response", se.fit=FALSE)
        prob0 <- dnbinom(df$bc.umi, mu=pred0, size=fit0$theta)
        pred1 <- predict(fit1, df, type = "response", se.fit=FALSE)
        prob1 <- dnbinom(df$tt.umi - df$bc.umi, mu=pred1, size=fit1$theta)

        prob0[df$bc.umi < pred0] = 1

        comp0 <- pi0 * prob0
        comp1 <- pi1 * prob1
        comp.sum <- comp0 + comp1

        df$pred0 <- pred0
        df$prob0 <- prob0
        df$pred1 <- pred1
        df$prob1 <- prob1
        df$p0 <- comp0/comp.sum
        df$p1 <- comp1/comp.sum

        glist[[1]] <- plot.em.diagnostics(df)
        sum.finite <- function(x) {
            sum(x[is.finite(x)])
        }

        Q <- 0
        Q[2] <- sum.finite(log(pi0)+log(df$prob0)) + sum.finite(log(pi1)+log(df$prob1))

        k <- 2

        while (abs(Q[k]-Q[k-1])>=converge.threshold & !fail.fit.I & k <= iter.max) {
            print(abs(Q[k]-Q[k-1]))
            # E step
            comp0 <- pi0 * df$prob0
            comp1 <- pi1 * df$prob1
            comp.sum <- comp0 + comp1

            p0 <- comp0/comp.sum
            p1 <- comp1/comp.sum

            # M step
            pi0 <- sum.finite(p0) / nrow(df)
            pi1 <- sum.finite(p1) / nrow(df)

            mem.iter = as.numeric(p1 > 0.5)

            if(sum(mem.iter==1) < 10) {
                cat("Less than 10 poistive cells detected in initialization. Consider changing the initial cosine cutoff to lower, and check if the well contain stained cells.", fill=T)
                fail.fit.I = 1
            } else if(sum(mem.iter==0) < 10) {
                cat("Less than 10 negative cells detected in initialization. Consider changing the initial cosine cutoff to higher.", fill=T)
                fail.fit.I = 1
            } else {
                tryCatch({
                    fit0 = glm.nb("bc.umi~log(tt.umi)", df[mem.iter ==0, ], link = log)
                    fit1 = glm.nb("(tt.umi - bc.umi)~log(tt.umi)", df[mem.iter ==1, ], link = log)
                }, error = function(e) {
                    cat("Fitting failed.", fill=T)
                    fail.fit.I <<- 1
                })
            }

            if(fail.fit.I) break
            pred0 <- predict(fit0, df, type = "response", se.fit=FALSE)
            prob0 <- dnbinom(df$bc.umi, mu=pred0, size=fit0$theta)
            pred1 <- predict(fit1, df, type = "response", se.fit=FALSE)
            prob1 <- dnbinom(df$tt.umi - df$bc.umi, mu=pred1, size=fit1$theta)

            prob0[df$bc.umi < pred0] = 1

            k <- k + 1
            Q[k] <- sum.finite(log(pi0)+log(prob0)) + sum.finite(log(pi1)+log(prob1))

            # Plotting
            df$prob0 <- prob0
            df$pred0 <- pred0
            df$p0 <- p0
            df$prob1 <- prob1
            df$pred1 <- pred1
            df$p1 <- p1
            glist[[k-1]] <-plot.em.diagnostics(df)
        }
    }

    if(fail.fit.I){
        # Refine later
        df$prob0 <- 1
        df$pred0 <- NA
        df$p0 <- 1
        df$prob1 <- 0
        df$pred1 <- NA
        df$p1 <- 0
        fit0 <- list(coefficients = NA, theta = NA)
        fit1 <- list(coefficients = NA, theta = NA)
    }

    # Compute some resduals
    #assign("df", df, env = .GlobalEnv)
    df$pearson_residual <- (df[['bc.umi']]-df$pred0)/sqrt(df$pred0 + df$pred0^2/fit0$theta)
    df$rqr <- rqr.nb(df, y="bc.umi",fit="pred0", model = fit0)

    return(list(
        df = df,
        plots = glist,
        fit0 = c(fit0$coefficients, theta = fit0$theta),
        fit1 = c(fit1$coefficients, theta = fit1$theta)
    ))
}


plot.em.diagnostics <- function(df) {
    g1 <- ggplot(df, aes_string(color = "prob0")) +
        geom_point_rast(aes_string("log(tt.umi)", "log(bc.umi)"), size = 1, stroke = 0) +
        geom_point_rast(aes_string("log(tt.umi)", "log(pred0)"), color = "grey", size = 1, stroke = 0) +
        scale_color_gradientn(colors =get_gradient_color("BlueGreenRed"))

    g2 <- ggplot(df, aes_string(color = "prob1")) +
        geom_point_rast(aes_string("log(tt.umi)", "log(tt.umi-bc.umi)"), size = 1, stroke = 0) +
        geom_point_rast(aes_string("log(tt.umi)", "log(pred1)"), color = "grey", size = 1, stroke = 0) +
        scale_color_gradientn(colors =get_gradient_color("BlueGreenRed"))

    g3 <- ggplot(df, aes_string(color = "p0")) +
        geom_point_rast(aes_string("log(tt.umi)", "log(bc.umi)"), size = 1, stroke = 0) +
        geom_point_rast(aes_string("log(tt.umi)", "log(pred0)"), color = "grey", size = 1, stroke = 0) +
        scale_color_gradientn(colors =get_gradient_color("BlueGreenRed"))

    g4 <- ggplot(df, aes_string(color = "p1")) +
        geom_point_rast(aes_string("log(tt.umi)", "log(tt.umi-bc.umi)"), size = 1, stroke = 0) +
        geom_point_rast(aes_string("log(tt.umi)", "log(pred1)"), color = "grey", size = 1, stroke = 0) +
        scale_color_gradientn(colors =get_gradient_color("BlueGreenRed"))

    return(arrangeGrob(grobs = list(g1,g2,g3,g4), ncol = 2))
}




