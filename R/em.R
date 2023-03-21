

# Code for EM, used by demutiplexTags in classify.R

sum.finite <- function(x) {
    sum(x[is.finite(x)])
}

fit.em <- function(df, init.cos.cut = .9, converge.threshold = 1e-3, max.iter = 1e2, min.cell.fit = 10, max.cell.fit = 1e4, min.quantile.fit = .05, max.quantile.fit = .95, plot = FALSE) {

    # Initialization
    mem.init = as.numeric(df$cos.umi > init.cos.cut)
    m.res <- m.step(df, posterior.prob = NULL, mem.init = mem.init, min.cell.fit = min.cell.fit, max.cell.fit = max.cell.fit, min.quantile.fit = min.quantile.fit, max.quantile.fit = max.quantile.fit)

    glist <- list()

    if(!m.res$fail.fit.I) {
        e.res <- e.step(df, m.res$fit0, m.res$fit1, m.res$pi.vector, plot = plot)
        glist[[1]] <- e.res$plot

        Q <- 0
        k <- 2
        Q[k] <- e.res$loglik

        message(paste0("Iteration ", k-1,"; Q diff: ",abs(Q[k]-Q[k-1])))


        while (abs(Q[k]-Q[k-1])>=converge.threshold & !m.res$fail.fit.I & k <= max.iter) {
            # M step
            m.res <- m.step(df, posterior.prob = e.res$posterior.prob, min.cell.fit = min.cell.fit, max.cell.fit = max.cell.fit, min.quantile.fit = min.quantile.fit, max.quantile.fit = max.quantile.fit)
            if(m.res$fail.fit.I) break

            # E step
            e.res <- e.step(df, m.res$fit0, m.res$fit1, m.res$pi.vector, plot = plot)
            glist[[k-1]] <- e.res$plot
            k <- k + 1
            Q[k] <- e.res$loglik
            Q_diff <- abs(Q[k]-Q[k-1])
            message(paste0("Iteration ", k-1,"; Q diff: ",Q_diff))
            if(is.infinite(Q_diff)) m.res$fail.fit.I = 1
        }

        if(k > max.iter && abs(Q[k]-Q[k-1]) >= converge.threshold) {
            message("Max number of iteration reached.")
            #cat("Warning: max number of iteration reached but algorithm did not converge to specified range. This is normally ok but please check the diagnostic plots. Consider increasing max.iter or max.cell.fit if fit is not good.", fill=T)
        }
    }

    if(m.res$fail.fit.I){
        # Refine later
        fit0 <- list(coefficients = NA, theta = NA)
        fit1 <- list(coefficients = NA, theta = NA)
        df$prob0 <- 1
        df$pred0 <- NA
        df$post0 <- 1
        df$prob1 <- 0
        df$pred1 <- NA
        df$post1 <- 0
    } else {
        fit0 <- m.res$fit0
        fit1 <- m.res$fit1
        df$pred0 <- predict(fit0, df, type = "response", se.fit=FALSE)
        df$prob0 <- dnbinom(df$bc.umi, mu=df$pred0, size=fit0$theta)
        df$pred1 <- predict(fit1, df, type = "response", se.fit=FALSE)
        df$prob1 <- dnbinom(df$tt.umi - df$bc.umi, mu=df$pred1, size=fit1$theta)
        #assign("df1",df, env=.GlobalEnv)
        #assign("fit0", fit0, env=.GlobalEnv)
        # Which is best to use? more robust?
        #df$prob0[df$bc.umi < df$pred0] = 1
        df$prob0[df$bc.umi < df$pred0] = dnbinom(ceiling(df$pred0[df$bc.umi < df$pred0]), mu=df$pred0[df$bc.umi < df$pred0], size=fit0$theta)
        df$post0 <- e.res$posterior.prob[,1]
        df$post1 <- e.res$posterior.prob[,2]
    }

    # Compute some resduals
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

    g3 <- ggplot(df, aes_string(color = "post0")) +
        geom_point_rast(aes_string("log(tt.umi)", "log(bc.umi)"), size = 1, stroke = 0) +
        geom_point_rast(aes_string("log(tt.umi)", "log(pred0)"), color = "grey", size = 1, stroke = 0) +
        scale_color_gradientn(colors =get_gradient_color("BlueGreenRed"))

    g4 <- ggplot(df, aes_string(color = "post1")) +
        geom_point_rast(aes_string("log(tt.umi)", "log(tt.umi-bc.umi)"), size = 1, stroke = 0) +
        geom_point_rast(aes_string("log(tt.umi)", "log(pred1)"), color = "grey", size = 1, stroke = 0) +
        scale_color_gradientn(colors =get_gradient_color("BlueGreenRed"))

    return(arrangeGrob(grobs = list(g1,g2,g3,g4), ncol = 2))
}



e.step <- function(df, fit0, fit1, pi.vector, plot = FALSE) {
    pred0 <- predict(fit0, df, type = "response", se.fit=FALSE)
    prob0 <- dnbinom(df$bc.umi, mu=pred0, size=fit0$theta)
    pred1 <- predict(fit1, df, type = "response", se.fit=FALSE)
    prob1 <- dnbinom(df$tt.umi - df$bc.umi, mu=pred1, size=fit1$theta)
    prob0[df$bc.umi < pred0] = 1

    pi0 = pi.vector[1]
    pi1 = pi.vector[2]
    prob0 = prob0
    prob1 = prob1

    comp0 <- pi0 * prob0
    comp1 <- pi1 * prob1
    comp.sum <- comp0 + comp1

    post0 <- comp0/comp.sum
    post1 <- comp1/comp.sum

    comp.sum.ln <- log(comp.sum, base = exp(1))
    comp.sum.ln.sum <- sum(comp.sum.ln)


    # Plotting
    df$pred0 <- pred0
    df$prob0 <- prob0
    df$pred1 <- pred1
    df$prob1 <- prob1
    df$post0 <- post0
    df$post1 <- post1
    if(plot) {
        g1 <- plot.em.diagnostics(df)
    } else {
        g1 <- NULL
    }

    list("loglik" = comp.sum.ln.sum,
         "posterior.prob" = cbind(post0, post1),
         "plot" = g1)
}


m.step <- function(df, posterior.prob, mem.init = NULL, min.cell.fit = 10, max.cell.fit = 1e4, min.quantile.fit = .05, max.quantile.fit = .95) {
    if(!is.null(mem.init)) {
        df$mem.iter = mem.init
        pi0 <- sum.finite(df$mem.iter == 0) / length(df$mem.iter)
        pi1 <- sum.finite(df$mem.iter == 1) / length(df$mem.iter)
    } else {
        post0 <- posterior.prob[,1]
        post1 <- posterior.prob[,2]
        df$mem.iter = as.numeric(post1 > 0.5)
        pi0 <- sum.finite(post0) / length(df$mem.iter)
        pi1 <- sum.finite(post1) / length(df$mem.iter)
    }

    fail.fit.I = 0
    fit0 <- NULL
    fit1 <- NULL

    if(any(is.na(df$mem.iter))) {
        cat(paste0("Something went wrong in fitting"), fill=T)
        fail.fit.I = 1
    } else if(sum(df$mem.iter==1) < min.cell.fit) {
        cat(paste0("Less than ", min.cell.fit,
                   " poistive cells detected in initialization. Consider decreasing the init.cos.cut, and check if the well contain stained cells."), fill=T)
        fail.fit.I = 1
    } else if(sum(df$mem.iter==0) < min.cell.fit) {
        cat("Less than ", min.cell.fit,
            " negative cells detected in initialization. Consider increasing the init.cos.cut.", fill=T)
        fail.fit.I = 1
    } else {
        min.tt.umi = quantile(df$tt.umi, min.quantile.fit)
        max.tt.umi = quantile(df$tt.umi, max.quantile.fit)
        df <- df[df$tt.umi >= min.tt.umi & df$tt.umi <= max.tt.umi, ,drop=F]
        if(sum(df$mem.iter==0) > max.cell.fit) {
            df_fit0 = df[sample(which(df$mem.iter ==0), max.cell.fit), ]
        } else {
            df_fit0 = df[df$mem.iter == 0, ]
        }
        if(sum(df$mem.iter==1) > max.cell.fit) {
            df_fit1 = df[sample(which(df$mem.iter ==1), max.cell.fit), ]
        } else {
            df_fit1 = df[df$mem.iter ==1, ]
        }
        tryCatch({
            fit0 = glm.nb("bc.umi~log(tt.umi)", df_fit0, link = log)
            fit1 = glm.nb("(tt.umi - bc.umi)~log(tt.umi)", df_fit1, link = log)
        }, error = function(e) {
            cat("Fitting failed.", fill=T)
            fail.fit.I <<- 1
        })
    }

    list("fit0" = fit0,
         "fit1" = fit1,
         "pi.vector" = c(pi0,pi1),
         "fail.fit.I" = fail.fit.I)
}



