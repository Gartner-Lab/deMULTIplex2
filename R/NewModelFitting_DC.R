wd <- "/Users/dannyconrad/Library/CloudStorage/Box-Box/Data/deMULTIplex2/Testing/EvePilot2/"
setwd(wd)

barTable <- readRDS("barTable.rds")

bc_mtx <- barTable[barTable$nUMI > 250, colnames(barTable) != "nUMI"]

# Generate canonical vector matrix for calculating cosine similarities
mtx_canon <- diag(x = 1, nrow = ncol(bc_mtx), ncol = ncol(bc_mtx))
colnames(mtx_canon) <- colnames(bc_mtx)
rownames(mtx_canon) <- colnames(bc_mtx)

# Calculate cosine similarity matrix
bc_mtx <- Matrix(as.matrix(bc_mtx), sparse = T)
vec_norm2 <- apply(bc_mtx, 1, function(x) norm(x, type = "2")) # for each cell, calculate its L2 norm (Euclidean distance from origin)
cos_mtx_raw <- bc_mtx / vec_norm2 # divide matrix by L2 norms to get cosine distances

### GLM-NB Fitting ###
barcodes <- colnames(bc_mtx)
model_list <- list()

# selection of these thresholds needs to be optimized (what about picking inflection points of sigmoid curve?)
posThresh <- 0.90
negThresh <- 1 - posThresh

glist <- list()

for (bc in 1:length(barcodes)) {
  
  cat("Fitting GLM-NB for ", bc, sep = "", fill = T)
  
  # iterate in case multiple thresholds are to be tested
  plist <- list()
  
  for (pos in 1:length(posThresh)) {
    
    df <- data.frame(Count = bc_mtx[, bc], 
                     nUMI = rowSums(bc_mtx),
                     CosSim = cos_mtx_raw[,bc])
    
    # pos_cells <- which(df$CosSim > max(df$CosSim) * (1-posThresh[pos]))
    # neg_cells <- which(df$CosSim < max(df$CosSim) * (1-posThresh[pos]))
    
    pos_cells <- which(df$CosSim > max(df$CosSim) * posThresh[pos])
    neg_cells <- which(df$CosSim < max(df$CosSim) * negThresh[pos])
    
    df$Initial <- NA
    df$Initial[pos_cells] <- "Pos"
    df$Initial[neg_cells] <- "Neg"
    
    ## Unsure how to best model the positive & negative populations:
    # model_pos <- glm.nb(Count ~ nUMI, data = df[pos_cells,], link = identity)
    # model_pos <- glm.nb(Count ~ offset(1*nUMI), data = df[pos_cells,], link = log)
    model_pos <- glm.nb(Count ~ nUMI, data = df[pos_cells,], link = identity)
    
    # model_neg <- glm.nb(Count ~ log(nUMI), data = df[neg_cells,], link = log)
    # model_neg <- glm.nb(Count ~ log(nUMI), data = df[neg_cells,], link = identity)
    model_neg <- glm.nb(Count ~ nUMI, data = df[neg_cells,], link = identity)
    
    # To plot regression lines just calculate coordinates of endpoints for geom_line() 
    lim <- data.frame(nUMI = c(min(df$nUMI),max(df$nUMI)))
    predict_pos <- predict(model_pos, lim, type = "response", se.fit=TRUE)
    predict_neg <- predict(model_neg, lim, type = "response", se.fit=TRUE)
    predict_pos <- cbind(lim, predict_pos)
    predict_neg <- cbind(lim, predict_neg)
    
    p1 <- ggplot(df, aes(x = nUMI, y = Count)) + geom_point(color = 'grey') +
      scale_x_log10() + scale_y_log10() +
      geom_line(data = predict_pos, aes(x = nUMI, y = fit), color = "blue", size = 1) +
      geom_line(data = predict_neg, aes(x = nUMI, y = fit), color = "red", size = 1) + 
      geom_abline(slope = 1, intercept = 0, color = "black", size = 1) + theme_classic()
    p2 <- ggplot(df, aes(x = nUMI, y = Count, color = Initial)) + geom_point() +
      scale_x_log10() + scale_y_log10() + theme_classic()
    p3 <- ggplot(df, aes(x = nUMI, y = Count, color = CosSim)) + geom_point() +
      scale_x_log10() + scale_y_log10() + theme_classic() +
      scale_color_gradientn(colors = get_numeric_color("BlueGreenRed"))
    
    plist[[pos]] <- plot_grid(p1,p2,p3, ncol = 3)
  }
  
  glist[[bc]] <- plot_grid(plotlist = plist, ncol = length(posThresh))
}

plot_grid(plotlist = glist, nrow = length(barcodes))
  
  

### Below is the code leftover from trying to run the GMM myself from that tutorial (this was previously inside the "bc in 1:length(barcodes)" loop)


  # Define Positive & Negative Classes
  df <- data.frame(Count = bc_mtx[, bc],
                   nUMI = rowSums(bc_mtx)) %>% as.matrix()
  df <- log10(df)
  df[!is.finite(df)] <- 0
  df <- as.data.frame(df)
  
  df$CosSim <- cos_mtx_raw[, bc]
  
  df$Initial <- "Pos"
  df$Initial[which(df$CosSim < max(df$CosSim) * (1-posThresh))] <- "Neg"
  
  ### Initialize parameters
  
  # Find mean of each variable for each separate class
  mu <- split(df[, c(1,2)], df$Initial)
  mu <- t(sapply(mu, colMeans))
  
  # Covariance Matrix for each initial class.
  cov <- replicate(diag(ncol(mu)), n = nclass, simplify = F)
  
  # Mixing Components
  a <- runif(nclass)
  a <- a/sum(a)
  
  # Visualize pos/neg classes and their means
  p <- ggplot(df, aes(x = nUMI, y = Count, color = Initial)) + 
    geom_point() + 
    # scale_y_log10() + scale_x_log10() + 
    geom_point(data = as.data.frame(mu), color = "black")
  
  ### Set up Loop to Iterate 100 times ###
  bic <- list()
  
  for (i in 1:10) {
    
    # Calculate PDF with class means and covariances.
    z <- sapply(1:nclass, function(x) {
      mvpdf(x = df[, c(1,2)], 
            mu = mu[x, ], 
            sigma = cov[[x]])
    })
    
    # Expectation Step for each class.
    r <- sapply(1:nclass, function(x) {
      (a[x] * z[, x])/rowSums(t((t(z) * a)))
    })
    
    # ggplot(df, aes(x = nUMI, y = Count, color = r[,1])) +
    #   geom_point() +
    #   geom_point(data = as.data.frame(mu), color = "black")
    # ggplot(df, aes(x = nUMI, y = Count, color = r[,2])) + 
    #   geom_point() + 
    #   geom_point(data = as.data.frame(mu), color = "black")
    
    # Choose the highest rowwise probability
    eK <- factor(apply(r, 1, which.max))
    
    table(eK) %>% print()
    
    ### Update parameters for next iteration
    
    # Total Responsibility
    mc <- colSums(r)
    
    # Update Mixing Components.
    a <- mc/nrow(df)
    
    # Update our Means
    mu <- sapply(1:nclass, function(x) {
      colSums(df[, 1:2] * r[, x]) * 1/mc[x]
    }) %>% t
    
    p <- p + geom_point(data = as.data.frame(mu), color = "black")
    
    # Update Covariance matrix.
    cov <- lapply(1:nclass, function(y) {
      t(r[, y] * t(apply(df[, 1:2], 1, function(x) x - mu[y, ]))) %*%
        (r[, y] * t(apply(df[, 1:2], 1, function(x) x - mu[y, ]))) * 1/mc[y]
    })
    
    ### Determine convergence with log-likelihood maximization
    
    # Compute the sum of the mixture densities, take the log, and add the column vector.
    loglik <- sum(log(apply(t(t(z) * a), 1, sum)))
    
    # Bayesian Information Criterion (BIC) is calculated using the equation
    # unsure how to determine the degrees of freedom k parameter (currently set at 10)
    bic[[i]] <- -2 * loglik + 10 * log(nrow(iris))
    
  }
  
  Final <- factor(eK, levels = 1:2, labels = c("Neg","Pos"))
  
  p2 <- ggplot(df, aes(x = nUMI, y = Count, color = Final)) + 
    geom_point() + geom_point(data = as.data.frame(mu), color = "black")
  p3 <- ggplot(df, aes(x = nUMI, y = Count, color = r[,2])) + 
    geom_point() + labs(color = "ProbPos")
  
  
  plot_grid(p, p2, p3, ncol = 3)
  
}



