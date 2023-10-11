
# Simulate tag count matrix for benchmarking
#' @importFrom MASS rnegbin
#' @export
simulateTags <- function(n.cell = 1000, # If supplied as a vector, then has to be same length as barcode, apply to each of the barcode
                         n.bc = 10,
                         seed = 1,
                         nb.theta = 10,
                         cell.staining.meanlog = 6, # Initial average log level of cell staining for each tag (before pooling), if supplied as vector must equal to number of tags
                         cell.staining.sdlog = 1, # Variance in initial staining (log scale), assume same sdlog across tags for simplicity
                         ambient.meanlog = 4, # Ambient average log level of each tag (before pooling), if supplied as vector must equal to number of tags
                         b0 = -3, # Controls surface bound contamination rate, exp(-3) ~ 5% tag reads are surface bound contamination
                         doublet.rate = .1,
                         dropout.lambda = Inf,
                         #separate.sim = TRUE,
                         cell.contam = T,
                         ambient.contam = T,
                         return.all = FALSE) {
    set.seed(seed)
    bcs = paste0("bc", 1:n.bc)

    # 1. Simulation of true initial staining
    cell.true.umi <- list()
    if(length(n.cell) > 1){ # If supplied as a vector, then has to be same length as barcode, apply to each of the barcode
        if(length(n.cell) != n.bc) stop("If n.cell is a vector, then its length must be equal to the number of tags.")
        n.cell.bc <- n.cell
    } else {
        n.cell.bc <- rep(n.cell, n.bc)
    }
    names(n.cell.bc) <- bcs

    if(length(cell.staining.meanlog) > 1) {
        if(length(cell.staining.meanlog) != n.bc) stop("If cell.staining.meanlog is a vector, then its length must be equal to the number of tags.")
        cell.staining.meanlog.bc <- cell.staining.meanlog
    } else {
        cell.staining.meanlog.bc <- rep(cell.staining.meanlog, n.bc)
    }
    names(cell.staining.meanlog.bc) <- bcs

    for(bc in bcs) {
        df = data.frame(ceiling(rlnorm(n.cell.bc[bc], meanlog = cell.staining.meanlog.bc[bc], sdlog = cell.staining.sdlog)))
        rownames(df) <- paste0(bc, "_", 1:nrow(df))
        colnames(df) = bc
        cell.true.umi[[bc]] <- df
    }
    cell.true.umi.mtx = do.call(plyr::rbind.fill, cell.true.umi)
    rownames(cell.true.umi.mtx) <- as.character(unlist(lapply(cell.true.umi, function(x) rownames(x))))
    cell.true.umi.mtx[is.na(cell.true.umi.mtx)] = 0

    # 2. Simulation of cell surface contamination
    cell.surface.size= log(rowSums(cell.true.umi.mtx)) # Note here direct sum because 0 on other entries

    if(cell.contam) {
        mu = exp(1*cell.surface.size + b0)
        cell.contam.umi.mtx = sapply(bcs, function(bc) {
            rnegbin(length(cell.surface.size), mu = mu, theta = nb.theta)
        })
    } else {
        cell.contam.umi.mtx = cell.true.umi.mtx *0
    }

    # 3. Simulation of floating barcodes contaminating ambient
    if(ambient.contam) {
        if(length(ambient.meanlog) > 1) {
            if(length(ambient.meanlog) != n.bc) stop("If ambient.meanlog is a vector, then its length must be equal to the number of tags.")
            ambient.meanlog.bc <- ambient.meanlog
        } else {
            ambient.meanlog.bc <- rep(ambient.meanlog, n.bc)
        }
        names(ambient.meanlog.bc) <- bcs
        if(any(ambient.meanlog.bc > cell.staining.meanlog.bc)) stop("Ambient contamination level cannot be greater than initial staining level.")
        ambient.contam.umi.mtx = sapply(bcs, function(bc) {
            rnegbin(n = length(cell.surface.size), mu = exp(ambient.meanlog.bc[bc]), theta = nb.theta)
        })
    } else {
        ambient.contam.umi.mtx = cell.true.umi.mtx *0
    }

    actual.umi.mtx = cell.true.umi.mtx + cell.contam.umi.mtx + ambient.contam.umi.mtx

    # Simulate doublet
    if(doublet.rate > 0) {
        ndoublet = nrow(actual.umi.mtx) * doublet.rate
        doub.idx1 = sample(1:nrow(actual.umi.mtx), ndoublet)
        doub.idx2 = sample(1:nrow(actual.umi.mtx), ndoublet)
        # doub.umi.mtx = actual.umi.mtx[doub.idx1, ] + actual.umi.mtx[doub.idx2, ]
        ambient.contam.idx = sample(1:nrow(ambient.contam.umi.mtx), ndoublet)
        doub.umi.mtx = cell.true.umi.mtx[doub.idx1, ] + cell.contam.umi.mtx[doub.idx1, ] +
            cell.true.umi.mtx[doub.idx2, ] + cell.contam.umi.mtx[doub.idx2, ] + ambient.contam.umi.mtx[ambient.contam.idx, ]
        rownames(doub.umi.mtx) <- paste0('doublet_', 1:nrow(doub.umi.mtx), "|", rownames(actual.umi.mtx)[doub.idx1], "_", rownames(actual.umi.mtx)[doub.idx2])
        actual.umi.mtx <- rbind(actual.umi.mtx, doub.umi.mtx)
    }

    # Simulate drop out events
    if(!is.infinite(dropout.lambda)) {
        # Compute drop out probablity based on actual umi count
        dropout.prob.mtx <- as.matrix(exp(-dropout.lambda * actual.umi.mtx^2))
        # Determine for each barcode and each cell if the value should be dropped out (create indicator)
        dropout.indicator <- matrix(rbinom(length(dropout.prob.mtx),prob=dropout.prob.mtx,size=1),nrow=nrow(dropout.prob.mtx))
        final.umi.mtx <- actual.umi.mtx * (1-dropout.indicator)
    } else {
        dropout.prob.mtx <- matrix(0,nrow(actual.umi.mtx), ncol(actual.umi.mtx))
        dropout.indicator <- dropout.prob.mtx
        final.umi.mtx <- actual.umi.mtx
    }



    if(return.all) {
        res <- list(
            input.parameters = list(
                n.cell = n.cell, # If supplied as a vector, then has to be same length as barcode, apply to each of the barcode
                n.bc = n.bc,
                seed = seed,
                nb.theta = nb.theta,
                cell.staining.meanlog = cell.staining.meanlog, # Initial average log level of cell staining for each tag (before pooling), if supplied as vector must equal to number of tags
                cell.staining.sdlog = cell.staining.sdlog, # Variance in initial staining (log scale), assume same sdlog across tags for simplicity
                ambient.meanlog = ambient.meanlog, # Ambient average log level of each tag (before pooling), if supplied as vector must equal to number of tags
                b0 = b0, # Controls surface bound contamination rate, exp(-3) ~ 5% tag reads are surface bound contamination
                doublet.rate = doublet.rate,
                dropout.lambda = dropout.lambda,
                cell.contam =cell.contam,
                ambient.contam = ambient.contam
            ),

            cell.true.umi.mtx = cell.true.umi.mtx,
            cell.contam.umi.mtx = cell.contam.umi.mtx,
            ambient.contam.umi.mtx = ambient.contam.umi.mtx,
            actual.umi.mtx = actual.umi.mtx,
            dropout.prob.mtx = dropout.prob.mtx,
            dropout.indicator = dropout.indicator,
            final.umi.mtx = final.umi.mtx
        )
    } else {
        res <- final.umi.mtx
    }
    return(res)
}



#' @export
downsample_mtx_umi <- function(mtx, ratio = .1, seed = 1) {
    set.seed(seed)
    ds_mtx <- t(apply(mtx, 1, function(x) {
        n <- floor(sum(x) * ratio)
        ds_reads <- sort(sample(seq_len(sum(x)), n))
        read_breaks <- c(0, cumsum(x))
        hist(ds_reads, breaks = read_breaks, plot = FALSE)$count
    }))
    colnames(ds_mtx) <- colnames(mtx)
    return(ds_mtx)
}

