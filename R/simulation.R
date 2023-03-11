


#' @export
simulateTags <- function(n.cell = 1000,
                         n.bc = 10,
                         seed = 1,
                         nb.theta = 10,
                         min.log.size = 1,
                         max.log.size = 7,
                         b0 = -3,
                         doublet.rate = .1,
                         separate.sim = TRUE) {
    set.seed(seed)
    bcs = paste0("bc", 1:n.bc)

    bc.contam.count = round(exp(runif(n.bc, min = 0, max = 5)))
    names(bc.contam.count) <- bcs

    # 1. Simulation of true initial staining
    cell.true.umi <- list()
    for(bc in bcs) {
        log.bc.umi = runif(n.cell, min = min.log.size, max = max.log.size)
        df = data.frame(ceiling(exp(log.bc.umi)))
        rownames(df) <- paste0(bc, "_", 1:nrow(df))
        colnames(df) = bc
        cell.true.umi[[bc]] <- df
    }
    cell.true.umi.mtx = do.call(plyr::rbind.fill, cell.true.umi)
    rownames(cell.true.umi.mtx) <- as.character(unlist(lapply(cell.true.umi, function(x) rownames(x))))
    cell.true.umi.mtx[is.na(cell.true.umi.mtx)] = 0

    if(separate.sim ) {
        # 2. Simulation of floating barcodes contaminating cell surface
        cell.surface.size= log(rowSums(cell.true.umi.mtx)) # Note here direct sum because 0 on other entries
        mu = exp(1*cell.surface.size + b0) # Assume same contamination level for all barcodes
        cell.contam.umi.mtx = sapply(bcs, function(bc) {
            rnegbin(length(cell.surface.size), mu = mu, theta = nb.theta)
        })

        # 3. Simulation of floating barcodes contaminating beads
        bead.contam.umi.mtx = sapply(bcs, function(bc) {
            rnegbin(n = length(cell.surface.size), mu = bc.contam.count[bc], theta = nb.theta)
        })

        final.umi.mtx = cell.true.umi.mtx + cell.contam.umi.mtx + bead.contam.umi.mtx
    }
    else {
        ## Combine 2 and 3 into a single NB? This can lead to a simpler model
        cell.surface.size= log(rowSums(cell.true.umi.mtx))
        mu_c = round(exp(1*cell.surface.size + b0))
        mu_list <- sapply(bcs, function(bc) {
            mu_c + bc.contam.count[bc]
        })
        cbn.contam.umi.mtx = sapply(bcs, function(bc) {
            rnegbin(length(cell.surface.size), mu = mu_list[,bc], theta = nb.theta)
        })
        final.umi.mtx = cell.true.umi.mtx + cbn.contam.umi.mtx
    }

    # Simulate doublet, inside cell.true.umi.mtx? TO BE DETERMINED
    if(doublet.rate > 0) {
        ndoublet = nrow(final.umi.mtx) * doublet.rate
        doub_idx1 = sample(1:nrow(final.umi.mtx), ndoublet)
        doub_idx2 = sample(1:nrow(final.umi.mtx), ndoublet)
        doub.umi.mtx = final.umi.mtx[doub_idx1, ] + final.umi.mtx[doub_idx2, ]
        rownames(doub.umi.mtx) <- paste0('doublet_', 1:nrow(doub.umi.mtx))
        final.umi.mtx <- rbind(final.umi.mtx, doub.umi.mtx)
    }

    return(final.umi.mtx)
}


