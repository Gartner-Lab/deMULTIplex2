


#' @export
simulateTags <- function(n.cell = 1000,
                         n.bc = 10,
                         seed = 1,
                         nb.theta = 10,
                         min.size.log = 1,
                         max.size.log = 7,
                         min.ambient.log = 0,
                         max.ambient.log = 5,
                         b0 = -3,
                         doublet.rate = .1,
                         separate.sim = TRUE,
                         cell.contam = T,
                         ambient.contam = T,
                         return.all = FALSE) {
    set.seed(seed)
    bcs = paste0("bc", 1:n.bc)

    bc.contam.count = round(exp(runif(n.bc, min = min.ambient.log, max = max.ambient.log)))
    names(bc.contam.count) <- bcs

    # 1. Simulation of true initial staining
    cell.true.umi <- list()
    for(bc in bcs) {
        log.bc.umi = runif(n.cell, min = min.size.log, max = max.size.log)
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

        if(cell.contam) {
            mu = exp(1*cell.surface.size + b0) # Assume same contamination level for all barcodes
            cell.contam.umi.mtx = sapply(bcs, function(bc) {
                rnegbin(length(cell.surface.size), mu = mu, theta = nb.theta)
            })
        } else {
            cell.contam.umi.mtx = cell.true.umi.mtx *0
        }

        # 3. Simulation of floating barcodes contaminating ambients
        if(ambient.contam) {
            ambient.contam.umi.mtx = sapply(bcs, function(bc) {
                rnegbin(n = length(cell.surface.size), mu = bc.contam.count[bc], theta = nb.theta)
            })
        } else {
            ambient.contam.umi.mtx = cell.true.umi.mtx *0
        }

        final.umi.mtx = cell.true.umi.mtx + cell.contam.umi.mtx + ambient.contam.umi.mtx
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
        doub.idx1 = sample(1:nrow(final.umi.mtx), ndoublet)
        doub.idx2 = sample(1:nrow(final.umi.mtx), ndoublet)
        # doub.umi.mtx = final.umi.mtx[doub.idx1, ] + final.umi.mtx[doub.idx2, ]
        ambient.contam.idx = sample(1:nrow(ambient.contam.umi.mtx), ndoublet)
        doub.umi.mtx = cell.true.umi.mtx[doub.idx1, ] + cell.contam.umi.mtx[doub.idx1, ] +
            cell.true.umi.mtx[doub.idx2, ] + cell.contam.umi.mtx[doub.idx2, ] + ambient.contam.umi.mtx[ambient.contam.idx, ]
        rownames(doub.umi.mtx) <- paste0('doublet_', 1:nrow(doub.umi.mtx), "|", rownames(final.umi.mtx)[doub.idx1], "_", rownames(final.umi.mtx)[doub.idx2])
        final.umi.mtx <- rbind(final.umi.mtx, doub.umi.mtx)
    }

    if(return.all) {
        if(separate.sim) {
            res <- list(
                cell.true.umi.mtx = cell.true.umi.mtx,
                cell.contam.umi.mtx = cell.contam.umi.mtx,
                ambient.contam.umi.mtx = ambient.contam.umi.mtx,
                final.umi.mtx = final.umi.mtx
            )
        } else {
            res <- list(
                cell.true.umi.mtx = cell.true.umi.mtx,
                cbn.contam.umi.mtx = cbn.contam.umi.mtx,
                final.umi.mtx = final.umi.mtx
            )
        }
    } else {
        res <- final.umi.mtx
    }
    return(res)
}


