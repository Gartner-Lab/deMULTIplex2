

#' @export
tags <- readxl::read_excel("data-raw/MULTI-seq_oligos_Nov2019.xlsx")$`Barcode Sequence`

#' @export
'%ni%' <- Negate('%in%')

### Some helper functions ###
#' @export
factor_color_opt <- function() {
    allowed_pals <- c("Set1", "Set2", "Paired", "Dark2", "Accent")
    return(allowed_pals)
}

#' @export
get_factor_color <-function (labels, pal = "Set1", maxCol = 9, nogrey = T)
{
    unq <- unique(labels)
    hmcol <- RColorBrewer::brewer.pal(maxCol, pal)
    if(nogrey) {
        hmcol <- hmcol[!hmcol %in% c("#999999","#B3B3B3")]
    }
    colv <- rep(NA, length(labels))
    #if (length(unq) > maxCol) {
    cp <- colorRampPalette(hmcol)
    hmcol <- cp(length(unq))
    #}
    for (i in 1:length(unq)) {
        colv[labels == unq[i]] <- hmcol[i]
    }
    return(colv)
}


#' @export
gg_color_hue2 <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}

#' @export
floor_dec <- function(x, level=1) round(x - 5*10^(-level-1), level)

#' @export
ceiling_dec <- function(x, level=1) round(x + 5*10^(-level-1), level)



#' @export
numeric_color_opt <- function() {
    allowed_pals <- c("BlueGreenRed", 'RdYlBu', 'RdBu','RdOgYl', 'viridis', 'magma', 'plasma', 'inferno', 'rainbow', 'gg_color_hue', 'grey&red')
    return(allowed_pals)
}

#' @export
get_gradient_color <- function(palette = NULL, cnum = 100) {
    if(is.null(palette)) stop("please specify a palette")
    allowed_pals <- numeric_color_opt()

    if(!palette %in% allowed_pals) stop(paste0("Please specify one of '", paste(allowed_pals, collapse = "', '"), "'."))
    if(palette %in% get_brewer_set(c("sequential", "diverging"))) {
        colorRampPalette(rev(RColorBrewer::brewer.pal(9,palette)))(cnum)
    } else if(palette %in% list("viridis" = "viridis", "magma" = "magma", "plasma" = "plasma", "inferno" = "inferno")) {
        viridis(n = cnum, option = palette)
    } else if(palette == "diverge_hcl") {
        colorRampPalette(colorspace::diverge_hcl(7))(cnum)
    } else if(palette == "redgreen") {
        rev(gplots::redgreen(75))
    } else if(palette == "rainbow") {
        colorRampPalette(rev(rainbow(10)))(cnum)
    } else if(palette == "grey&red") {
        colorRampPalette(c("grey", "#b2182b"))(cnum)
    } else if(palette == "RdOgYl") {
        colorRampPalette(c("grey85", "red", "orange", "yellow"))(cnum)
    } else if(palette == "gg_color_hue") {
        gg_color_hue2(cnum)
    } else if(palette == "blue_green_gold"){
        colorRampPalette(c("grey85", "blue", "green", "#FFD200", "gold"))(cnum)
    } else if(palette == "black_red_gold"){
        colorRampPalette(c("grey85", "black", "red", "#FFD200"))(cnum)
    } else if(palette == "black_red") {
        colorRampPalette(c("grey85", "black", "red"))(cnum)
    } else if(palette == "red_yellow") {
        colorRampPalette(c("grey85",  "red", "yellow"))(cnum)
    } else if(palette == "black_yellow") {
        colorRampPalette(c("grey85",  "black", "yellow"))(cnum)
    } else if(palette == "black_yellow_gold") {
        colorRampPalette(c("grey85",  "black", "yellow", "gold"))(cnum)
    } else if(palette == "BlueGreenRed") {
        colorRampPalette(c("midnightblue", "dodgerblue", "seagreen", "#00C000", "gold2", "darkorange1", "red1"))(cnum)
    }
}

#' @export
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}


#' @export
numeric_rgb_range <- function(col = NULL, zgrey=F) {
    rgb_scale <- apply(col2rgb(col),2,function(x){paste0(x,collapse =",")})
    rgb_scale_text<-paste0("rgb(", rgb_scale, ")")
    if(zgrey) {
        cus_range <- seq(1e-10, 1, length.out =length(rgb_scale))
        rgb_scale_list<- list(c(0, "rgb(178, 178, 178)"))
    } else {
        cus_range <- seq(0, 1, length.out =length(rgb_scale))
        rgb_scale_list<- list()
    }
    for(i in 1:length(rgb_scale)) {
        rgb_scale_list[[length(rgb_scale_list) + 1]] <- c(cus_range[i], rgb_scale_text[i])
    }
    return(rgb_scale_list)
}

#' @export
get_brewer_set <- function(palette = c("sequential", "diverging", "qualitative")) {
    match.arg(palette,
              several.ok = TRUE)

    sequential_palette <- c('Blues', 'BuGn', 'BuPu', 'GnBu', 'Greens', 'Greys',
                            'Oranges', 'OrRd', 'PuBu', 'PuBuGn', 'PuRd', 'Purples', 'RdPu', 'Reds',
                            'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd')
    names(sequential_palette) <- sequential_palette
    diverging_palette <- c("BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral")
    names(diverging_palette) <- diverging_palette
    qualitative_palette <- c('Accent','Dark2','Paired', 'Pastel1', 'Pastel2','Set1', 'Set2', 'Set3')
    names(qualitative_palette) <- qualitative_palette
    return_palette = list()
    if("qualitative" %in% palette) {
        return_palette <- c(return_palette, as.list(qualitative_palette))
    }
    if("diverging" %in% palette) {
        return_palette <- c(return_palette, as.list(diverging_palette))
    }
    if("sequential" %in% palette) {
        return_palette <- c(return_palette, as.list(sequential_palette))
    }

    return(return_palette)
}

#' @export
compute_umap <- function(input_data, use_dim = 50,
                         n_component=2,
                         metric = "cosine",
                         min_dist = 0.1,
                         n_neighbors = 15L,
                         fast_sgd = FALSE,
                         nn_method = "annoy",
                         cores=1,
                         verbose=T, ...) {
    umap_proj <- uwot::umap(as.matrix(input_data[, 1:use_dim]),
                            n_components = n_component,
                            metric = metric,
                            min_dist = min_dist,
                            n_neighbors = n_neighbors,
                            fast_sgd = fast_sgd,
                            n_threads=cores,
                            verbose=verbose,
                            nn_method = nn_method,
                            ...)
    colnames(umap_proj) <- paste0("UMAP_", 1:n_component)
    rownames(umap_proj) <- rownames(input_data)
    return(umap_proj)
}




#' @export
floor_dec <- function(x, level=1) round(x - 5*10^(-level-1), level)

#' @export
ceiling_dec <- function(x, level=1) round(x + 5*10^(-level-1), level)


#' @export
read_excel_allsheets <- function(filename, tibble = FALSE) {
    # I prefer straight data.frames
    # but if you like tidyverse tibbles (the default with read_excel)
    # then just pass tibble = TRUE
    sheets <- readxl::excel_sheets(filename)
    x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
    if(!tibble) x <- lapply(x, as.data.frame)
    names(x) <- sheets
    x
}


#' @export
leiden_clus <- function (embedding, k = 30) {
    require(RANN)
    require(leiden)
    snn <- RANN::nn2(embedding, k=k)$nn.idx
    adjacency_matrix <- matrix(0L, nrow(embedding), nrow(embedding))
    rownames(adjacency_matrix) <- colnames(adjacency_matrix) <- rownames(embedding)
    for(ii in 1:nrow(embedding)) {
        adjacency_matrix[ii,rownames(embedding)[snn[ii,]]] <- 1L
    }
    partition <- leiden(adjacency_matrix)
    return(partition)
}
