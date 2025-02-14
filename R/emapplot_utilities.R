##' Get the similarity matrix
##'
##' @param y A data.frame of enrichment result
##' @param geneSets A list, the names of geneSets are term ids,
##' and every object is a vertor of genes.
##' @param method Method of calculating the similarity between nodes,
##' one of "Resnik", "Lin", "Rel", "Jiang" , "Wang"  and
##' "JC" (Jaccard similarity coefficient) methods
##' @param semData GOSemSimDATA object
##' @noRd
get_similarity_matrix <- function(y, geneSets, method, semData = NULL) {
    id <- y[, "ID"]
    geneSets <- geneSets[id]
    y_id <- unlist(strsplit(y$ID[1], ":"))[1]
    ## Choose the method to calculate the similarity
    if (method == "JC") {
        w <- .cal_jc_similarity(geneSets, id = id, name = y$Description)
        return(w)
    }

    if (y_id == "GO") {
        if(is.null(semData)) {
            stop("The semData parameter is missing,
                and it can be obtained through godata function in GOSemSim package.")
        }
        w <- GOSemSim::mgoSim(id, id, semData=semData, measure=method,
                              combine=NULL)
    }

    if (y_id == "DOID") w <- DOSE::doSim(id, id, measure=method)
    rownames(y) <- y$ID
    rownames(w) <- colnames(w) <- y[colnames(w), "Description"]
    return(w)
}


##' Check whether the similarity matrix exists
##'
##' @param x result of enrichment analysis
##'
##' @noRd
has_pairsim <- function(x) {
    if (length(x@termsim) == 0) {
        error_message <- paste("Term similarity matrix not available.",
            "Please use pairwise_termsim function to",
            "deal with the results of enrichment analysis.")
        stop(error_message)
    }

}


#' Get graph_from_data_frame() result
#'
#' @importFrom igraph graph.empty
#' @importFrom igraph graph_from_data_frame
#' @param enrichDf A data.frame of enrichment result.
#' @param geneSets A list gene sets with the names of enrichment IDs
#' @param color a string, the column name of y for nodes colours
#' @param cex_line Numeric, scale of line width
#' @param min_edge The minimum similarity threshold for whether 
#' two nodes are connected, should between 0 and 1, default value is 0.2.
#' @param pair_sim Semantic similarity matrix.
#' @param method Method of calculating the similarity between nodes,
#' one of "Resnik", "Lin", "Rel", "Jiang" , "Wang"  and
#' "JC" (Jaccard similarity coefficient) methods
#' @return result of graph_from_data_frame()
#' @importFrom igraph V
#' @importFrom igraph 'V<-'
#' @importFrom igraph E
#' @importFrom igraph 'E<-'
#' @importFrom igraph add_vertices
#' @importFrom igraph delete.edges
#' @noRd
build_emap_graph <- function(enrichDf, geneSets, color, cex_line, min_edge,
                             pair_sim, method) {

    if (!is.numeric(min_edge) | min_edge < 0 | min_edge > 1) {
    	stop('"min_edge" should be a number between 0 and 1.')
    }

    if (is.null(dim(enrichDf)) | nrow(enrichDf) == 1) {  # when just one node
        g <- graph.empty(0, directed=FALSE)
        g <- add_vertices(g, nv = 1)
        V(g)$name <- as.character(enrichDf$Description)
        V(g)$color <- "red"
        return(g)
    } else {
        w <- pair_sim[as.character(enrichDf$Description), 
            as.character(enrichDf$Description)]
    }

    wd <- reshape2::melt(w)
    wd <- wd[wd[,1] != wd[,2],]
    # remove NA
    wd <- wd[!is.na(wd[,3]),]
    if (method != "JC") {
        # map id to names
        wd[, 1] <- enrichDf[wd[, 1], "Description"]
        wd[, 2] <- enrichDf[wd[, 2], "Description"]
    }

    g <- graph_from_data_frame(wd[, -3], directed=FALSE)
    E(g)$width <- sqrt(wd[, 3] * 5) * cex_line
    # Use similarity as the weight(length) of an edge
    E(g)$weight <- wd[, 3]
    g <- delete.edges(g, E(g)[wd[, 3] < min_edge])
    idx <- unlist(sapply(V(g)$name, function(x) which(x == enrichDf$Description)))
    cnt <- sapply(geneSets[idx], length)
    V(g)$size <- cnt
    if (color %in% names(enrichDf)) {
        colVar <- enrichDf[idx, color]
    } else {
        colVar <- color
    }
    
    V(g)$color <- colVar
    return(g)
}



##' Get an iGraph object
##'
##' @param x Enrichment result.
##' @param nCategory Number of enriched terms to display.
##' @param color variable that used to color enriched terms, e.g. 'pvalue',
##' 'p.adjust' or 'qvalue'.
##' @param cex_line Scale of line width.
##' @param min_edge The minimum similarity threshold for whether 
##' two nodes are connected, should between 0 and 1, default value is 0.2.
##'
##' @return an iGraph object
##' @noRd
get_igraph <- function(x, nCategory, color, cex_line, min_edge){
    y <- as.data.frame(x)
    geneSets <- geneInCategory(x) ## use core gene for gsea result
    if (is.numeric(nCategory)) {
        y <- y[1:nCategory, ]
    } else {
        y <- y[match(nCategory, y$Description),]
        nCategory <- length(nCategory)
    }

    if (nCategory == 0) {
        stop("no enriched term found...")
    }

    build_emap_graph(enrichDf = y, geneSets = geneSets, color = color,
             cex_line = cex_line, min_edge = min_edge,
             pair_sim = x@termsim, method = x@method)
}


##' Merge the compareClusterResult file
##'
##' @param yy A data.frame of enrichment result.
##'
##' @return a data.frame
##' @noRd
merge_compareClusterResult <- function(yy) {
    yy_union <- yy[!duplicated(yy$ID),]
    yy_ids <- lapply(split(yy, yy$ID), function(x) {
        ids <- unique(unlist(strsplit(x$geneID, "/")))
        cnt <- length(ids)
        list(ID=paste0(ids, collapse="/"), cnt=cnt)
    })

    ids <- vapply(yy_ids, function(x) x$ID, character(1))
    cnt <- vapply(yy_ids, function(x) x$cnt, numeric(1))

    yy_union$geneID <- ids[yy_union$ID]
    yy_union$Count <- cnt[yy_union$ID]
    yy_union$Cluster <- NULL
    yy_union
}

##' add alpha attribute to ggraph edges
##'
##' @param g ggraph object.
##' @param hilight_category category nodes to be highlight.
##' @param alpha_hilight alpha of highlighted nodes.
##' @param alpha_nohilight alpha of unhighlighted nodes.
##' @noRd
edge_add_alpha <- function(g, hilight_category, alpha_nohilight, alpha_hilight) {
    if (!is.null(hilight_category) && length(hilight_category) > 0) {
        edges <- attr(E(g), "vnames")
        E(g)$alpha <- rep(alpha_nohilight, length(E(g)))
        hilight_edge <- grep(paste(hilight_category, collapse = "|"), edges)
        E(g)$alpha[hilight_edge] <- min(0.8, alpha_hilight)
        # E(g)$alpha[hilight_edge] <- alpha_hilight
        } else {
            E(g)$alpha <- rep(min(0.8, alpha_hilight), length(E(g)))
    }
    return(g)
}

##' add alpha attribute to ggraph nodes
##'
##' @param p ggraph object.
##' @param hilight_category category nodes to be highlight.
##' @param hilight_gene gene nodes to be highlight.
##' @param alpha_hilight alpha of highlighted nodes.
##' @param alpha_nohilight alpha of unhighlighted nodes.
##' @noRd
node_add_alpha <- function(p, hilight_category, hilight_gene, alpha_nohilight, alpha_hilight) {
    alpha_node <- rep(1, nrow(p$data))
    if (!is.null(hilight_category)) {
        alpha_node <- rep(alpha_nohilight, nrow(p$data)) 
        hilight_node <- c(hilight_category, hilight_gene)
        alpha_node[match(hilight_node, p$data$name)] <- alpha_hilight
    }
    p$data$alpha <- alpha_node
    return(p)
}







##' Get the location of group label
##'
##' @param ggData data of a ggraph object
##' @param label_format A numeric value sets wrap length, alternatively a
##' custom function to format axis labels.
##' @return a data.frame object.
##' @noRd
get_label_location <- function(ggData, label_format) {
    label_func <- default_labeller(label_format)
    if (is.function(label_format)) {
        label_func <- label_format
    }
    label_x <- stats::aggregate(x ~ color2, ggData, mean)
    label_y <- stats::aggregate(y ~ color2, ggData, mean)
    data.frame(x = label_x$x, y = label_y$y,
        label = label_func(label_x$color2))
}



##' Cluster similar nodes together by k-means
##' 
##' @param p a ggraph object.
##' @param enrichDf data.frame of enrichment result.
##' @param nWords Numeric, the number of words in the cluster tags.
##' @param clusterFunction function of Clustering method, such as stats::kmeans, cluster::clara,
##' cluster::fanny or cluster::pam.
##' @param nCluster Numeric, the number of clusters, 
##' the default value is square root of the number of nodes.
##' @noRd
groupNode <- function(p, enrichDf, nWords, clusterFunction =  stats::kmeans, nCluster) {
    ggData <- p$data
    wrongMessage <- paste("Wrong clusterFunction parameter or unsupported clustering method;",
         "set to default `clusterFunction = kmeans`")
    if (is.character(clusterFunction)) {
        clusterFunction <- eval(parse(text=clusterFunction))
    }   
    if (!"color2" %in% colnames(ggData)) {
        dat <- data.frame(x = ggData$x, y = ggData$y)
        nCluster <- ifelse(is.null(nCluster), floor(sqrt(nrow(dat))), 
                           min(nCluster, nrow(dat)))
        ggData$color2 <- tryCatch(expr  = clusterFunction(dat, nCluster)$cluster, 
                                  error = function(e) {
                                      message(wrongMessage)
                                      stats::kmeans(dat, nCluster)$cluster
                                  })
        if (is.null(ggData$color2)) {
            message(wrongMessage)
            ggData$color2 <- stats::kmeans(dat, nCluster)$cluster
        }
    }
    goid <- enrichDf$ID
    cluster_color <- unique(ggData$color2)
    clusters <- lapply(cluster_color, function(i){goid[which(ggData$color2 == i)]})
    cluster_label <- sapply(cluster_color, get_wordcloud, ggData = ggData,
                            nWords=nWords)
    names(cluster_label) <- cluster_color
    ggData$color2 <- cluster_label[as.character(ggData$color2)] 
    return(ggData)
}

#' add ellipse to group the node
#'
#' @param p ggplot2 object
#' @param group_legend Logical, if TRUE, the grouping legend will be displayed.
#' The default is FALSE.
#' @param label_style style of group label, one of "shadowtext" and "ggforce".
#' @param ellipse_style style of ellipse, one of "ggforce" an "polygon".
#' @param ellipse_pro numeric indicating confidence value for the ellipses
#' @param alpha the transparency of ellipse fill.
#' @importFrom rlang check_installed
#' @importFrom ggplot2 scale_fill_discrete
#' @noRd
add_ellipse <- function(p, group_legend, label_style, 
    ellipse_style = "ggforce", ellipse_pro = 0.95, alpha = 0.3, ...) {
    show_legend <- c(group_legend, FALSE)
    names(show_legend) <- c("fill", "color") 
    ellipse_style <- match.arg(ellipse_style, c("ggforce", "polygon"))
    
    check_installed('ggforce', 'for `add_ellipse()`.')
    
    if (ellipse_style == "ggforce") {
        if (label_style == "shadowtext") {
            p <- p + ggforce::geom_mark_ellipse(
                        aes(x = .data$x, y = .data$y, fill = .data$color2), 
                        alpha = alpha, color = NA, show.legend = show_legend)
        } else {
            p <- p + ggforce::geom_mark_ellipse(
                        aes(x = .data$x, y = .data$y, fill = .data$color2, label = .data$color2), 
                        alpha = alpha, color = NA, show.legend = show_legend)
        }
        if (group_legend) p <- p + scale_fill_discrete(name = "groups")
     } 
    
    if (ellipse_style == "polygon") {
        p <- p + ggplot2::stat_ellipse(aes_(x =~ x, y =~ y, fill =~ color2),
                                       geom = "polygon", level = ellipse_pro,
                                       alpha = alpha,
                                       show.legend = group_legend, ...)
    }

    return(p)
}



list2df <- ggtangle:::list2df
