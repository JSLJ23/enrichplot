#' cnetplot
#' 
#' category-gene-network plot
#' @rdname cnetplot
#' @param x input object
#' @param layout network layout
#' @param showCategory selected category to be displayed
#' @param color_category color of category node
#' @param size_category relative size of the category
#' @param color_item color of item node
#' @param size_item relative size of the item (e.g., genes)
#' @param color_edge color of edge
#' @param size_edge relative size of edge
#' @param node_label one of 'all', 'none', 'category', 'item', 'exclusive' or 'share'
#' @param foldChange numeric values to color the item (e.g, foldChange of gene expression values)
#' @param hilight selected category to be highlighted
#' @param hilight_alpha transparent value for not selected to be highlight
#' @param ... additional parameters
#' @importFrom ggtangle cnetplot
#' @method cnetplot enrichResult
#' @export
#' @seealso
#' [cnetplot][ggtangle::cnetplot]
cnetplot.enrichResult <- function(
        x, layout = igraph::layout_with_kk,
        showCategory = 5,
        color_category= "#E5C494", size_category = 1, 
        color_item = "#B3B3B3", size_item = 1, 
        color_edge = "grey", size_edge=.5,
        node_label = "all", 
        foldChange = NULL,
        hilight = "none",
        hilight_alpha = .3,
        ...) {

    geneSets <- extract_geneSets(x, showCategory)
    foldChange <- fc_readable(x, foldChange)    

    p <- cnetplot(geneSets, 
        layout = layout, 
        showCategory = showCategory, 
        foldChange = foldChange, 
        color_category = color_category,
        size_category = size_category,
        color_item = color_item,
        size_item = size_item,
        color_edge = color_edge,
        size_edge = size_edge,
        node_label = node_label,
        hilight = hilight,
        hilight_alpha = hilight_alpha,
        ...
    )

    p <- p + set_enrichplot_color(colors = get_enrichplot_color(3), name = "fold change")
    if (!is.null(foldChange)) {
        p <- p + guides(size  = guide_legend(order = 1), 
                        color = guide_colorbar(order = 2))
    }

    return(p + guides(alpha = "none"))
}

#' @rdname cnetplot
#' @method cnetplot gseaResult
#' @export
cnetplot.gseaResult <- cnetplot.enrichResult

#' @rdname cnetplot
#' @param pie one of 'equal' or 'Count' to set the slice ratio of the pies
#' @method cnetplot compareClusterResult
#' @export
cnetplot.compareClusterResult <- function(
        x, layout = igraph::layout_with_kk,
        showCategory = 5,
        color_category= "#E5C494", size_category = 1, 
        color_item = "#B3B3B3", size_item = 1, 
        color_edge = "grey", size_edge=.5,
        node_label = "all", 
        foldChange = NULL,
        hilight = "none",
        hilight_alpha = .3,
        pie = "equal",
        ...) {

    d <- tidy_compareCluster(x, showCategory)
    y <- split(d$geneID, d$Description)
    gs <- lapply(y, function(item) unique(unlist(strsplit(item, split="/"))))

    if (node_label == "all") {
        node_label = "item"
        add_category_label <- TRUE
    } else if (node_label == "category") {
        node_label = "none"
        add_category_label <- TRUE
    } else {
        add_category_label <- FALSE
    }

    p <- cnetplot(gs, layout = layout, 
        showCategory=length(gs), 
        foldChange = foldChange, 
        color_category = color_category,
        size_category=0, 
        color_item = color_item,
        size_item = size_item,
        color_edge = color_edge,
        size_edge = size_edge,
        node_label = node_label,
        hilight = hilight,
        hilight_alpha = hilight_alpha,        
        ...)

    p <- add_node_pie(p, d, pie, pie_scale=size_category)

    if (add_category_label) {
        p <- p + geom_cnet_label(node_label='category')
    }
    return(p)
}

#' @importFrom ggplot2 coord_fixed
add_node_pie <- function(p, d, pie = "equal", pie_scale = 1) {
    dd <- d[,c('Cluster', 'Description', 'Count')]
    pathway_size <- sapply(split(dd$Count, dd$Description), sum)
    if (pie == "equal") dd$Count <- 1
    dd <- tidyr::pivot_wider(dd, names_from="Cluster", values_from="Count", values_fill=0)
    # dd$pathway_size <- sqrt(pathway_size[dd$Description]/sum(pathway_size))
    dd$pathway_size <- pathway_size[dd$Description]/sum(pathway_size)

    p <- p %<+% dd +
        scatterpie::geom_scatterpie(aes(x=.data$x, y=.data$y, r=.data$pathway_size), 
            cols=as.character(unique(d$Cluster)), 
            legend_name = "Cluster", color=NA, pie_scale = pie_scale) +
        scatterpie::geom_scatterpie_legend(
            dd$pathway_size, x=min(p$data$x), y=min(p$data$y), n=3,
            # labeller=function(x) round(sum(pathway_size) * x^2)
            labeller=function(x) round(x * sum(pathway_size))
        ) +
        coord_fixed() +
        guides(size = "none") 

    return(p)
}

tidy_compareCluster <- function(x, showCategory) {
    d <- fortify(x, showCategory = showCategory, includeAll = TRUE, split = NULL)
    d$Cluster <- sub("\n.*", "", d$Cluster)

    if ("core_enrichment" %in% colnames(d)) { ## for GSEA result
        d$geneID <- d$core_enrichment
    }
    return(d)
}
