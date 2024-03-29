################################################################################
# nice theme
################################################################################

#' @title nice ggplot theme
#' @description nice theme with square border
#' @importFrom ggplot2 theme element_blank element_rect
#' @export

theme_rect <-function() {
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_rect(color = "black", size = 1, fill = NA),
          aspect.ratio = 1)
}

################################################################################
# feature plots with predefined marker list
################################################################################

#' @title Seurat feature plot
#' @description create and save a nice Seurat feature plot in folder `featureplot`
#' @param path path to markers.csv
#' @param object Seurat object
#' @param par column name in markers.csv
#' @param reduction a character string specifying the dimension reduction
#' @param width width of output plot (default: 16)
#' @param height height of output plot (default: length of genes divided by four, ceiling, times three)
#' @param order should the feature plot be ordered in order of expression 
#' @return save feature plot to folder `/results/featureplot/`
#' @importFrom ggplot2 theme element_blank element_rect ggsave
#' @examples \dontrun{fPlot(sc_merge, path = file.path("lookup", "markers.csv"), par = "main", reduction = "umap")}
#' @export

fPlot <- function(path, object, par, reduction,  width = 16, height = ceiling(length(genes_found) / 4) * 3, order) {
  if (!methods::is(object) == "Seurat") {
    stop("Object must be a Seurat object")
  }
  dir.create(file.path("results", "featureplot"), showWarnings = FALSE)
  markers <- readr::read_csv(path) |>
    as.list(markers) |>
    lapply(function(x) x[!is.na(x)])
  genes <- markers[[par]]
  if (is.null(genes)) {
    stop("No genes were found. Make sure that `par` exists in markers.csv")
  }
  available_genes <- rownames(object)
  genes_found <- genes[genes %in% available_genes]
  object_parse <- deparse(substitute(object))
  fp <- Seurat::FeaturePlot(object = object, features = unique(genes), cols = c("#F0F0F0", "#CB181D"), reduction = reduction, pt.size = .1, order = order, coord.fixed = TRUE, ncol = 4, raster = FALSE) &
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_rect(color = "black", size = 1, fill = NA)
    )
  ggsave(filename = file.path("results", "featureplot", glue::glue("fp_{object_parse}_{par}.png")), width = width, height = height, limitsize = FALSE)
}


################################################################################
# feature plots with custom marker list
################################################################################

#' @title Seurat feature plot
#' @description create and save a Seurat feature plot in folder `featureplot`
#' @param object Seurat object
#' @param markers a data frame with a column called `cell_source` that represents the cell population and its source and a column `gene`
#' @param par a character string representing the cell_source to plot
#' @param reduction a character string specifying the dimension reduction
#' @param width width of output plot (default: 16)
#' @param height height of output plot (default: length of genes divided by four, ceiling, times three)
#' @return save feature plot to folder `/results/featureplot/`
#' @importFrom ggplot2 theme element_blank element_rect ggsave
#' @examples \dontrun{fPlot(sc_merge, par = "main", filepath = file.path("results", "featureplot", glue::glue("fp_")))}
#' @export

fPlotCustom <- function(object, markers, par, reduction, width = 16, height = ceiling(length(genes_found) / 4) * 3) {
  if (!methods::is(object) == "Seurat") {
    stop("Object must be a Seurat object")
  }
  dir.create(file.path("results", "featureplot"), showWarnings = FALSE)
  genes <- markers[markers$cell_source == par, ]$gene
  genes <- genes[!is.na(genes)]
  available_genes <- rownames(object)
  genes_found <- genes[genes %in% available_genes]
  object_parse <- deparse(substitute(object))
  fp <- Seurat::FeaturePlot(object = object, features = unique(genes), cols = c("#F0F0F0", "#CB181D"), reduction = reduction, pt.size = .1, order = FALSE, coord.fixed = TRUE, ncol = 4, raster = FALSE, alpha = 0.2) &
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_rect(color = "black", size = 1, fill = NA)
    )
  ggsave(filename = file.path("results", "featureplot", glue::glue("fp_{object_parse}_{par}.png")), width = width, height = height, limitsize = FALSE)
}

################################################################################
# dot plots
################################################################################

#' @title nice Seurat dot plot
#' @description create and save a nice Seurat dot plot
#' @param path path to markers.csv
#' @param object Seurat object
#' @param par column name in markers.csv
#' @param dot_min minimal dot size
#' @param width width of output plot (default: 10)
#' @param height height of output plot (default: 10)
#' @param ortho convert to orthologues? Allowed values: `none`, `mouse2human` or `human2mouse`
#' @param scale should the values be scaled? (default: TRUE)
#' @return save dot plot to folder `results/dotplot/`
#' @importFrom ggplot2 ggplot scale_size theme xlab ylab element_text ggsave
#' @examples
#' \dontrun{
#' dotPlot(object = sc_merge, par = "cellmarkers_covid", dot_min = 0.1)
#' }
#' @export

dotPlot <- function(path, object, par, dot_min, scale = TRUE, ortho = "none", width = 10, height = 10) {
    if(!methods::is(object) == "Seurat") {
        stop("Object must be a Seurat object")
    }
    dir.create(file.path("results", "dotplot"), showWarnings = FALSE)
    markers <- readr::read_csv(path) |>
        as.list(markers) |>
        lapply(function(x) x[!is.na(x)])
    genes <- markers[[par]]

    if(is.null(genes)) {
        stop("No genes were found. Make sure that `par` exists in markers.csv")
}
    if(!(ortho %in% c("none", "mouse2human", "human2mouse"))) {
        stop("ortho must take values: `none`, `mouse2human` or `human2mouse`")
    }
    if (ortho == "mouse2human") {
        genes <- homologene::mouse2human(genes, db = homologene::homologeneData2)$humanGene
        message("genes converted from mouse to human")
    } else if (ortho == "human2mouse") {
        genes <- homologene::human2mouse(genes, db = homologene::homologeneData2)$mouseGene
        message("genes converted from human to mouse")
    } else if (ortho == "none") {
        message("no genes were converted")
    }
    object_parse <- deparse(substitute(object))
    dp <- Seurat::DotPlot(object, features = unique(genes), dot.scale = 10, scale.by = "size", dot.min = dot_min, scale = scale) +
        viridis::scale_color_viridis(option = "viridis") +
        scale_size(range=c(0,10))+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "italic"))+
        xlab(NULL)+
        ylab(NULL)
    ggsave(file.path("results", "dotplot", glue::glue("dp_{object_parse}_{par}.pdf")), width = width, height = height, limitsize = FALSE)
}


################################################################################
# pheatmap
################################################################################

#' @title nice pheatmap
#' @description create and save a nice pheatmap in the folder `heatmap` using color breaks
#' @param matrix numeric matrix of the values to be plotted
#' @param scale should the values be centered and scaled in row, column direction? Allowed: `row`, `column`, `none` (default: `none`)
#' @param width width of output plot (default: number of columns divided by two)
#' @param height height of output plot (default: number of rows divided by three)
#' @param cellwidth individual cell width (default: 10)
#' @param cellheight individual cell height (default: 10)
#' @param treeheight_row height of a tree for rows (default: 10)
#' @param treeheight_col height of a tree for columns (default: 10)
#' @param fontsize fontsize (default: 10)
#' @param cluster_rows cluster rows? (default: true)
#' @param cluster_cols cluster columns? (default: true)
#' @param annotation_row data frame that contains the annotations. Rows in the data and in the annotation are matched using row names. (default: NA)
#' @return save heatmap to folder `/results/heatmap`
#' @examples \dontrun{pHeatmap(szabo_tc_tc_avg, scale = "row", cluster_cols = FALSE)}
#' @export

pHeatmap <- function(matrix, scale = "none", height = ceiling(nrow(matrix)/3), width = ceiling(ncol(matrix)/2), cellwidth = 10, cellheight = 10, treeheight_row = 10, treeheight_col = 10, fontsize = 10, cluster_rows = TRUE, cluster_cols = TRUE, annotation_row = NA) {
    dir.create(file.path("results", "heatmap"), showWarnings = FALSE)
    matrix_parse <- deparse(substitute(matrix))
    matrix <- matrix[!rowSums(matrix) == 0,] # filter rows with only zeros
    break_max <- round(max(abs(c(max(scale_mat(matrix, scale = scale)), min(scale_mat(matrix, scale = scale)))))-0.1,1) #use internal function to get scaled matrix and max value for color legend
    break_min <- -break_max
    phmap <- pheatmap::pheatmap(matrix,
                                color = viridis::viridis(100),
                                scale = scale,
                                cellwidth = cellwidth,
                                cellheight = cellheight,
                                fontsize = fontsize,
                                treeheight_row = treeheight_row,
                                treeheight_col = treeheight_col,
                                cluster_rows = cluster_rows,
                                cluster_cols = cluster_cols,
                                clustering_method = "complete",
                                border_color = NA,
                                legend_breaks = seq(break_min, break_max, length = 3),
                                annotation_row = annotation_row
                                )
    pdf(file.path("results", "heatmap", glue::glue("hm_{matrix_parse}.pdf")), width = width, height = height)
    print(phmap)
    dev.off()
}

################################################################################
# stackedPlot
################################################################################

#' @title abundance stacked bar plot
#' @description create and save an abundance stacked bar plot in the folder `abundance` 
#' @param object Seurat object
#' @param x_axis variable in meta data that is used for the y axis
#' @param y_axis variable in meta data that is used for the x axis
#' @param x_order vector determining the order of the x axis
#' @param y_order vector determining the order of the y axis
#' @param color color palette
#' @param width width of output plot (default: 10)
#' @param height height of output plot
#' @return save stacked abundance barplot to folder `/results/abundance`
#' @examples
#' \dontrun{
#' stackedPlot(
#' object = sc_merge,
#' x_axis = "pool",
#' y_axis = "cluster",
#' x_order = unique(sc_merge$pool),
#' y_order = cluster_order,
#' color = col_vector,
#' width = 4)
#' }
#' @export

stackedPlot <- function(object, x_axis, y_axis, x_order, y_order, color, width, height = 10) {
    if(!methods::is(object) == "Seurat") {
        stop("Object must be a Seurat object")
    }
    dir.create(file.path("results", "abundance"), showWarnings = FALSE)
    object_parse <- deparse(substitute(object))
    result_wide <- as.data.frame.matrix(table(object@meta.data[[y_axis]], object@meta.data[[x_axis]])) |>
        rownames_to_column("cell") |>
        mutate(across(where(is.numeric), function(x) x/sum(x)*100))
    result_long <- result_wide |>
        pivot_longer(!cell, names_to = "type", values_to = "count") |>
        mutate(cell = factor(cell, levels = y_order)) |>
        mutate(type = factor(type, levels = x_order)) |>
        dplyr::filter(count != 0)
    sbp <- ggplot(data = result_long)+
        geom_col(aes(x = type, y = count, fill = cell), color = "black", size = 0.1, position = "fill")+
        scale_fill_manual(values = color)+
        guides(fill = guide_legend(title = NULL))+ #remove guide label
        theme_classic()+ #remove background
        ylab("Proportion of cells")+
        xlab("")+
        theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.3))
    ggsave(file.path("results", "abundance", glue::glue("stacked_barplot_{object_parse}_{x_axis}.pdf")), sbp, width = width, height = height)
}

################################################################################
# abundance Volcano Plot
################################################################################

#' @title abundance volcano plot
#' @description create and save an abundance volcano bar plot in the folder `abundance` 
#' @param object Seurat object
#' @param cluster_idents variable in meta data with cluster names
#' @param sample variable in meta data for each sample
#' @param cluster_order vector determining the order of the clusters
#' @param group_by variable in meta data that categorize samples in groups
#' @param group1 first group (nominator)
#' @param group2 second group (denominator)
#' @param color color palette
#' @param width width of output plot (default: 5)
#' @param height height of output plot (default: 5)
#' @param min_cells remove all clusters that have less than minimal amount of cells (default = 10)
#' @param paired logical indicating whether you want a paired test (default FALSE)
#' @return save volcano abundance plot to folder `/results/abundance`
#' @examples
#' \dontrun{
#' abVolPlot(object = aie_sct,
#          cluster_idents = "predicted.id",
#          sample = "sample",
#          cluster_order = unique(aie_sct$predicted.id),
#          group_by  = "AIE_type",
#          group1 = "LGI1",
#          group2 = "control",
#          color = dittoColors(), 
#          min_pct = 0.5)
#' }
#' @export


abVolPlot <- function(object, cluster_idents, sample, cluster_order, group_by, group1, group2, color, width = 5, height = 5, min_cells = 10, paired = FALSE) {
    if(!methods::is(object) == "Seurat") {
        stop("Object must be a Seurat object")
    }
    dir.create(file.path("results", "abundance"), showWarnings = FALSE)
    object_parse <- deparse(substitute(object))

    cl_size_ind <- as.data.frame.matrix(table(object@meta.data[[cluster_idents]], object@meta.data[[sample]])) |>
        rownames_to_column("cluster") |>
        mutate(across(where(is.numeric), function(x) x/sum(x)*100)) |>
        pivot_longer(!cluster, names_to = "sample", values_to = "count") |>
        left_join(unique(tibble(sample = object@meta.data[[sample]], group_by = object@meta.data[[group_by]]))) |>
        dplyr::filter(group_by == group1 | group_by == group2)

    pvalue_res <- vector("double") # define output

    for (i in cluster_order) {
        out1 <- cl_size_ind[cl_size_ind$cluster == i,]
        pvalue_res[i] <- wilcox.test(count ~ group_by, data = out1, paired = paired)$p.value
    }

                                        #wilcox_res <- p.adjust(wilcox_res, "BH")
    pvalue_cl <- data.frame(cluster = cluster_order, pvalue = pvalue_res)

    cl_size <- as.data.frame.matrix(table(object@meta.data[[cluster_idents]], object@meta.data[[group_by]])) |>
        rownames_to_column("cluster") |> 
        mutate(cluster = factor(cluster, levels = cluster_order)) |>
        dplyr::filter(.data[[group1]] > min_cells & .data[[group2]] > min_cells) |>
        mutate(across(where(is.numeric), function(x) x/sum(x)*100)) |>
        mutate(logratio = log2(.data[[group1]]/.data[[group2]])) |>
        left_join(pvalue_cl, by = "cluster") |>
        mutate(log_pvalue = -log10(pvalue))
      

    p1 <- ggplot(cl_size, aes(x = logratio, y = log_pvalue, color = cluster, size = 3, label = cluster))+
        geom_point()+
        scale_color_manual(values = color)+
        theme_classic()+
        ggrepel::geom_text_repel(nudge_y = 0.07)+
        geom_hline(yintercept = -log10(0.05), color = "blue", linetype = "dashed")+ #horizontal line p unadjusted
        geom_hline(yintercept = -log10(0.05/nrow(cl_size)), color = "blue")+
        geom_vline(xintercept = 0, color = "red", linetype = "solid")+
        geom_vline(xintercept = -1, color = "red", linetype = "dashed")+
        geom_vline(xintercept = 1, color = "red", linetype = "dashed")+ 
        xlab(bquote(~Log[2]~ 'fold change'))+
        ylab(bquote(-Log[10]~ "p value")) +
        theme(legend.position = "none") #remove guide
    ggsave(file.path("results", "abundance", glue::glue("volcano_plot_{cluster_idents}_{object_parse}_{group1}_{group2}.pdf")), width = width, height = height)
}

################################################################################
# internal helper function to calculcate significance for boxplots
################################################################################

#' @title significance for boxplot
#' @description calculcate significance for ggsignif
#' @param x_var numeric variable (x in the formula `x ~ group`)
#' @param group grouping variable (group in the formula `x ~ group`)
#' @param data data frame containing the variables
#' @param paired logical value. Do you want to do perform pair testing?
#' @return data frame with adjusted signficance values
#' @examples
#' \dontrun{compStat(x_var = "pct", group = "type", data = bp_data)}

compStat <- function(x_var, group, data, paired) {
# initalize stats
  stats <- vector("list")

  # for character run pairwise fisher test for all parameters, only keep important columns so they match
  for (par in x_var) {
    f_str <- paste0(par, "~", group)
    if (length(unique(data[[par]])) > 2) {
      if(paired == FALSE) {
        stats[[par]] <- rstatix::dunn_test(as.formula(f_str), data = data, p.adjust.method = "none") |>
          dplyr::select(.y., group1, group2, p, p.adj, p.adj.signif)
      }
      if (paired == TRUE) {
        stats[[par]] <- rstatix::wilcox_test(as.formula(f_str), data = data, p.adjust.method = "none", paired = paired) |>
          dplyr::select(.y., group1, group2, p) |>
          dplyr::mutate(p.adj = NA,
                        p.adj.signif = NA)
      }
    } else {
      stats[[par]] <- rstatix::wilcox_test(as.formula(f_str), data = data, p.adjust.method = "none", paired = paired) |>
        dplyr::select(.y., group1, group2, p) |>
        dplyr::mutate(p.adj = NA,
                      p.adj.signif = NA)
    }
  }
  # combine these lists into a dataframe, do p value adjustment with BH, and then extract only significant values
  stats_df <-
    do.call(rbind, stats) |>
    dplyr::mutate(p.adj = p.adjust(p, method = "BH")) |>
    dplyr::filter(p.adj < 0.05) |>
    dplyr::mutate(p.adj.signif = as.character(symnum(p.adj, corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " "))))
  return(stats_df)
}

################################################################################
# abundance boxplot
################################################################################

#' @title abundance boxplot plot
#' @description create and save an abundance boxplot in the folder `abundance`
#' @param object Seurat object
#' @param cluster_idents variable in meta data with cluster names
#' @param sample variable in meta data for each sample
#' @param cluster_order vector determining the order of the clusters
#' @param group_by variable in meta data that categorize samples in groups
#' @param group_order vector determining the order of the samples
#' @param color color palette
#' @param width width of output plot (default: 9)
#' @param height height of output plot (default: length of cluster_idents divided by four, ceiling, times three)
#' @param paired logical indicating whether you want a paired test (default FALSE)
#' @return save abundance box plot in the folder `/results/abundance`
#' @examples
#' \dontrun{
#'abBoxPlot(object = aie_pbmc,
#'          cluster_idents = "cluster",
#'          sample = "sample",
#'        cluster_order = cluster_order,
#'          group_by = "AIE_type",
#'          group_order = c("control", "CASPR2", "LGI1"),
#'          color = my_cols)
#' }
#' @export


abBoxPlot <- function(object, cluster_idents, sample, cluster_order, group_by, group_order, color, width = 9, height = ceiling(length(unique(object@meta.data[[cluster_idents]]))/4)*3, paired = FALSE, number_of_tests) {
  if(!methods::is(object) == "Seurat") {
    stop("Object must be a Seurat object")
  }
  dir.create(file.path("results", "abundance"), showWarnings = FALSE)
  object_parse <- deparse(substitute(object))
  bp_data <-
    table(object@meta.data[[cluster_idents]], object@meta.data[[sample]]) |>
    as.data.frame.matrix() |>
    dplyr::select(where(function(x) any(x != 0))) |> # filter out columns with only zeros
    tibble::rownames_to_column("cluster") |>
    dplyr::mutate(across(where(is.numeric), function(x) x/sum(x)*100)) |>
    tidyr::pivot_longer(where(is.numeric), names_to = "sample", values_to = "pct") |>
    dplyr::mutate(cluster = factor(cluster, levels = cluster_order)) |>
    dplyr::left_join(unique(tibble(sample = object@meta.data[[sample]], type = object@meta.data[[group_by]]))) |>
    dplyr::mutate(type = factor(type, levels = group_order)) |>
    tidyr::pivot_wider(names_from = cluster, values_from = pct)

  ## |>
  ##   dplyr::group_split(cluster) |>
  ##   setNames(cluster_order)

  # calculate stats for all tests and adjust
  bp_stats <- compStat(x_var = names(bp_data)[-c(1,2)], group = "type", data = bp_data, paired = paired)

  bp_data_plot <- bp_data |>
    dplyr::mutate(patient = gsub(sample, pattern = "(\\w+)_(\\d+)" , replacement = "\\1" ))

  bp_plot <- vector("list")


  for (var in names(bp_data)[-c(1,2)]) {
    #extract stats and data for each plot
    stats_df <- dplyr::filter(bp_stats, .y. == var)
    stats_list <- vector("list")
    if(nrow(stats_df) != 0) {
      stats_list$annotation <- stats_df$p.adj.signif
      for (i in 1:nrow(stats_df)) {
        stats_list$comparisons[[i]] <- c(stats_df$group1[i], stats_df$group2[i])
      }
    }
    # create plot for each variable
    bp_plot[[var]] <-
      bp_data_plot |>
      ggplot(aes(x = type, y = .data[[var]])) +
      ggsignif::geom_signif(comparisons = stats_list$comparisons, annotation = stats_list$annotation, textsize = 5, step_increase = 0.05, vjust = 0.7)+
      geom_boxplot(aes(fill = type)) +
      geom_point() +
      geom_line(aes(group = patient)) +
      theme_bw()+
      ggtitle(var) +
      theme(legend.position = "none")+
      xlab("") +
      ylab("percentage")+
      scale_fill_manual(values = color)
  }
  patchwork::wrap_plots(bp_plot, ncol = 4)
  ggsave(file.path("results", "abundance", glue::glue("boxplot_{cluster_idents}_{object_parse}_{group_by}.pdf")), width = width, height = height)
}

################################################################################
# violin/boxplot module plot
################################################################################

#' @title module plot
#' @description create combined violin and boxplot of the module score
#' @param x_var variable in meta data which represents the x-axis
#' @param module variable in meta data of the module score
#' @param object Seurat object
#' @param color color palette
#' @return plot module plot
#' @examples
#' \dontrun{ModulePlot(object = aie_csf, x_var = "AIE_type", module = "TCRVG1", color = my_cols)}
#' @export

ModulePlot <- function(x_var, module, object, color) {
data_module <- tibble(x_axis = object@meta.data[[x_var]], module = object@meta.data[[module]])
signif <- vector("list")
f_str <- paste0("x_axis" ~ "module")
if(length(unique(data_module$x_axis)) > 2) {
    stats <- rstatix::dunn_test(as.formula(f_str), data = data_module, p.adjust.method = "BH") |>
        dplyr::filter(p.adj < 0.05) |>
        mutate(p.adj.signif = as.character(symnum(p.adj, corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " "))))

    signif$annotation <- stats$p.adj.signif
} else {
        stats <- rstatix::wilcox_test(as.formula(f_str), data = data_module) |>
            dplyr::filter(p < 0.05) |>
            mutate(p.signif = as.character(symnum(p, corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " "))))
        signif$annotation <- stats$p.signif
}

if(nrow(stats) != 0) {
    for (i in 1:nrow(stats)) {
        signif$comparisons[[i]] <-c(stats$group1[i], stats$group2[i])
    }
} else {
    signif <- list()
}

module_plot <-
    ggplot(data_module, aes(x = x_axis, y = module))+
    geom_violin(aes(fill = x_axis))+
    geom_boxplot(width = .15)+
#    stat_summary(fun = mean, geom = "point")+
    scale_fill_manual(values = color)+
    theme_bw()+
    xlab("")+
    ylab("")+
    ggtitle(module)+
    theme(legend.position = "none")
    ## ggsignif::geom_signif(comparisons = signif$comparisons, 
    ##                       annotation = signif$annotation, textsize = 5, 
    ##                       step_increase = 0.05, vjust = 0.7)

return(module_plot)
}


################################################################################
# enrichr plot
################################################################################

#' @title enrichr plot
#' @description plot enrichr results in a bar plot
#' @param filename name of file without extension and without `enrichr_`
#' @param sheet name of the excel sheet
#' @param width width of output plot
#' @param height height of output plot
#' @return save enrichr plot in the folder `/results/enrichr`
#' @importFrom ggplot2 ggplot theme labs ggsave aes geom_col theme_classic
#' @examples
#' \dontrun{plotEnrichr("de_ALZ_Naive_CSF_neg_pDC", sheet = "GO_Biological_Process_2021", width = 10, height = 5)}
#' @export
plotEnrichr <- function(filename, sheet, width, height) {
    dir.create(file.path("results", "enrichr"), showWarnings = FALSE)
    colors <- RColorBrewer::brewer.pal(5, "Set2")
    color <- ifelse(grepl(x = filename, pattern = "pos"), colors[[1]], colors[[2]])
    enrichr <- readxl::read_excel(file.path("results", "enrichr", glue::glue("enrichr_{filename}.xlsx")), sheet = sheet) |>
        dplyr::filter(Adjusted.P.value < 0.05) |>
        dplyr::slice_min(order_by = Adjusted.P.value, n = 10, with_ties = FALSE) |>
        tidyr::separate(Overlap, into = c("overlap1", "overlap2")) |>  # separate overlap in two columns
        dplyr::mutate(overlap = as.numeric(overlap1)/as.numeric(overlap2)) |> # calculcate overlap
        ggplot(aes(y = reorder(Term, -log10(Adjusted.P.value)), x = -log10(Adjusted.P.value))) +
        geom_col(fill = color)+
        labs(x = "-Log10 Adjusted P value",
             y = "")+
        theme_classic()+
        theme(legend.position = "none")
    ggsave(file.path("results", "enrichr", glue::glue("barplot_enrichr_{filename}_{sheet}.pdf")), width = width, height = height)
}

################################################################################
# abundance propeller plot volcano
################################################################################
#' @title plot propeller results
#' @description The function creates a volcano plot of the propeller results and saves the plot in results/abundance folder
#' @param data A dataframe containing the results from propeller calculation
#' @param color A vector of colors for the clusters in the plot
#' @param filename A character representing the file name of the plot
#' @param width The width of the plot
#' @param height The height of the plot
#' @param FDR The FDR threshold for the plot
#' @importFrom ggplot2 ggplot theme labs ggsave aes geom_col theme_classic
#' @examples
#' \dontrun{plotPropeller(data = pnp_ctrl_csf_sex_age, color = cluster_col, filename = "pnp_ctrl_csf_sex_age")}
#' @export

plotPropeller <- function(data, color, filename, width = 5, height = 5, FDR) {
  dir.create(file.path("results", "abundance"), showWarnings = FALSE)
  ggplot(data, aes(x = log2ratio, y = FDR_log, color = cluster, size = 3, label = cluster)) +
    geom_point() +
    scale_color_manual(values = color) +
    theme_classic() +
    ggrepel::geom_text_repel(nudge_y = 0.07, max.overlaps = 20) +
    geom_hline(yintercept = -log10(FDR), color = "blue", linetype = "dashed") + # horizontal line p unadjusted
    geom_vline(xintercept = 0, color = "red", linetype = "solid") + # vertical line
    geom_vline(xintercept = -1, color = "red", linetype = "dashed") + # vertical line
    geom_vline(xintercept = 1, color = "red", linetype = "dashed") + # vertical line
    xlab(bquote(~ Log[2] ~ "fold change")) +
    ylab(bquote(~ -Log[10] ~ "adjusted p value")) +
    theme(legend.position = "none") # remove guide
  ggsave(file.path("results", "abundance", glue::glue("propeller_{filename}.pdf")), width = width, height = height)
}

################################################################################
# abundance propeller plot barplot
################################################################################
#' @title plot propeller results in a barplot
#' @description The function creates a barplot of the propeller results and saves the plot in results/abundance folder
#' @param data A dataframe containing the results from propeller calculation
#' @param color A vector of colors for the clusters in the plot
#' @param filename A character representing the file name of the plot
#' @param width The width of the plot
#' @param height The height of the plot
#' @importFrom ggplot2 ggplot theme labs ggsave aes geom_col theme_classic
#' @examples
#' \dontrun{dotplotPropeller(data = pnp_ctrl_csf_sex_age, color = cluster_col, filename = "pnp_ctrl_csf_sex_age")}
#' @export

dotplotPropeller <- function(data, color, filename, width = 5, height = 5) {
  dir.create(file.path("results", "abundance"), showWarnings = FALSE)
  ggplot(data, aes(x = log2ratio, y = fct_reorder(cluster, log2ratio), color = cluster)) +
    geom_point(size = 5) +
    theme_classic() +
    geom_vline(xintercept = 0, color = "red", linetype = "solid") + # vertical line
    # geom_vline(xintercept = -1, color = "red", linetype = "dashed") + # vertical line
    # geom_vline(xintercept = 1, color = "red", linetype = "dashed") + # vertical line
    scale_color_manual(values = color) +
    xlab("Log2 fold change") +
    ylab(NULL) +
    theme(legend.position = "none") # remove legend
  ggsave(file.path("results", "abundance", glue::glue("propeller_dotplot_{filename}.pdf")), width = width, height = height)
}

################################################################################
# slingshot visualization
################################################################################
#' @title plot slingshot results
#' @description The function creates a colored UMAP plot and curves split by lineage

#' @param object A Seurat object
#' @param lineage A character string indicating the lineage of interest.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{sds_plots_list <- lapply(colnames(pt), slingshotPlot, object = bcells)}
#' @importFrom ggplot2 aes geom_point geom_path ggtitle theme_classic element_blank element_rect
#' @export

plotSlingshot <- function(object, lineage) {
    if(!methods::is(object) == "Seurat") {
        stop("Object must be a Seurat object")
    }
    my_curves <- dplyr::filter(curves, Lineage == lineage)
sds_plot <-
    Embeddings(object, "umap") |>
    tibble::as_tibble() |> 
    dplyr::mutate(color = pt[,lineage])|>
    tidyr::drop_na(color) |>
    ggplot(aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = color), size = 0.1)+
    geom_path(data = my_curves, aes(group = Lineage))+
    viridis::scale_color_viridis(name = "pseudotime")+
    theme_classic() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_rect(color = "black", linewidth = 1, fill = NA),
          aspect.ratio = 1)+
    ggtitle(lineage)
return(sds_plot)
}


################################################################################
# PCA of cluster abundances
################################################################################
#' @title plot PCA of cluster abundances
#' @description The function creates a PCA plot of the cluster abundances and saves the plot in results/abundance folder

#' @param object A Seurat object
#' @param cluster A character string indicating the cluster column in the metadata of the Seurat object
#' @param sample A character string indicating the sample column in the metadata of the Seurat object
#' @param condition A character string indicating the condition column in the metadata of the Seurat object
#' @param width The width of the plot
#' @param height The height of the plot
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' pcaSeurat(
#'   object = sc_final,
#'   cluster = "cluster",
#'   sample = "patient",
#'   condition = "condition")
#' }
#' @importFrom ggplot2 aes geom_point geom_path ggtitle theme_classic element_blank element_rect
#' @export

pcaSeurat <-function(object, cluster, sample, condition, width = 20, height = 5) {
  dir.create(file.path("results", "pca"), showWarnings = FALSE)
  object_parse <- deparse(substitute(object))
  cl_size <-
    as.data.frame.matrix(table(object@meta.data[[cluster]], object@meta.data[[sample]])) |>
    t()

  colnames(cl_size) <- levels(object@meta.data[[cluster]])

  pca_result <-FactoMineR::PCA(cl_size, scale.unit = TRUE, ncp = 30, graph = FALSE)
  pca_eigen <- factoextra::fviz_eig(pca_result, addlabels = TRUE, ylim = c(0,50), ncp = 7)

pca_var_plot <-
  factoextra::fviz_pca_var(
    pca_result,
    col.var = "contrib",
    gradient.cols = viridis::viridis(100),
    repel = TRUE
  ) +
    labs(title = "")+
    theme_classic()
  # select.var = list(contrib = 15)) # top 15

  lookup_pre <-
    data.frame(
      cluster = object@meta.data[[sample]],
      condition = object@meta.data[[condition]]
    ) |>
    distinct()

  lookup <-
    data.frame(cluster = rownames(cl_size)) |>
    left_join(lookup_pre)

factoextra::fviz_pca_ind(pca_result)

  pca_plot_ind <-
    factoextra::fviz_pca_ind(
      pca_result,
      pointsize = 5,
      pointshape = 21,
      fill = "#E7B800",
      col.ind = "black",
      palette = "Set2",
      axes.linetype = "solid"
    )

  pca_ggplot_ind <-
    ggpubr::ggpar(
    pca_plot_ind,
    title = "",
    xlab = "PC1",
    ylab = "PC2",
    ggtheme = theme_bw() +
      theme(axis.title.x = element_text(size=15),
            axis.title.y = element_text(size=15),
            plot.title = element_text(size=25))) +
    theme_classic()

  pca_plot_group <-
    factoextra::fviz_pca_ind(
      pca_result,
      pointsize = 5,
      pointshape = 21,
      geom.ind = "point",
      fill.ind = lookup$condition,
      col.ind = "black",
      palette = "Set2",
      addEllipses = TRUE,
      ellipse.type = "confidence",
      legend.title = "group",
      axes.linetype = "solid"
    )

  pca_ggplot_group <-
    ggpubr::ggpar(
      pca_plot_group,
      title = "",
      xlab = "PC1",
      ylab = "PC2",
      ggtheme = theme_bw() +
        theme(axis.title.x = element_text(size=15),
              axis.title.y = element_text(size=15),
              plot.title = element_text(size=25))) +
    theme_classic()

  pca_plots <- patchwork::wrap_plots(pca_eigen, pca_var_plot, pca_ggplot_ind, pca_ggplot_group, ncol = 4)
  ggsave(file.path("results", "pca", paste0(object_parse, "_", condition, "_", cluster, ".pdf")), width = width, height = height,
         plot = pca_plots)
}
