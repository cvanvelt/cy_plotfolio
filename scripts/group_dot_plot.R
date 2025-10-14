#' Dot-plot Heatmap plots of group summary statistics
#' 
#' @param data A data frame containing gene expression values. The first column should be sample_name
#' @param anno Sample annotations. The first column should be sample_name, and each annotation should have \_id, \_label, and \_color columns
#' @param genes A character vector containing gene symbols to be plotted. 
#' @param grouping A character vector specifying the desc base that should be used to group cells
#' @param group_order Optional: Explicit specification of group order by supplying a vector of group_ids.
#' @param fill_stat The statistic to apply to each group for use as dot fill color. Default = "median". Options are: 
#' \itemize{
#'   \item "median"
#'   \item "mean"
#'   \item "tmean" (25\% trimmed mean)
#'   \item "nzmean" (mean of non-zero values)
#'   \item "nzmedian" (median of non-zero values)
#'   \item "prop_gt0" (proportion of samples > 0)
#'   \item "prop_gt1" (proportion of samples > 1)
#'   \item "min"
#'   \item "max"
#'   }
#' @param size_stat The statistic to apply to each group for scaling dot size. Same options as fill_stat. Default = "prop_gt0".
#' @param max_size Maximum size of dots, in pts.
#' @param log_scale Logical , determines if data is log scaled before plotting. Default = FALSE.
#' @param normalize_rows Logical, whether or not to rescale data within each row of the plot. Default = FALSE.
#' @param font_size numeric object, the font size (in pts) used to make the plot.
#' @param label_height numeric object, Percent of the plot height that should be used for the labels (0 to 100). Default is 25.
#' @param show_counts Logical, whether or not to display sample counts at the top of labels. Default = TRUE.
#' @param rotate_counts Logical, whether or not to rotate sample counts by 90 degrees. Default = FALSE.
#' @param max_width numeric object, percent of plot width that should be used for maximum expression values (0 to 100). Default is 10.
#' @param return_type What values to return - can be "plot", "data", or "both". Default is "plot".
#' 
#' @return a ggplot2 plot object
#'
group_dot_plot <- function(data,
                           anno,
                           genes,
                           grouping,
                           group_order = NULL,
                           fill_stat = "median",
                           size_stat = "prop_gt0",
                           max_size = 10,
                           log_scale = TRUE,
                           normalize_rows = FALSE,
                           colorset = NULL,
                           font_size = 7, 
                           label_height = 25,
                           show_counts = TRUE, 
                           rotate_counts = FALSE,
                           max_width = 10,
                           return_type = "plot",
                           show_legend = F,
                           rel.legend.width = 0.15 
                           ) {
  # Reverse so that genes go from top to bottom
  # instead of bottom to top.
  genes <- rev(genes)
  
  group_cols <- group_columns(grouping)
  
  # Filter data to genes and samples in anno
  gene_data <- filter_gene_data(data, 
                                genes, 
                                anno, 
                                group_cols,
                                group_order, 
                                "sample_name")
  
  # Filter annotations if group_order is provided
  if(!is.null(group_order)) {
    anno <- anno[anno[[group_cols$id]] %in% group_order,]
  }
  
  gene_fill_stats <- group_stats(gene_data,
                                 value_cols = genes,
                                 anno = anno,
                                 grouping = group_cols$label,
                                 stat = fill_stat)
  
  # Get maximum values for each gene before rescaling to plot space.
  max_vals_unscaled <- max_gene_vals(gene_fill_stats, genes)
  
  gene_size_stats <- group_stats(gene_data,
                                 value_cols = genes,
                                 anno = anno,
                                 grouping = group_cols$label,
                                 stat = size_stat)
  
  if(log_scale) {
    gene_fill_stats <- scale_gene_data(gene_fill_stats, genes, scale_type = "log10")
  }
  
  # Convert the data values to heatmap colors
  gene_fill_data <- data_df_to_colors(gene_fill_stats,
                                      value_cols = genes,
                                      per_col = normalize_rows,
                                      colorset = colorset)
  
  names(gene_fill_data)[match(genes, names(gene_fill_data))] <- paste0(genes, "_fill")
  names(gene_size_stats)[match(genes, names(gene_size_stats))] <- paste0(genes, "_size")
  
  # Left-join data to anno. This will ensure that data is filtered for the cells provided in anno
  plot_anno <- anno %>%
    select(one_of(group_cols$id, group_cols$label, group_cols$color)) %>%
    unique()
  
  group_counts <- anno %>%
    group_by_(group_cols$id) %>%
    summarise(group_n = n())
  
  plot_data <- plot_anno %>%
    left_join(gene_fill_data, by = group_cols$label) %>%
    left_join(gene_size_stats, by = group_cols$label) %>%
    left_join(group_counts, by = group_cols$id)
  
  # Add x-positions for each group
  plot_data <- add_group_xpos(plot_data,
                              group_cols = group_cols,
                              group_order = group_order)
  
  # Compute basic count stats that are used downstream
  # n_stats$genes, n_stats$samples, and n_stats$groups
  n_stats <- get_n_stats(plot_data, group_cols, genes)
  
  header_labels <-build_header_labels(data = plot_data, 
                                      grouping = grouping,
                                      group_order = group_order,
                                      ymin = n_stats$genes + 1, 
                                      label_height = label_height, 
                                      label_type = "simple")
  
  label_y_size <- max(header_labels$ymax) - min(header_labels$ymin)
  
  group_data <- plot_data %>%
    select(xpos, group_n) %>%
    mutate(label_y = n_stats$genes + label_y_size * 0.05,
           group_n_y = max(header_labels$ymax) - 0.1 * label_y_size)
  
  # Plot setup
  p <- ggplot() +
    scale_fill_identity() +
    scale_size_area(max_size = pt2mm(max_size)) +
    scale_y_continuous("", 
                       breaks = 1:length(genes) + 0.45, 
                       labels = genes, 
                       expand = c(0, 0)) +
    scale_x_continuous("", 
                       expand = c(0, 0)) +
    theme_classic(font_size) +
    theme(axis.text = element_text(size = rel(1), face = "italic"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none"
          ) +
    geom_hline(aes(yintercept = 1:(n_stats$genes)), size = 0.2)
  
  # plot the dots for each gene
  for(i in 1:length(genes)) {
    gene <- genes[[i]]
    gene_fill <- paste0(gene, "_fill")
    gene_size <- paste0(gene, "_size")
    p <- p + 
      geom_point(data = plot_data,
                 color = "white",
                 aes_string(x = "xpos",
                            y = i + 0.5, 
                            fill = paste0("`", gene_fill, "`"), 
                            size = paste0("`", gene_size, "`")),
                 pch = 21,
                 stroke=0)
  }
  
  
  
  # Cluster labels
  p <- ggplot_header_labels(p,
                            header_labels = header_labels,
                            header_polygons = NULL,
                            font_size = font_size)
  
  ### Maximum value labels at the right edge of the plot
  max_val_dfs <- build_max_dfs(n_stats, 
                               width_stat = "groups", 
                               max_vals_unscaled, 
                               max_width)
  
  p <- ggplot_max_vals(p,
                       n_stats = n_stats,
                       width_stat = "groups",
                       max_val_dfs = max_val_dfs,
                       font_size = font_size)
  
  # Cluster counts
  if (show_counts) {
    if (rotate_counts) {
      p <- p + geom_text(data = group_data,
                         aes(x = xpos,
                             y = group_n_y, 
                             label = group_n),
                         angle = 90,
                         hjust = 1, vjust = 0.35, 
                         size = pt2mm(font_size))
    } else {
      p <- p + geom_text(data = group_data,
                         aes(x = xpos,
                             y = group_n_y, 
                             label = group_n),
                         size = pt2mm(font_size))
    }
  }
  
  if(isTRUE(show_legend)){
    
    df <- expand.grid(x = 0:5, y = 0:5) 
    df$z <- runif(nrow(df),min = 0, max = 1)
    df$z[1] <- 1
    df$z[2] <- 0
    #plot
    dot.col.legend <- ggplot(df, aes(x, y, fill = z, size=z), shape=19) + geom_point() + 
      scale_fill_gradientn(colours=colorset,
                           na.value = "transparent",
                           breaks=c(0,1),
                           labels=c("Min","Max"),
                           limits=c(0,1)) +
      scale_size_area(max_size = pt2mm(max_size),
                      breaks = c(0.05, 0.25, 0.5, 0.75, 1),
                      labels = c(5, 25, 50, 75, 100)) +
      labs(fill = "Gene expression",
           size = "Cell expressing (%)")
    
    leg =cowplot::get_legend(dot.col.legend)
    
    p = cowplot::plot_grid(p, leg, ncol =2,rel_widths = c(1,rel.legend.width),greedy = F)
    
  } 
  
  if(return_type == "plot") {
    return(p)
  } else if(return_type == "data") {
    return(list(plot_data = plot_data,
                header_labels = header_labels,
                header_polygons = header_polygons,
                max_val_dfs = max_val_dfs,
                n_stats = n_stats))
  } else if(return_type == "both") {
    return(list(plot = p,
                plot_data = plot_data,
                header_labels = header_labels,
                group_counts = group_data,
                max_val_dfs = max_val_dfs,
                n_stats = n_stats))
  }
}





# data_df_to_colors <- function (df, value_cols = NULL, colorset = NULL, scale = "linear", 
#                                per_col = FALSE, min_val = 0, max_val = NULL) 
# {
#   library(purrr)
#   if (is.null(value_cols)) {
#     value_cols <- names(df)[-1]
#   }
#   if (scale == "log2") {
#     df[[value_cols]] <- log2(df[[value_cols]])
#   }   else if (scale == "log10") {
#     df[[value_cols]] <- log10(df[[value_cols]])
#   }
#   if (is.null(max_val) & per_col == FALSE) {
#     max_val <- max(unlist(df[, value_cols]), na.rm = TRUE)
#   }
#   df[, value_cols] <- purrr::map(value_cols, function(x) {
#     print(x)
#     vals <- unlist(df[[x]])
#     if (is.null(colorset)) {
#       values_to_colors(vals, min_val = min_val, max_val = max_val)
#     }    else {
#       values_to_colors(vals, min_val = min_val, max_val = max_val, 
#                        colorset = colorset)
#     }
#   })
#   return(df)
# }
