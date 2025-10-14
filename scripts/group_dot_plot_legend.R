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
  
  
  gf <- gene_fill_data
  colnames(gf) <- gsub("_fill","", colnames(gf))
  gf <- pivot_longer(gf,
                     cols = -group_cols$label,
                     names_to = "genes", 
                     values_to = "fill")
  
  gs <- gene_size_stats
  colnames(gs) <- gsub("_size","", colnames(gs))
  gs <- pivot_longer(gs,
                     cols = -group_cols$label,
                     names_to = "genes", 
                     values_to = "size")
  
  plot_dat <- plot_anno %>%
    right_join(gf, by = group_cols$label,relationship = "many-to-many") %>% 
    left_join(gs) 
  
  plot_dat <- plot_dat %>% left_join(select(plot_data,group_cols$label, xpos ))
  gid <- as.data.frame(genes)
  gid$id <- 1:nrow(gid)
  plot_dat$i <- gid$id[match(plot_dat$genes, gid$genes)]
  
  # Plot setup
  p <- ggplot() +
    scale_fill_identity() +
    scale_size_area(max_size = pt2mm(max_size),
                    breaks = c(0.05, 0.25, 0.5, 0.75, 1),
                    labels = c(5, 25, 50, 75, 100)) +
    scale_y_continuous("", 
                       breaks = 1:length(genes) + 0.45, 
                       labels = genes, 
                       expand = c(0, 0)) +
    scale_x_continuous("", 
                       expand = c(0, 0)) +
    theme_classic(font_size) +
    geom_point(data = plot_dat,
               color = "white",
               aes(x = xpos,
                   y = i + 0.5, 
                   fill = fill, 
                   size = size),
               pch = 21,
               stroke=0) +
    theme(#legend.position = "none",
      axis.text = element_text(size = rel(1), face = "italic"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank() ) +
    geom_hline(aes(yintercept = 1:(n_stats$genes)), size = 0.2) +
    guides(size = guide_legend(override.aes=list(shape = 19, color="black"), keyheight = 0.75)) +
    geom_line(data = plot_dat, aes(x=0,y=i, color=i), size=0) +
    scale_color_gradientn(colours=colorset,
                          na.value = "transparent",
                          breaks=c(0,1),
                          labels=c("Min","Max"),
                          limits=c(0,1),
                          guide = guide_colorbar(barwidth = 0.75,
                                                 barheight = 4)) +
    labs(color = "Gene expression",
         size = "Cell expressing (%)")
  
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
    p <- p
  } else{
    p <- p + theme(legend.position = "none")
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



addSmallLegend <- function(myPlot, pointSize = 0.75, textSize = 6, spaceLegend = 0.5, ncol=1) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize)),
           fill=guide_legend(ncol=ncol)) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}



group_split_dot_plot <- function (data, anno, genes, grouping, group_order = NULL, split_by = NULL, 
    split_order = NULL, max_size = 10, fill_stat = "median", 
    size_stat = "prop_gt0", log_scale = TRUE, normalize_rows = FALSE, 
    colorset = NULL, font_size = 7, label_height = 25, show_counts = TRUE, 
    rotate_counts = FALSE, max_width = 10, row_padding = 0.1, offset=-2.5, missing_size=0.05, show_legend = F,rel.legend.width = 0.15 ) 
{
    library(dplyr)
    library(ggplot2)
    library(purrr)
    genes <- rev(genes)
    group_id <- paste0(grouping, "_id")
    group_label <- paste0(grouping, "_label")
    group_color <- paste0(grouping, "_color")
    group_anno <- anno %>% select(one_of(group_id, group_label, 
                                        group_color)) %>% unique()
    if (!is.null(split_by)) {
        split_id <- paste0(split_by, "_id")
        anno[[split_id]] = as.factor(anno[[split_id]])
        split_label <- paste0(split_by, "_label")
        split_color <- paste0(split_by, "_color")
        split_ids <- unique(anno[[split_id]])
        split_anno <- anno %>% select(one_of(split_id, split_label, 
                                             split_color)) %>% unique()        
        split_samples <- map(split_ids, function(x) anno$sample_name[anno[[split_id]] == 
            x])
        names(split_samples) <- split_ids
    }
    gene_data <- data[, c("sample_name", genes)]
    gene_data <- gene_data[match(anno$sample_name, data$sample_name), 
        ]
    if (is.null(split_by)) {
        gene_fill_stats <- group_stats(gene_data, value_cols = genes, 
            anno = anno, grouping = group_id, stat = fill_stat)
    }
    else {
        gene_fill_stats <- group_stats(gene_data, value_cols = genes, 
            anno = anno, grouping = c(group_id, split_id), stat = fill_stat)
    }
    if (is.null(split_by)) {
        max_vals <- map_dbl(genes, function(x) {
            max(gene_fill_stats[[x]])
        })
        names(max_vals) <- genes
        gene_size_stats <- group_stats(gene_data, value_cols = genes, 
            anno = anno, grouping = group_id, stat = size_stat)
    }
    else {
        max_vals <- map_dfr(split_ids, function(x) {
            split_data <- gene_fill_stats[gene_fill_stats[[split_id]] == 
                x, ]
            max_vals <- map_dfc(genes, function(x) {
                max(split_data[[x]])
            })
            names(max_vals) <- genes
            max_vals
        })
        max_vals <- cbind(split_col = split_ids, max_vals)
        names(max_vals)[1] <- split_id
        gene_size_stats <- group_stats(gene_data, value_cols = genes, 
            anno = anno, grouping = c(group_id, split_id), stat = size_stat)
    }
    if (log_scale) {
        gene_fill_stats[, genes] <- log10(gene_fill_stats[, genes] + 
            1)
    }
    if (is.null(split_by)) {
        gene_fill_data <- data_df_to_colors(gene_fill_stats, 
            value_cols = genes, per_col = normalize_rows, colorset = colorset)
    }
    else {
        gene_fill_data <- map_dfr(split_ids, function(x) {
            split_fill_stats <- gene_fill_stats[gene_fill_stats[[split_id]] == 
                x, ]
            data_df_to_colors(split_fill_stats, value_cols = genes, 
                per_col = normalize_rows, colorset = colorset)
        })
    }
    names(gene_fill_data)[match(genes, names(gene_fill_data))] <- paste0(genes, 
        "_fill")
    names(gene_size_stats)[match(genes, names(gene_size_stats))] <- paste0(genes, 
        "_size")
    if (is.null(split_by)) {
        plot_anno <- group_anno
        group_n <- anno %>% group_by_(group_id) %>% summarise(group_n = n())
        plot_data <- plot_anno %>% left_join(gene_fill_data, 
            by = group_id) %>% left_join(gene_size_stats, by = group_id) %>% 
            left_join(group_n, by = c(group_id))
    }
    else {
        group_n <- anno %>% group_by_at(c(group_id, split_id),.drop=FALSE) %>% 
          summarise(group_n = n())
        plot_anno <- group_n %>% left_join(group_anno)  %>% left_join(split_anno)
        plot_data <- plot_anno %>% left_join(gene_fill_data, 
            by = c(group_id, split_id)) %>% left_join(gene_size_stats, 
            by = c(group_id, split_id)) 
    }
    if (!is.null(group_order)) {
        group_order_df<- data.frame(group = group_order) %>% 
            mutate(xpos = 1:n())
        names(group_order_df)[1] <- group_id
        plot_data <- plot_data %>% filter_(paste0(group_id, " %in% group_order")) %>% 
            left_join(group_order_df, by = group_id)
    }
    else {
      group_order_df <- group_anno %>% arrange_(group_id) %>% select(one_of(group_id))%>% mutate(xpos = 1:n())
      plot_data <- plot_data %>% left_join(group_order_df, 
                                           by = group_id)
    }
    if (is.null(split_by)) {
        plot_data <- plot_data %>% mutate(ypos_adj = 0)
    }
    else {
        adj_vals <- seq(1 - row_padding, row_padding, length.out = length(split_ids) * 
            2 + 1)[1:length(split_ids) * 2]
        if (!is.null(split_order)) {
            split_order <- split_order[split_order %in% anno[[split_id]]]
            split_order_df <- data.frame(split_val = split_order) %>% 
                mutate(ypos_adj = adj_vals)
            names(split_order_df)[1] <- split_id
            plot_data <- plot_data %>% filter_(paste0(split_id, 
                " %in% split_order")) %>% left_join(split_order_df, 
                by = split_id)
        }
        else {
          split_order_df<- split_anno %>% arrange_(split_id) %>% select(one_of(split_id)) %>% mutate(ypos_adj = adj_vals)
          plot_data <- plot_data %>% left_join(split_order_df, by = split_id)
                                               
        }
    }
    n_genes <- length(genes)
    n_groups <- length(unique(plot_data[[group_id]]))
    n_samples <- nrow(data)
    header_labels <- build_header_labels(data = plot_data, grouping = grouping, 
        group_order = group_order, ymin = n_genes + 1, label_height = label_height, 
        label_type = "simple")
    if (is.null(split_by)) {
        max_labels <- data.frame(x = (n_groups + 0.5) * 1.01, 
            y = 1:n_genes + 0.5, label = sci_label(max_vals))
    }
    else {
     max_labels = melt(max_vals)
     colnames(max_labels)[2] = "gene"
     max_labels$x = (n_groups + 0.5) * 1.01
     max_labels = max_labels %>% left_join(split_order_df)
     max_labels$y =  as.integer(max_labels[[2]]) + max_labels$ypos_adj
     max_labels$label = sci_label(max_labels$value)
   }
    split_gene_pos = max_labels %>% select(split_id, gene, y)
    
    max_header <- data.frame(x = (n_groups + 0.7) * 1.01, y = n_genes + 
        1, label = paste0("Max ", fill_stat, " value"))
    max_width <- n_groups * (max_width/100)/(1 - max_width/100)
    label_y_size <- max(header_labels$ymax) - min(header_labels$ymin)

    p <- ggplot() + 
    scale_fill_identity() + 
    scale_size_area(max_size = pt2mm(max_size),
                    breaks = c(0.05, 0.25, 0.5, 0.75, 1),
                    labels = c(5, 25, 50, 75, 100)) + 
    scale_y_continuous("", breaks = 1:length(genes) + 0.5, 
                labels = genes, expand = c(0, 0), limits= c(0, max(header_labels$ymax))) +
    scale_x_continuous("", limits = c(offset, n_groups + 0.5 + max_width), expand = c(0, 0)) +
    theme_classic(font_size) + 
    theme(axis.text.y = element_text(size = rel(2), 
        face = "italic"), axis.text.x = element_blank(), axis.ticks = element_blank(), 
        legend.position = "none", axis.line.y = element_blank()) + 
        geom_hline(aes(yintercept = 1:(n_genes)), size = 0.2)

    if (!is.null(split_by)) {
        split_line_adj <- seq(1 - row_padding, row_padding, length.out = length(split_ids) + 
            1)
        split_line_adj <- split_line_adj[-length(split_line_adj)]
        split_line_adj <- split_line_adj[-1]
        #p <- p + geom_hline(aes(yintercept = 1:n_genes + split_line_adj), 
                                        #    size = 0.2, linetype = "dashed")
        split_label_pos <- data.frame(gene = rep(1:n_genes, each = length(split_ids)), 
            split_col = rep(split_ids, n_genes))
        names(split_label_pos)[2] <- split_id
        split_label_pos <- split_label_pos %>% left_join(split_order_df, 
            by = split_id) %>% left_join(split_anno) %>% mutate(ypos = gene + ypos_adj)
        p <- p + geom_text(data = split_label_pos, aes_string(x = "0", 
            y = "ypos", label = split_label), hjust = 1, vjust = 0.3, 
            size = pt2mm(font_size))
      }

      
    for (i in 1:length(genes)) {
        gene <- genes[[i]]
        gene_fill <- paste0(gene, "_fill")
        gene_size <- paste0(gene, "_size")
        ypos <- paste0(i, " + ypos_adj")
        p <- p + geom_point(data = plot_data, colour="white",aes_string(x = "xpos", 
            y = ypos, fill = gene_fill, size = gene_size), pch = 21)
        missing_data = plot_data %>% filter(group_n==0)
        if(nrow(missing_data)>0){
          missing_data$label = "x"
          p <- p + geom_point(data = missing_data, colour="black",
                             aes_string(x = "xpos", y = ypos, size =missing_size),pch=4)
        }
      }
    p <- p + geom_rect(data = header_labels, aes(xmin = xmin, ymin = ymin + 0.05 ,
                         xmax = xmax, ymax = ymin + 0.15, fill = color)) +
             geom_text(data = header_labels, aes(x = (xmin + xmax)/2, 
                         y = ymin + 0.20, label = label), angle = 90, vjust = 0.35, 
                       hjust = 0, size = pt2mm(font_size)) +
             geom_text(data = max_header, 
                       aes(x = x, y = y, label = label), hjust = 0, 
                       vjust = 0.35, size = pt2mm(font_size)) + geom_text(data = max_labels, 
                                       aes(x = x, y = y, label = label), hjust = 0, vjust = 0.35, 
                                       size = pt2mm(font_size), parse = TRUE)
    if (show_counts) {
      group_data <- plot_data %>% select(xpos, group_n) %>%
        mutate(label_y = n_genes + 
               label_y_size * 0.05, group_n_y = max(header_labels$ymax) - 
               0.1 * label_y_size)
      if (rotate_counts) {
        p <- p + geom_text(data = group_data, aes(x = xpos, 
                             y = group_n_y, label = group_n), angle = 90, 
                           hjust = 1, vjust = 0.35, size = pt2mm(font_size))
      }
      else {
        p <- p + geom_text(data = group_data, aes(x = xpos, 
                             y = group_n_y, label = group_n), size = pt2mm(font_size))
      }
    }
    
    return(p)
}