#' Create hierarchical dataframe fro sunburst plotting for hierarchical visualization 
#'
#' @author Cindy van Velthoven
#'
#' @param cl.df level annotations
#' @param levels hierarchical levels to be visualized
#' @param rootname the name to be given to the center of the plot. Default is "total"
#'
#' 
#' @usage hierDF <- as.hierDF(cl.df, levels = c("class_id_label", "subclass_id_label","supertype_id_label"),rootname="wmb")
#' 
#' @return Dataframe with hierarchical representation of 

as.hierDF <- function(cl.df, 
                      levels=c("class_id_label","subclass_id_label","supertype_id_label"),
                      rootname="all") {
  
  sub.df <- cl.df[,levels]
  
  ## set root (center)
  sub.df <- data.frame(root = rootname, sub.df)
  
  
  ## add the first label set == all labels == root
  hierarchyList <- list()
  
  hierDF <- data.frame(parent = character(),
                       child = character())
  
  for(i in seq_along(levels)){
    print(i)
    
    currentCols <- c("root", levels)
    
    parentCols <- currentCols[i]
    idCols <- currentCols[(i+1)]
    
    
    if(is.numeric(sub.df[[idCols]])){
      sub.df[[idCols]] <- as.character(as.numeric(sub.df[[idCols]]))
      maxchar <- max(nchar(sub.df[[idCols]]))
      sub.df[[idCols]] <- stringr::str_pad(sub.df[[idCols]] ,maxchar, pad = "0")
    }
    
    currentDF <- sub.df %>%
      group_by(.data[[parentCols]],
               .data[[idCols]]) %>% 
      summarise(values=n()) %>%
      select(.data[[parentCols]],
             .data[[idCols]])
    
    colnames(currentDF) <- c("parent","child")
    hierDF <- rbind(hierDF,currentDF)
  }
  
  
  return(hierDF)
}


# library(ggraph)
# library(igraph)
# library(tidyverse)
# theme_set(theme_void())


# hierDF <- as.hierDF(cl.df, levels = c("class_id_label", "subclass_id_label"))
# 
# # Create a graph object
# graph <- graph_from_data_frame( hierDF)
# 
# ggraph(graph, layout = 'dendrogram', circular = FALSE) + 
#   geom_edge_diagonal()