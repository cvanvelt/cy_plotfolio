
add_group_xpos <- function (data, group_cols, group_order = NULL) 
{
  if (is.null(group_order)) {
    group_order_df <- data %>% dplyr::select(one_of(group_cols$id)) %>% 
      unique() %>% dplyr::arrange_(group_cols$id) %>% dplyr::mutate(xpos = 0:n())
    data <- data %>% dplyr::left_join(group_order_df, by = group_cols$id)
  }
  else {
    group_order_df <- data.frame(group = group_order) %>% 
      dplyr::mutate(xpos = seq(0,n()-1))
    names(group_order_df)[1] <- group_cols$id
    data <- data %>% dplyr::left_join(group_order_df, by = group_cols$id)
  }
  data
}
