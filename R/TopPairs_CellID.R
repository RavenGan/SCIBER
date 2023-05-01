
#------#------#------#------#------#------#------#------#------#------
# functions in this file:
# obtain_top_pairs: auxiliary function
# obtain_top_pairs_core: obtain top pairs
# obtain_cellID_from_top_pair: auxiliary function
# obtain_cellID_from_top_pair_core: obtain top-pairs' cellID
#------#------#------#------#------#------#------#------#------#------

#function: Choose top pairs---------------------------------------------------------------------------------------
# nRow means the number of clusters in dataset 2 (reference), while nCol means the number of clusters in dataset 1
obtain_top_pairs_core <- function(p_value_summary, nPairs, obtain_cluster_type_dataset1, obtain_cluster_type_dataset2){
  nCounts <- 1
  chosen_top_pairs <- data.frame(matrix(NA, nrow = nPairs, ncol = 2))
  colnames(chosen_top_pairs) <- c("row", "col")
  p_value_col <- as.data.frame(matrix(NA, nrow = nPairs, ncol = 1))
  colnames(p_value_col) <- c("p_values")
  while (nCounts <= nPairs) {
    minimal_location <- which(p_value_summary == min(p_value_summary, na.rm = TRUE), arr.ind = TRUE)%>%
      as.data.frame()
    chosen_top_pairs[nCounts, ] <- minimal_location[1, ]
    minimal_row <- minimal_location[1, 1]
    minimal_col <- minimal_location[1, 2]
    p_value_col[nCounts, 1] <- p_value_summary[minimal_row, minimal_col]
    p_value_summary[, minimal_col] <- 1 #Here choosing any value >= 1 should be fine
    nCounts <- nCounts + 1
  }
  chosen_top_pairs <- cbind(chosen_top_pairs, p_value_col) #Row stands for dataset2; column stands for dataset1

  row_type <- as.data.frame(matrix(NA, nrow = nrow(chosen_top_pairs), ncol = 1))
  colnames(row_type) <- c("row_cluster_types")
  col_type <- as.data.frame(matrix(NA, nrow = nrow(chosen_top_pairs), ncol = 1))
  colnames(col_type) <- c("col_cluster_types")
  matching_record <- as.data.frame(matrix(0, nrow = nrow(chosen_top_pairs), ncol = 1))
  colnames(matching_record) <- c("matching_decision")
  for (i in 1:nrow(chosen_top_pairs)) {
    row_index <- chosen_top_pairs[i, 1]
    col_index <- chosen_top_pairs[i, 2]

    row_type[i, 1] <- obtain_cluster_type_dataset2[nrow(obtain_cluster_type_dataset2),
                                                   which(as.integer(colnames(obtain_cluster_type_dataset2)) == row_index)]
    col_type[i, 1] <- obtain_cluster_type_dataset1[nrow(obtain_cluster_type_dataset1),
                                                   which(as.integer(colnames(obtain_cluster_type_dataset1)) == col_index)]
    if (row_type[i, 1] == col_type[i, 1]){
      matching_record[i, 1] <- 1
    } else {
      matching_record[i, 1] <- 0
    }
  }
  chosen_top_pairs <- chosen_top_pairs %>%
    cbind(row_type) %>%
    cbind(col_type) %>%
    cbind(matching_record)
  return(chosen_top_pairs)
}

# nRow means the number of clusters in dataset 2 (reference), while nCol means the number of clusters in dataset 1
obtain_top_pairs_core2 <- function(p_value_summary, sig_level, obtain_cluster_type_dataset1, obtain_cluster_type_dataset2){
  row_idx <- c()
  col_idx <- c()
  p_value <- c()

  while (min(p_value_summary, na.rm = TRUE) <= sig_level) {
    minimal_location <- which(p_value_summary == min(p_value_summary, na.rm = TRUE), arr.ind = TRUE)%>%
      as.data.frame()
    minimal_row <- minimal_location[1, 1]
    minimal_col <- minimal_location[1, 2]
    row_idx <- append(row_idx, minimal_location[1, 1])
    col_idx <- append(col_idx, minimal_location[1, 2])
    p_value <- append(p_value, p_value_summary[minimal_row, minimal_col])
    p_value_summary[, minimal_col] <- 1 #Here choosing any value >= 1 should be fine
  }
  chosen_top_pairs <- data.frame(row = row_idx,
                                 col = col_idx,
                                 p_values = p_value)
  row_type <- as.data.frame(matrix(NA, nrow = nrow(chosen_top_pairs), ncol = 1))
  colnames(row_type) <- c("row_cluster_types")
  col_type <- as.data.frame(matrix(NA, nrow = nrow(chosen_top_pairs), ncol = 1))
  colnames(col_type) <- c("col_cluster_types")
  matching_record <- as.data.frame(matrix(0, nrow = nrow(chosen_top_pairs), ncol = 1))
  colnames(matching_record) <- c("matching_decision")
  for (i in 1:nrow(chosen_top_pairs)) {
    row_index <- chosen_top_pairs[i, 1]
    col_index <- chosen_top_pairs[i, 2]

    row_type[i, 1] <- obtain_cluster_type_dataset2[nrow(obtain_cluster_type_dataset2),
                                                   which(as.integer(colnames(obtain_cluster_type_dataset2)) == row_index)]
    col_type[i, 1] <- obtain_cluster_type_dataset1[nrow(obtain_cluster_type_dataset1),
                                                   which(as.integer(colnames(obtain_cluster_type_dataset1)) == col_index)]
    if (row_type[i, 1] == col_type[i, 1]){
      matching_record[i, 1] <- 1
    } else {
      matching_record[i, 1] <- 0
    }
  }
  chosen_top_pairs <- chosen_top_pairs %>%
    cbind(row_type) %>%
    cbind(col_type) %>%
    cbind(matching_record)
  return(chosen_top_pairs)
}



obtain_top_pairs <- function(FisherExactTest, top_pairs_prop, obtain_cluster_type, ref_index, datasets_cluster){
  remain_cluster_type <- obtain_cluster_type[-ref_index]
  ref_cluster_type <- obtain_cluster_type[[ref_index]]
  k <- datasets_cluster$k
  top_pairs_with_ref_summary <- c()
  if(!is.numeric(top_pairs_prop)) {
    for (i in 1:length(FisherExactTest)) {
      nPairs <- ceiling(top_pairs_prop[[i]]*k)
      chosen_ref_top_pairs <- obtain_top_pairs_core(FisherExactTest[[i]],
                                                    nPairs,
                                                    remain_cluster_type[[i]],
                                                    ref_cluster_type)
      top_pairs_with_ref_summary[[i]] <- chosen_ref_top_pairs
    }
    return(top_pairs_with_ref_summary)
  } else {
    for (i in 1:length(FisherExactTest)) {
      chosen_ref_top_pairs <- obtain_top_pairs_core2(FisherExactTest[[i]],
                                                    top_pairs_prop,
                                                    remain_cluster_type[[i]],
                                                    ref_cluster_type)
      top_pairs_with_ref_summary[[i]] <- chosen_ref_top_pairs
    }
    return(top_pairs_with_ref_summary)
  }
}



#function: Obtain required datasets, a multilist containing the cell ID from the top pairs.----------------
obtain_cellID_from_top_pair_core <- function(dataset1_dataset2_top, new_meta_data_dataset1, new_meta_data_dataset2){
  #Row stands for dataset2; column stands for dataset1
  mlist <- list()
  for (i in 1:nrow(dataset1_dataset2_top)){
    row_index <- dataset1_dataset2_top[i, 1]
    col_index <- dataset1_dataset2_top[i, 2]
    dataset1_cluster <- new_meta_data_dataset1[which(new_meta_data_dataset1$cluster_assignment == col_index), , drop = FALSE]
    dataset2_cluster <- new_meta_data_dataset2[which(new_meta_data_dataset2$cluster_assignment == row_index), , drop = FALSE]

    rownames(dataset1_cluster) <- dataset1_cluster$cell_id
    rownames(dataset2_cluster) <- dataset2_cluster$cell_id

    mlist[[i]] <- list(ind1 = rownames(dataset1_cluster), ind2 = rownames(dataset2_cluster))
  }
  return(mlist)
}


obtain_cellID_from_top_pair <- function(top_pairs_with_ref_summary, new_meta_data, ref_index){
  ref_meta_data <- new_meta_data[[ref_index]]
  remain_meta_data <- new_meta_data[-ref_index]

  cellID_from_top_pair_summary <- c()
  for (i in 1:length(top_pairs_with_ref_summary)){
    cellID_from_top_pair <- obtain_cellID_from_top_pair_core(top_pairs_with_ref_summary[[i]],
                                                             remain_meta_data[[i]],
                                                             ref_meta_data)
    cellID_from_top_pair_summary[[i]] <- cellID_from_top_pair
  }
  return(cellID_from_top_pair_summary)
}




