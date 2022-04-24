#------#------#------#------#------#------#------#------#------#------
# functions in this file:
# obtain_anchor_matrices_core: obtain anchor matrices between a chosen batch and refernce batch
# obtain_anchor_matrices: auxiliary function
# obtain_proj_to_ref_core: obtain a projected batch
# obtain_proj_to_ref: auxiliary function
# obtain_projected_original_data: return projected batches and original batches
#------#------#------#------#------#------#------#------#------#------

#function: Obtain top pairs anchor matrix-----------------------------------------------------------------
obtain_anchor_matrices_core <- function(dataset1_dataset2_top_pairs_cellID, new_dataset1_clean, new_dataset2_clean){
  dataset1_anchor <- as.data.frame(matrix(0, nrow = nrow(new_dataset1_clean), ncol = length(dataset1_dataset2_top_pairs_cellID)))
  dataset2_anchor <- as.data.frame(matrix(0, nrow = nrow(new_dataset2_clean), ncol = length(dataset1_dataset2_top_pairs_cellID)))

  new_dataset1_clean <- as.matrix(new_dataset1_clean)
  new_dataset2_clean <- as.matrix(new_dataset2_clean)

  for (i in (1:length(dataset1_dataset2_top_pairs_cellID))) {
    dataset1_cluster_cellID <- dataset1_dataset2_top_pairs_cellID[[i]][["ind1"]]
    dataset2_cluster_cellID <- dataset1_dataset2_top_pairs_cellID[[i]][["ind2"]]

    dataset1_cluster <- new_dataset1_clean[, dataset1_cluster_cellID] %>%
      as.matrix()
    dataset2_cluster <- new_dataset2_clean[, dataset2_cluster_cellID] %>%
      as.matrix()


    dataset1_gene_avg <- rowSums(dataset1_cluster)/ncol(dataset1_cluster)

    dataset2_gene_avg <- rowSums(dataset2_cluster)/ncol(dataset2_cluster)

    dataset1_anchor[, i] <- dataset1_gene_avg
    dataset2_anchor[, i] <- dataset2_gene_avg
  }
  colnames(dataset1_anchor) <- c(1:length(dataset1_dataset2_top_pairs_cellID))
  rownames(dataset1_anchor) <- rownames(new_dataset1_clean)

  colnames(dataset2_anchor) <- c(1:length(dataset1_dataset2_top_pairs_cellID))
  rownames(dataset2_anchor) <- rownames(new_dataset2_clean)
  return(list("dataset1_anchor" = dataset1_anchor, "dataset2_anchor" = dataset2_anchor))
}

obtain_anchor_matrices <- function(cellID_from_top_pairs_summary, batches_clean, ref_index){
  ref_batch_clean <- batches_clean[[ref_index]]
  remain_batches_clean <- batches_clean[-ref_index]

  anchor_matrices_summary <- c()
  for (i in (1:length(cellID_from_top_pairs_summary))) {
    anchor_matrices_summary[[i]] <- obtain_anchor_matrices_core(cellID_from_top_pairs_summary[[i]],
                                                                remain_batches_clean[[i]],
                                                                ref_batch_clean)
  }
  return(anchor_matrices_summary)
}



#function: Project batch to reference--------------------------------------------------------------------

obtain_proj_to_ref_core <- function(projdataset_refdataset_anchor_matrices, new_projdataset_clean){
  projdataset_anchor <- projdataset_refdataset_anchor_matrices$dataset1_anchor%>%
    as.matrix()
  refdataset_anchor <- projdataset_refdataset_anchor_matrices$dataset2_anchor%>%
    as.matrix()

  hat_matrix <- solve(t(projdataset_anchor)%*%projdataset_anchor)%*%t(projdataset_anchor)

  projected_dataset <- refdataset_anchor%*%hat_matrix%*%as.matrix(new_projdataset_clean)%>%
    as.data.frame()
  colnames(projected_dataset) <- colnames(new_projdataset_clean)
  rownames(projected_dataset) <- rownames(new_projdataset_clean)

  return(projected_dataset)
}

obtain_proj_to_ref <- function(anchor_matrices_summary, batches_clean, ref_index){
  remain_batches_clean <- batches_clean[-ref_index]
  projected_datasets_summary <- c()
  for (i in 1:length(anchor_matrices_summary)) {
    projected_datasets_summary[[i]] <- obtain_proj_to_ref_core(anchor_matrices_summary[[i]],
                                                               remain_batches_clean[[i]])
  }
  return(projected_datasets_summary)
}

obtain_projected_original_data <- function(projected_datasets_summary, batches_clean, ref_index
                                           # ,combine = TRUE
                                           ){
  # if (combine == TRUE){
  #   ref_batch_clean <- batches_clean[[ref_index]]
  #   remain_batches_clean <- batches_clean[-ref_index]
  #
  #   projected_data_1 <- projected_datasets_summary[[1]]
  #   original_data_1 <- remain_batches_clean[[1]]
  #
  #   if (length(projected_datasets_summary) <= 1){
  #     total_projected_data <- cbind(projected_data_1, ref_batch_clean)
  #     total_original_data <- cbind(original_data_1, ref_batch_clean)
  #   } else if (length(projected_datasets_summary) == 2){
  #     total_projected_data <- cbind(projected_data_1, projected_datasets_summary[[2]])%>%
  #       cbind(ref_batch_clean)
  #     total_original_data <- cbind(original_data_1, remain_batches_clean[[2]])%>%
  #       cbind(ref_batch_clean)
  #   }
  #   else if (length(projected_datasets_summary) >= 3){
  #     total_projected_data <- cbind(projected_data_1, projected_datasets_summary[[2]])
  #     total_original_data <- cbind(original_data_1, remain_batches_clean[[2]])
  #     for (i in 3:length(projected_datasets_summary)) {
  #       total_projected_data <- cbind(total_projected_data, projected_datasets_summary[[i]])
  #       total_original_data <- cbind(total_original_data, remain_batches_clean[[i]])
  #     }
  #     total_projected_data <- as.matrix(total_projected_data %>% cbind(ref_batch_clean))
  #     total_original_data <- as.matrix(total_original_data %>% cbind(ref_batch_clean))
  #   }
  #   return(list("total_projected_data" = total_projected_data, "total_original_data" = total_original_data))
  # } else if (combine == FALSE){
  #   ref_batch_clean <- batches_clean[[ref_index]]
  #   projected_data_ls <- projected_datasets_summary
  #   projected_data_ls[[(length(projected_datasets_summary)+1)]] <- as.matrix(ref_batch_clean)
  #   return(projected_data_ls)
  # }
  ref_batch_clean <- batches_clean[[ref_index]]
  projected_data_ls <- c()
  for (i in 1:length(batches_clean)) {
    if (i < ref_index){
      projected_data_ls[[i]] <- projected_datasets_summary[[i]]
    } else if (i == ref_index){
      projected_data_ls[[i]] <- as.matrix(ref_batch_clean)
    } else if (i > ref_index){
      projected_data_ls[[i]] <- projected_datasets_summary[[i-1]]
    }
  }
  return(projected_data_ls)
}


















