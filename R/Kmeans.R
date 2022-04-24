#------#------#------#------#------#------#------#------#------#------
# functions in this file:
# obtain_clustered_data: K-means clustering
# obtain_new_meta_data: obtain new meta_data for each batch
#------#------#------#------#------#------#------#------#------#------

obtain_clustered_data <- function(ref_index, batches_clean, K, numCores = parapllel::detectCores()){
  numCores <- numCores
  num_bat <- length(batches_clean)

  bat_idxes <- c(1:num_bat)
  # obtain k
  if (is.null(K)){
    k <- floor(sqrt(ncol(batches_clean[[ref_index]]))) # Use the reference batch to choose k
  } else {
    k <- K
  }

  seed_ls <- c(1:15)
  clean_batch_cluster <- c()
  for (i in 1:num_bat) {
    chosen_batch <- batches_clean[[i]]
    fx <- function(seed){
      set.seed(seed)
      return(stats::kmeans(t(chosen_batch), k, nstart = 1))
    }
    cluster_results <- parallel::mclapply(seed_ls, fx, mc.cores = numCores)

    # collect totss of all seeds
    totss_all <- c()
    for (j in 1:length(seed_ls)) {
      totss_all <- rbind(totss_all, cluster_results[[j]][["totss"]])
    }
    best_seed <- which.min(totss_all)

    clean_batch_cluster[[i]] <- as.data.frame(cluster_results[[best_seed]][["cluster"]])
    colnames(clean_batch_cluster[[i]]) <- "cluster_assignment"
  }

  return(list("k" = k, "clean_batch_cluster" = clean_batch_cluster))
}

obtain_new_meta_data <- function(datasets_cluster, batches_meta_data){
  clean_batch_cluster <- datasets_cluster$clean_batch_cluster
  new_meta_data_batches <- c()
  for (i in 1:length(batches_meta_data)) {
    chosen_meta_data <- batches_meta_data[[i]]
    chosen_meta_data$cluster_assignment <- clean_batch_cluster[[i]]
    new_meta_data_batch_i <- chosen_meta_data[, c("cell_id", "cell_type", "cluster_assignment")]
    new_meta_data_batches[[i]] <- new_meta_data_batch_i
  }
  return(new_meta_data_batches)
}
