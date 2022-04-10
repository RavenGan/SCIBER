#------#------#------#------#------#------#------#------#------#------
# The following test considers unequal sample size with the same variance
# The formula is t = (X_bar-Y_bar)/(S_p*\sqrt(1/n_x + 1/n_y))
# where S_p = \sqrt(((n_x-1)*S_x^2 + (n_y-1)*S_y^2)/(n_x + n_y - 2))
# where S_x and S_y are the sample variance of sample X and Y.
#------#------#------#------#------#------#------#------#------#------

obtain_tstat_pvalue_core <- function(k, data_cluster, data_clean){
  nGenes <- nrow(data_clean)
  sample_size <- ncol(data_clean)

  X <- matrix(NA, nrow=nGenes, ncol=k)
  nx <- ny <- rep(NA, k)
  cluster_idx <- data_cluster$cluster_assignment
  for (i in 1 : k) {
    chosen_col <- data_clean[, cluster_idx==i]%>%as.matrix()
    if(ncol(chosen_col) == 1){
      X[, i] <- chosen_col
    } else {
      X[, i] <- rowSums(chosen_col)
    }
    nx[i] <- sum(cluster_idx==i)
    ny[i] <- sum(cluster_idx!=i)
  }
  Y <- rowSums(data_clean) - X
  Z2 <- rowSums(data_clean ^ 2)

  numerator <- scale(X, center=FALSE, scale=nx) - scale(Y, center=FALSE, scale=ny)
  denominator <- Z2 - scale(X ^ 2, center=FALSE, scale=nx) - scale(Y ^ 2, center=FALSE, scale=ny)
  # denominator <- scale(denominator, center = -(nx+ny)*0.01^2)
  denominator[which(abs(denominator)< 1e-10)] <- 0

  denominator_1 <- sqrt(denominator)
  denominator_1[which(denominator_1==0)] <- 1e10 # Assign an extreme large value when two samples are both 0



  t <- numerator/denominator_1
  t <- scale(t, center=FALSE, scale=sqrt(1 / nx + 1 / ny) / sqrt(sample_size - 2))
  rownames(t) <- rownames(data_clean)
  return(list("data_tstat" = t))
}


obtain_tstat_pvalue <- function(datasets_cluster, batches_clean, numCores = parallel::detectCores()){
  tstat_summary <- c()

  k <- datasets_cluster$k

  iteration_times <- c(1:length(batches_clean))
  numCores <- numCores

  ttest_function <- function(iteration_times){
    chosen_batch_cluster <- datasets_cluster[["clean_batch_cluster"]][[iteration_times]]
    chosen_batch_clean <- batches_clean[[iteration_times]]
    pvalue_tstat_summary <- obtain_tstat_pvalue_core(k, chosen_batch_cluster, chosen_batch_clean)
    return(pvalue_tstat_summary$data_tstat)
  }
  tstat_summary <- parallel::mclapply(iteration_times, ttest_function, mc.cores = numCores)
  return(list("tstat_summary" = tstat_summary))
}





