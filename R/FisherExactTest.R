
#------#------#------#------#------#------#------#------#------#------
# functions in this file:
# FisherExact_test: auxiliary function
# FisherExact_test_core: Fisher's exact test
#------#------#------#------#------#------#------#------#------#------

# need to specify the reference batch
# need to specify the top_genes


FisherExact_test_core <- function(dataset1_tstat,
                                  dataset2_tstat, # ref
                                  data1_cluster_types,
                                  data2_cluster_types, # ref
                                  top_genes,
                                  numCores = numCores){
  nClusters1 <- ncol(dataset1_tstat)
  nClusters2 <- ncol(dataset2_tstat)
  rownames_dataset1 <- rownames(dataset1_tstat)
  rownames_dataset2 <- rownames(dataset2_tstat)

  intersected_genes <- intersect(rownames_dataset1, rownames_dataset2)
  nIntersected_total <- length(intersected_genes)

  new_dataset1_tstat <- dataset1_tstat[intersected_genes, ]
  new_dataset2_tstat <- dataset2_tstat[intersected_genes, ]

  fx_conv1 <- function(col_idx1){ # convert values larger than the top_genes'th to 1 each column
    chosen_col1 <- new_dataset1_tstat[,col_idx1]
    cut1 <- sort(chosen_col1)[length(chosen_col1)-top_genes+1]
    chosen_col1[chosen_col1 > cut1] <- 1
    chosen_col1[chosen_col1 != 1] <- 0
    return(as.matrix(chosen_col1))
  }
  fx_conv2 <- function(col_idx2){ # convert values larger than the top_genes'th to 1 each column
    chosen_col2 <- new_dataset2_tstat[,col_idx2]
    cut2 <- sort(chosen_col2)[length(chosen_col2)-top_genes+1]
    chosen_col2[chosen_col2 > cut2] <- 1
    chosen_col2[chosen_col2 != 1] <- 0
    return(as.matrix(chosen_col2))
  }
  col_idx1 <- c(1:ncol(new_dataset1_tstat))
  dataset1_conv_list <- parallel::mclapply(col_idx1, fx_conv1, mc.cores = numCores)
  col_idx2 <- c(1:ncol(new_dataset2_tstat))
  dataset2_conv_list <- parallel::mclapply(col_idx2, fx_conv2, mc.cores = numCores)

  dataset1_indicator <- cbind(dataset1_conv_list[[1]],
                              dataset1_conv_list[[2]])
  dataset2_indicator <- cbind(dataset2_conv_list[[1]],
                              dataset2_conv_list[[2]])
  for (i in 3:length(dataset1_conv_list)) {
    dataset1_indicator <- cbind(dataset1_indicator,
                                dataset1_conv_list[[i]])
    dataset2_indicator <- cbind(dataset2_indicator,
                                dataset2_conv_list[[i]])
  }
  dataset1_indicator <- as.matrix(dataset1_indicator)
  dataset2_indicator <- as.matrix(dataset2_indicator)
  # colSums may not be exactly equal to top_genes
  # due to same values of t_stat

  intersection_counts <- t(dataset2_indicator)%*%dataset1_indicator

  # record row and columns of NAs
  row_na_idx <- which(is.na(intersection_counts[, 1]))
  col_na_idx <- which(is.na(intersection_counts[1, ]))
  intersection_counts[is.na(intersection_counts)] <- 0 # convert NA to 0


  # construct fisher exact test
  a <- as.vector(intersection_counts)
  b <- top_genes - a
  c <- b
  d <- nIntersected_total - a - b - c

  FisherExact_mat <- cbind(a, b)%>%
    cbind( c)%>%
    cbind( d)

  fx_fisher <- function(row_idx){
    test_table <- matrix(FisherExact_mat[row_idx, ], nrow = 2, ncol = 2, byrow = TRUE)
    p <- stats::fisher.test(test_table, alternative="two.sided")$p.value #"greater"
    return(p)
  }
  row_idx <- 1:nrow(FisherExact_mat)
  FisherTest_p <- parallel::mclapply(row_idx, fx_fisher, mc.cores = numCores)

  p_values <- c()
  for (j in 1:length(FisherTest_p)) {
    p_values[j] <- FisherTest_p[[j]]
  }
  FisherExact_summary <- matrix(p_values, nrow = nClusters1, byrow = FALSE)
  FisherExact_summary[row_na_idx, ] <- 1
  FisherExact_summary[, col_na_idx] <- 1
  return(FisherExact_summary)
}

FisherExact_test <- function(batches_p_tstat, obtain_cluster_type, ref_index,
                             top_genes, numCores = parallel::detectCores()){
  batches_tstat <- batches_p_tstat[["tstat_summary"]]
  remain_batches_tstat <- batches_tstat[-ref_index]
  remain_cluster_type <- obtain_cluster_type[-ref_index]

  ref_tstat <- batches_tstat[[ref_index]]
  ref_cluster_type <- obtain_cluster_type[[ref_index]]

  FisherExact_summary <- c()

  for (i in 1:length(remain_batches_tstat)){
    chosen_batch_tstat <- remain_batches_tstat[[i]]
    chosen_batch_cluster <- remain_cluster_type[[i]]
    FisherExact_chosen_ref <- FisherExact_test_core(chosen_batch_tstat,
                                                    ref_tstat,
                                                    chosen_batch_cluster,
                                                    ref_cluster_type,
                                                    top_genes,
                                                    numCores = numCores)
    FisherExact_summary[[i]] <- FisherExact_chosen_ref
  }
  return(FisherExact_summary)
}

