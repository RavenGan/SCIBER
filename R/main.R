
globalVariables(c("cluster_assignment" # used in the newly created metadata
                  ))

#' Batch effect removal with SCIBER
#'
#' @param input_batches A list contains all the pre-processed matrices with dimension of n_genes*n_cells.
#' @param ref_index The index of the reference batch in the object "input_batches"
#' @param batches_meta_data A list contains the meta data for all the batches. The order should be consistent with that in "input_batches". Each meta data contains three columns, "cell_id", "cell_type", and "dataset". "dataset" indicates which batch the data comes from. The row names of meta data should match the column names of batch.
#' @param omega A list of proportion of matched clusters or a single value between 0 and 1 applied to all query batches.
#' @param alpha The significance level for all clusters to choose the number of matched clusters. The default is 0.05.
#' @param h_fisher The number of marker genes used for Fisher exact test.
#' @param n_core Specify the number of cores otherwise use all the available cores.
#' @param seed random seed.
#' @param k Number of clusters used for K-means. If not provided, the default is k = the square root of n_0, where n_0 is the number of cells in the reference batch.
#'
#' @return A list which contains the reference and batch-effect-corrected batches. The order is the same as that in input_batches.
#' @importFrom dplyr %>%
#' @export
#'
#' @examples
#' data(HumanDC)
#' exp <- HumanDC[["exp"]]
#' meta <- HumanDC[["metadata"]]
#' omega <- c()
#' omega[[1]] <- 0.5
#' res <- SCIBER(input_batches = exp, ref_index = 1,
#' batches_meta_data = meta, omega = omega, n_core = 1)

SCIBER <- function(input_batches,
                       ref_index = NULL,
                       batches_meta_data = NULL,
                       omega = 0.5,
                       alpha = 0.05,
                       h_fisher = 75,
                       n_core = parallel::detectCores(),
                       seed = 7,
                       k = NULL
                       ) {
  set.seed(seed)

  # Check ref_index----
  # If ref_index is not provided, use the batch with the largest number of cells as the reference.
  if (is.null(ref_index)){
    print("ref_index is not provided and use the batch with the largest number of cells as the reference.")
    n_cells <- c()
    for (i in 1:length(input_batches)) {
      n_cells <- append(n_cells, ncol(input_batches[[i]]))
    }
    ref_index <- which.max(n_cells)
  }

  # Check batches_meta_data----
  # 1. First check the dimensions
  if (!is.null(batches_meta_data)){
    # If batches_meta_data is provided, first check the number of batches
    if (length(input_batches) != length(batches_meta_data)){
      stop("The number of batches is not equal to the number of meta data")
    } else {
      # Check whether the number of cells matches
      batch_n_cells <- c()
      meta_n_cells <- c()
      for (i in 1:length(input_batches)) {
        batch_n_cells <- append(batch_n_cells, ncol(input_batches[[i]]))
        meta_n_cells <- append(meta_n_cells, nrow(batches_meta_data[[i]]))
      }
      if (!all(batch_n_cells == meta_n_cells)){
        stop("The number of cells in some bacth(es) do(es) not match the number of cells in some meta data")
      }

      # Check whether the cell IDs match.
      for (i in 1:length(input_batches)) {
        if(!all(colnames(input_batches[[i]]) == rownames(meta_n_cells[[i]]))){
          stop(paste0("The ", i, "-th batch has different cell IDs from the ", i, "-th meta data"))
        }
      }
    }
  } else if (is.null(batches_meta_data)){
    print("The meta data for each batch is not provided and hence a pseudo meta data is created for each batch.")
    # 2. If batches_meta_data is not provided, create a pseudo meta_data
    batches_meta_data <- c()
    for (i in 1:length(input_batches)) {
      batches_meta_data[[i]] <- data.frame(cell_id = colnames(input_batches[[i]]),
                                           cell_type = "cell_type",
                                           dataset = paste0("Batch", i))
    }
  }

  # Check the number of cores used
  core_avail <- parallel::detectCores()
  if (core_avail < n_core) {
    print(paste0("The available number of cores is ", core_avail, " which is smaller than the specified ", n_core, " cores. SCIBER uses ", core_avail, " cores instead."))
    n_core <- core_avail
  } else {
    print(paste0("The available number of cores is ", core_avail, ". SCIBER uses ", n_core, " to perform batch effect removal."))
  }


  # Check the top_pairs_prop
  # If omega is not provided
  if (is.null(omega)){
    if((0 < alpha) & (alpha < 1)){
      top_pairs_prop <- alpha # Use the significance level to choose top pairs.
    } else {
      stop("The provided alpha is not within 0 and 1")
    }
  } else if (is.numeric(omega)){
    # If omega is provided, use omega to choose top pairs.
    if ((0 < omega) & (omega < 1)){
    omega_ls <- c()
    for (i in 1:(length(input_batches) - 1)) { # If only a single value of omega is provided, all the query batches use the same omega
      omega_ls <- append(omega_ls, omega)
    }
    top_pairs_prop <- omega_ls
    } else {
      stop("The provided omega is not within 0 and 1")
    }
  } else { # If omega is given as a list object.
    omega_ls <- c()
    for (i in 1:length(omega)) {
      omega_ls <- append(omega_ls, omega[[i]])
    }

    if ((length(omega_ls)+1) != length(input_batches)){
      stop("The length of ", omega, " does not match the number of query batches.")
    } else if (!all(0 < omega_ls & omega_ls < 1)){
      stop("Some omega(s) is(are) not within the range (0, 1)")
    } else {
      top_pairs_prop <- omega_ls
    }
  }

  # Check whether k is provided.
  if (!is.null(k)){
    if (k%%1 == 0){
      K <- k
    } else {
      stop("The provided k is not an integer")
    }
  } else {
    K = NULL
  }




  datasets_cluster <- obtain_clustered_data(ref_index, input_batches, K, numCores = n_core)
  new_meta_data <- obtain_new_meta_data(datasets_cluster, batches_meta_data)
  obtain_cluster_type <- obtain_cluster_type(new_meta_data)
  batches_p_tstat <- obtain_tstat_pvalue(datasets_cluster, input_batches, numCores = n_core)
  FisherExactTest <- FisherExact_test(batches_p_tstat, obtain_cluster_type, ref_index, h_fisher,
                                      numCores = n_core)
  top_pairs_with_ref_summary <- obtain_top_pairs(FisherExactTest,
                                                 top_pairs_prop,
                                                 obtain_cluster_type,
                                                 ref_index,
                                                 datasets_cluster)
  cellID_from_top_pairs_summary <- obtain_cellID_from_top_pair(top_pairs_with_ref_summary,
                                                               new_meta_data,
                                                               ref_index)
  anchor_matrices_summary <- obtain_anchor_matrices(cellID_from_top_pairs_summary, input_batches, ref_index)
  projected_datasets_summary <- obtain_proj_to_ref(anchor_matrices_summary, input_batches, ref_index)
  projected_original_data <- obtain_projected_original_data(projected_datasets_summary, input_batches,
                                                            ref_index)

  return(projected_original_data)
}



