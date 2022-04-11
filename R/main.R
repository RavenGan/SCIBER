
# input:
# 1. batches_clean: a list contains all the pre-processed matrices.
# 2. ref_index: the index of the reference batch in the object "batches_clean"
# 3. batches_meta_data: a list contains the meta data for all the batches.
#    The order should be consistent with that in "batches_clean"
# 4. top_pairs_with_ref: proportion of matched clusters.
#    Default uses the 0.05 significance level.
# 5. n_core: specify the number of cores otherwise use all the available cores.
# 6. combine: TURE returns both raw and integrated batches with all batched combined.
#    FALSE returns a list which contains all the integrated batches. The default is TRUE.

globalVariables(c("cluster_assignment" # used in the newly created metadata
                  ))

#' Batch effect removal with SCIBER
#'
#' @param batches_clean A list contains all the pre-processed matrices with dimension of n_genes*n_cells.
#' @param ref_index The index of the reference batch in the object "batches_clean"
#' @param batches_meta_data A list contains the meta data for all the batches. The order should be consistent with that in "batches_clean"
#' @param top_pairs_prop A list of proportion of matched clusters. Default uses the 0.05 significance level for all clusters.
#' @param top_genes The number of marker genes used for Fisher exact test.
#' @param n_core Specify the number of cores otherwise use all the available cores.
#' @param combine TURE returns both raw and integrated batches with all batched combined. FALSE returns a list which contains all the integrated batches. The default is TRUE.
#'
#' @return A list of both raw and integrated batches with all batched combined or a list which contains all the integrated batches.
#' @importFrom dplyr %>%
#' @export
#'
# @examples
#SCIBER_int(batches_clean, ref_index, batches_meta_data, top_pairs_with_ref, n_core = 8, combine = TRUE)
SCIBER_int <- function(batches_clean, ref_index,
                       batches_meta_data, top_pairs_prop, top_genes = 50,
                       n_core = parallel::detectCores(), combine = TRUE
                       ) {
  datasets_cluster <- obtain_clustered_data(ref_index, batches_clean, numCores = n_core)
  new_meta_data <- obtain_new_meta_data(datasets_cluster, batches_meta_data)
  obtain_cluster_type <- obtain_cluster_type(new_meta_data)
  batches_p_tstat <- obtain_tstat_pvalue(datasets_cluster, batches_clean, numCores = n_core)
  FisherExactTest <- FisherExact_test(batches_p_tstat, obtain_cluster_type, ref_index, top_genes,
                                      numCores = n_core)
  top_pairs_with_ref_summary <- obtain_top_pairs(FisherExactTest,
                                                 top_pairs_prop,
                                                 obtain_cluster_type,
                                                 ref_index,
                                                 datasets_cluster)
  cellID_from_top_pairs_summary <- obtain_cellID_from_top_pair(top_pairs_with_ref_summary,
                                                               new_meta_data,
                                                               ref_index)
  anchor_matrices_summary <- obtain_anchor_matrices(cellID_from_top_pairs_summary, batches_clean, ref_index)
  projected_datasets_summary <- obtain_proj_to_ref(anchor_matrices_summary, batches_clean, ref_index)
  projected_original_data <- obtain_projected_original_data(projected_datasets_summary, batches_clean,
                                                            ref_index, combine = combine)

  return(projected_original_data)
}



