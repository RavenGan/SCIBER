#------#------#------#------#------#------#------#------#------#------
# functions in this file:
# obtain_cluster_type: obtain cluster types based on the maximum voting
#------#------#------#------#------#------#------#------#------#------

obtain_cluster_type <- function(new_meta_data){
  cluster_type_summary <- c()
  for (i in 1:length(new_meta_data)){
    new_meta_data_i <- new_meta_data[[i]]%>%
      as.matrix()
    Cluster_CellType_summary <- table(new_meta_data_i[,2], new_meta_data_i[,3])%>%
      as.data.frame.matrix()
    nCluster <- ncol(Cluster_CellType_summary)
    Cluster_type <- c()
    for (j in 1:nCluster) {
      max_index <- which.max(Cluster_CellType_summary[, j])
      Cluster_type[j] <- rownames(Cluster_CellType_summary)[max_index]
    }
    Cluster_CellType_summary <- rbind(Cluster_CellType_summary, Cluster_type)
    row.names(Cluster_CellType_summary)[nrow(Cluster_CellType_summary)] <- "Cluster_type"
    cluster_type_summary[[i]] <- Cluster_CellType_summary
  }
  return(cluster_type_summary)
}

