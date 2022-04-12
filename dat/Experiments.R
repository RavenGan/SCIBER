rm(list = ls())
set.seed(7)
batches_clean <- readRDS("./dat/HumanDC_exp.rds")
meta <- readRDS("./dat/HumanDC_meta.rds")
meta_1 <- meta[[1]]
meta_2 <- meta[[2]]
meta_data <- rbind(meta_1, meta_2)
rownames(meta_data) <- meta_data$cell_id

batches_meta_data <- meta

ref_index <- 1

top_genes <- 50

# Specified proportion for the query batch
top_pairs_prop <- c()
top_pairs_prop[[1]] <- 0.6

test <- SCIBER_int(batches_clean, ref_index, batches_meta_data,
                   top_pairs_prop, top_genes, n_core = parallel::detectCores(),
                   combine = TRUE)


# Use significance level
top_pairs_prop <- 0.05
test2 <- SCIBER_int(batches_clean, ref_index, batches_meta_data,
                   top_pairs_prop, top_genes, n_core = parallel::detectCores(),
                   combine = TRUE)


m <- c(1, 2, 3)
n <- c(1, 2)
if(length(m) != length(n)){
  stop(paste0("Stop here ", 1, " LOL"))
} else {
  print(m + n)
}

