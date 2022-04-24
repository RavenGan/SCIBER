rm(list = ls())
set.seed(7)
library(SCIBER)
input_batches <- readRDS("./dat/HumanDC_exp.rds")
meta <- readRDS("./dat/HumanDC_meta.rds")
meta_1 <- meta[[1]]
meta_2 <- meta[[2]]
meta_data <- rbind(meta_1, meta_2)
rownames(meta_data) <- meta_data$cell_id

batches_meta_data <- meta

ref_index <- 1

h_fisher <- 50

# Specified proportion for the query batch
omega <- c()
omega[[1]] <- 0.6

test <- SCIBER_int(input_batches, ref_index, batches_meta_data,
                   omega, h_fisher, n_core = parallel::detectCores())


# Remove ref_index
test2 <- SCIBER_int(input_batches = input_batches,
                    # ref_index,
                    batches_meta_data = batches_meta_data,
                   omega = omega,
                   h_fisher = h_fisher,
                   n_core = parallel::detectCores())

# Remove meta_data
test3 <- SCIBER_int(input_batches = input_batches,
                    # ref_index,
                    # batches_meta_data = batches_meta_data,
                    omega = omega,
                    h_fisher = h_fisher,
                    n_core = parallel::detectCores())

# Remove omega
test3 <- SCIBER_int(input_batches = input_batches,
                    # ref_index,
                    # batches_meta_data = batches_meta_data,
                    # omega = omega,
                    h_fisher = h_fisher,
                    n_core = parallel::detectCores())

# Check number of cores
test4 <- SCIBER_int(input_batches = input_batches,
                    # ref_index,
                    # batches_meta_data = batches_meta_data,
                    # omega = omega,
                    h_fisher = h_fisher,
                    n_core = 2)
