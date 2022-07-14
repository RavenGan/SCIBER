library(SCIBER)
data("HumanDC")

exp <- HumanDC[["exp"]]
meta <- HumanDC[["metadata"]]
omega <- c()
omega[[1]] <- 0.5

res <- SCIBER(input_batches = exp, ref_index = 1,
                  batches_meta_data = meta, omega = omega, n_core = 1)

test_that('Dimensions of the input and output data match', {
  expect_equal(dim(res[[1]]), dim(exp[[1]]))
  expect_equal(dim(res[[2]]), dim(exp[[2]]))
})

test_that('There are no null values in the corrected embedding', {
  expect_true(all(!is.na(res[[1]])))
  expect_true(all(!is.na(res[[2]])))
})


test_that('Error messages work', {
  expect_error(
    SCIBER_int(input_batches = exp, ref_index = 1,
               batches_meta_data = list(meta[[1]]), omega = omega, n_core = 1)
  )
  expect_error(
    SCIBER_int(input_batches = exp, ref_index = 1,
               batches_meta_data = meta, omega = list(omega[[1]], omega[[1]]), n_core = 1)
  )
})
