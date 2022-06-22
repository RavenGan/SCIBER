
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SCIBER

<!-- badges: start -->
<!-- badges: end -->

SCIBER is a simple method that outputs the batch-effect corrected
expression data in the original space/dimension. These expression data
of individual genes can be directly used for all follow-up analyses.
SCIBER has four steps; each step has a clear biological meaning, and the
algorithms used for them are k-means clustering, t-test, Fisher’s exact
test, and linear regression, respectively, all of which are easily
comprehensible

## Installation

You can install the development version of SCIBER with the following
instructions:

``` r
# install.packages("devtools")
devtools::install_github("RavenGan/SCIBER")
```

## Example

The following example uses the pre-processed Human dendritic cell
dataset \[1\] to perform batch integration.

``` r
library(SCIBER)
rm(list = ls())
set.seed(7)
data(HumanDC)
exp <- HumanDC[["exp"]]
meta <- HumanDC[["metadata"]]

# Specify the proportion for each query batch to integrate batches.
omega <- c()
omega[[1]] <- 0.6

res <- SCIBER(input_batches = exp, ref_index = 1,
batches_meta_data = meta, omega = omega, n_core = 1)
#> [1] "The available number of cores is 10. SCIBER uses 1 to perform batch effect removal."
```

## Dataset reference

1.  Villani, A. C., Satija, R., Reynolds, G., Sarkizova, S., Shekhar,
    K., Fletcher, J., … & Hacohen, N. (2017). Single-cell RNA-seq
    reveals new types of human blood dendritic cells, monocytes, and
    progenitors. Science, 356(6335), eaah4573.

<!-- What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so: -->
<!-- ```{r cars} -->
<!-- summary(cars) -->
<!-- ``` -->
<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. You could also use GitHub Actions to re-render `README.Rmd` every time you push. An example workflow can be found here: <https://github.com/r-lib/actions/tree/v1/examples>. -->
<!-- You can also embed plots, for example: -->
<!-- ```{r pressure, echo = FALSE} -->
<!-- plot(pressure) -->
<!-- ``` -->
<!-- In that case, don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN. -->
