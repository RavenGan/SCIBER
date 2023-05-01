## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=5, fig.height=7
)

## ----eval=FALSE---------------------------------------------------------------
#  install.packages('SCIBER')

## ----eval=FALSE---------------------------------------------------------------
#  # install.packages("devtools")
#  devtools::install_github("RavenGan/SCIBER")

## -----------------------------------------------------------------------------
library(SCIBER)

## -----------------------------------------------------------------------------
data("HumanDC")
exp <- HumanDC[["exp"]]
meta <- HumanDC[["metadata"]]

## -----------------------------------------------------------------------------
omega <- c()
omega[[1]] <- 0.5

ref_index <- 1
n_core <- 1

## -----------------------------------------------------------------------------
res <- SCIBER(input_batches = exp, ref_index = ref_index,
              batches_meta_data = meta, omega = omega, n_core = n_core)

## -----------------------------------------------------------------------------
library(stats)
library(Matrix)
library(uwot)

do_PCA <- function(dat, PCs){
  dat_pca_embeddings <- prcomp(t(as.matrix(dat)), scale. = F)
  dat_pca_embeddings <- dat_pca_embeddings$x
  dat_pca_embeddings <- dat_pca_embeddings[, 1:as.numeric(PCs)]

  return(dat_pca_embeddings)
}

do_umap <- function(V) {
  umap(
    X = V,
    n_threads = 6,
    n_neighbors = 30L,
    n_components = 2L,
    metric = 'cosine',
    n_epochs = NULL,
    learning_rate = 1.0,
    min_dist = 0.3,
    spread = 1.0,
    set_op_mix_ratio = 1.0,
    local_connectivity = 1L,
    repulsion_strength = 1,
    negative_sample_rate = 1,
    a = NULL,
    b = NULL,
    fast_sgd = FALSE,
    verbose = FALSE
  )
}

meta_data <- rbind(meta[[1]], meta[[2]])
rownames(meta_data) <- meta_data$cell_id

projected_dat <- cbind(res[[1]], res[[2]])

all(rownames(meta_data) == colnames(projected_dat))

SCIBER_pca <- do_PCA(projected_dat, PCs = 20)
SCIBER_umap <- do_umap(SCIBER_pca)

## -----------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(ggthemes)
library(cowplot)

obtain_plot <- function(
  umap_use,
  meta_data,
  label_name,
  palette_use = tableau_color_pal()(10),
  pt_size = 4, point_size = 0.5, pt_shape = '.',
  base_size = 12,
  do_points = TRUE,
  do_density = FALSE,
  legend_position = "top"
){
  plt_df <- umap_use %>% data.frame() %>% cbind(meta_data) %>%
    sample_frac(1L)
  plt <- plt_df %>%
    ggplot(aes_string("X1", "X2", col = label_name,fill = label_name)) +
    theme_tufte(base_size = base_size) +
    theme(panel.background = element_rect(fill = NA, color = "black")) +
    guides(color = guide_legend(override.aes = list(stroke = 1,
                                                    alpha = 1, shape = 16, size = 4)), 
           alpha = FALSE) +
    scale_color_manual(values = palette_use, guide = "none") +
    scale_fill_manual(values = palette_use, guide = "none") +
    theme(plot.title = element_text(hjust = 0.5, family = "sans"),
          legend.text = element_text(family = "sans"),
          legend.title = element_text(family = "sans"),
          legend.position= as.character(legend_position)) +
    labs(x = "UMAP 1", y = "UMAP 2")

  if (do_points)
    plt <- plt + geom_point(shape = pt_shape, size = point_size)
  if (do_density)
    plt <- plt + geom_density_2d()

  return(plt)
}


## -----------------------------------------------------------------------------
colors_cell <- tableau_color_pal("Classic 20", 
                                 direction = 1)(length(unique(meta_data$cell_type)))
colors_batch <- tableau_color_pal("Classic Green-Orange 6", 
                                  direction = 1)(length(unique(meta_data$dataset)))

## -----------------------------------------------------------------------------
SCIBER_plt1 <- obtain_plot(SCIBER_umap, meta_data, "dataset", palette_use = colors_batch,
                           pt_shape = 19, pt_size = .4, legend_position = "top")
SCIBER_plt2 <- obtain_plot(SCIBER_umap, meta_data, "cell_type", palette_use = colors_cell,
                           pt_shape = 19, pt_size = .4, legend_position = "top")


plot_grid(SCIBER_plt1, SCIBER_plt2, nrow = 2)

## -----------------------------------------------------------------------------
sessionInfo()

