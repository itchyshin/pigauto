## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = requireNamespace("torch", quietly = TRUE) &&
         isTRUE(try(torch::torch_is_installed(), silent = TRUE))
)

## ----birdtree-download, eval = FALSE------------------------------------------
# install.packages("megatrees")
# bird_trees <- megatrees::get_tree_bird_n100()

## ----canonical, eval = FALSE--------------------------------------------------
# library(pigauto)
# data(avonet300, trees300)
# df <- avonet300
# rownames(df) <- df$Species_Key
# df$Species_Key <- NULL
# 
# mi <- multi_impute_trees(df, trees = trees300, m_per_tree = 1L)
# # gnn = FALSE by default; the baseline responds to each tree
# 
# mass_by_tree <- vapply(mi$datasets, function(dat) dat$Mass, numeric(nrow(df)))
# apply(mass_by_tree, 1L, stats::sd)  # descriptive sensitivity, not an MI SE

## ----opt_out, eval = FALSE----------------------------------------------------
# mi_slow <- multi_impute_trees(df, trees300, m_per_tree = 1L,
#                                gnn = TRUE, share_gnn = FALSE)
# # Fits a separate GNN on each tree; runtime depends on tree size and hardware.

