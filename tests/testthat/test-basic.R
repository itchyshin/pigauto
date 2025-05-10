
test_that('placeholder imputation fills numeric NAs', {
  library(ape)
  set.seed(1)
  tree <- rtree(5); tree$tip.label <- paste0('sp',1:5)
  traits <- data.frame(mass = c(1,NA,3,NA,4))
  env <- data.frame(temp=rnorm(5))
  res <- impute_phylo(traits, tree, env, species_id=tree$tip.label)
  expect_true(all(!is.na(res$completed_data$mass)))
})
