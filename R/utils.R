# util helpers

normalize_adj <- function(A) {
  Dinv <- diag(1 / sqrt(rowSums(A)))
  Dinv[!is.finite(Dinv)] <- 0
  Dinv %*% A %*% Dinv
}

one_hot <- function(f) {
  m <- diag(nlevels(f))[f, , drop = FALSE]
  colnames(m) <- paste0(attr(f, "name", exact = TRUE), "_", levels(f))
  m
}

prepare_input_matrix <- function(traits, env = NULL) {
  X_num <- traits[sapply(traits, is.numeric)]
  X_cat <- traits[sapply(traits, is.factor)]
  cat_index <- list()
  if (length(X_cat)) {
    mats <- lapply(names(X_cat), function(nm) {
      x <- X_cat[[nm]]
      attr(x, "name") <- nm
      m <- one_hot(x)
      cat_index[[nm]] <<- seq_len(ncol(m)) + ncol(X_num) + sum(sapply(cat_index, length))
      m
    })
    X_all <- cbind(X_num, do.call(cbind, mats))
  } else {
    X_all <- as.matrix(X_num)
  }
  if (!is.null(env)) X_all <- cbind(X_all, env)
  mask <- !is.na(X_all)
  X_all[is.na(X_all)] <- 0
  attr(X_all, "mask") <- mask
  attr(X_all, "cat_index") <- cat_index
  X_all
}

fill_missing_from_matrix <- function(orig, mat, cat_index) {
  out <- orig
  num_cols <- sapply(orig, is.numeric)
  out[num_cols] <- mat[, seq_len(sum(num_cols))]
  for (nm in names(cat_index)) {
    cols <- cat_index[[nm]]
    lev  <- sub(paste0("^", nm, "_"), "", colnames(mat)[cols])
    class_imputed <- factor(lev[max.col(mat[, cols])], levels = lev)
    out[[nm]][is.na(out[[nm]])] <- class_imputed[is.na(out[[nm]])]
  }
  out
}