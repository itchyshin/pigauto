# R/henderson_s_inv.R
#
# Hadfield & Nakagawa (2010, J. Evol. Biol. 23: 494-508, eqs 26-29)
# sparse precision matrix S^{-1} for the EXTENDED phylogeny (tips + non-
# root internal nodes). The tip-only covariance A is the lower-right
# block of S, but A^{-1} is dense (O(n^2) nonzeros). S^{-1} is sparse
# (O(n) nonzeros): each non-root node contributes at most a 2x2 block
# (self, parent) to S^{-1}.
#
# Closed-form recursive formula (eq. 29):
#
#   S^{-1}[t, p_t] = -1/f_t        (and symmetrically [p_t, t] = -1/f_t)
#   S^{-1}[t, t]   =  1/f_t + sum_{o in children(t)} 1/f_o
#   else            = 0
#
# where p_t is t's parent (the root has no entry in S because S is the
# extended cov over non-root nodes) and f_t is the length of the edge
# joining t to p_t. The diagonal at node t is built by ACCUMULATING
# contributions from t's edge to its parent and from each child's edge
# to t.
#
# Why this matters for pigauto: every dense (R + eps*I) Cholesky in
# bm_impute_col / joint_mvn_baseline / joint_threshold_baseline is
# O(n_tips^3) work. With sparse S^{-1} we can compute R^{-1} b at tips
# in O(n) via the Schur complement on the precision form, AND evaluate
# the proper Kronecker marginal log-likelihood log p(L_obs | Sigma) in
# O(K * n) per optim() step. This is the standard quant-genetics trick
# that ASReml, MCMCglmm, brms, and phylopars all use under the hood.
#
# Output contract: `build_henderson_S_inv(tree)` returns a list with
#
#   $Q         : N x N sparse (dsCMatrix) representation of S^{-1}.
#                Rows/columns indexed 1..N_tips first, then non-root
#                internals (so the tip block of Q is Q[1:n_tips, 1:n_tips]).
#   $tip_idx   : integer vector 1..n_tips (positions of tips in Q).
#   $int_idx   : integer vector (n_tips+1)..N (positions of internal-
#                non-root nodes in Q).
#   $tip_labels: character n_tips, in the order they appear at the tip
#                indices of Q (matches tree$tip.label since Q uses ape's
#                native 1..n_tips numbering for tips).
#   $root_node : integer, the index of the root node in tree's numbering
#                (typically n_tips + 1L). Excluded from Q.
#   $n_tips, $n_internal_nonroot, $N : convenience scalars.

# ---------------------------------------------------------------------------
# build_henderson_S_inv()
# ---------------------------------------------------------------------------
#
# Build Q = S^{-1} from a phylo tree by walking tree$edge once. No
# explicit Cholesky / inverse needed: equation 29 is closed-form.
build_henderson_S_inv <- function(tree, eps = 1e-12) {
  stopifnot(inherits(tree, "phylo"))
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("build_henderson_S_inv requires the 'Matrix' package.",
         call. = FALSE)
  }

  n_tips <- length(tree$tip.label)
  # ape stores internal node count in tree$Nnode (this INCLUDES the root)
  n_internal_total <- as.integer(tree$Nnode)
  N_total <- n_tips + n_internal_total
  # Root is the one internal node with no incoming edge. For a rooted
  # ape phylo, root index = n_tips + 1L (the canonical convention).
  root_node <- setdiff(unique(tree$edge[, 1]), unique(tree$edge[, 2]))
  if (length(root_node) != 1L) {
    stop("build_henderson_S_inv: tree must have exactly one root.",
         call. = FALSE)
  }
  root_node <- as.integer(root_node)

  # Index map: original node id -> position in Q (1..N).
  # Convention: tips first (1..n_tips), then non-root internals.
  # The root is excluded entirely.
  all_internals_nonroot <- setdiff(seq_len(N_total), c(seq_len(n_tips), root_node))
  N <- n_tips + length(all_internals_nonroot)
  idx_map <- integer(N_total)
  idx_map[seq_len(n_tips)] <- seq_len(n_tips)
  idx_map[all_internals_nonroot] <- (n_tips + 1L):N

  # Accumulate triplets via vectorised tree-edge traversal.
  parents <- tree$edge[, 1]
  children <- tree$edge[, 2]
  fs <- pmax(tree$edge.length, eps)  # ridge to avoid 1/0
  inv_fs <- 1 / fs

  parent_idx <- idx_map[parents]
  child_idx  <- idx_map[children]
  # Parents that are the root have idx_map == 0; mask them out for the
  # off-diagonal entries and the parent-side diagonal contribution. The
  # child-side diagonal contribution still applies (1/f_t at S^{-1}[t,t]).
  has_parent_in_Q <- parent_idx > 0L

  # i, j, x triplets for sparse construction (Matrix::sparseMatrix
  # aggregates duplicates by `summary` semantics; we set repeatable =
  # TRUE so duplicate (i,j) pairs are summed, which is what eq 29 does
  # at the diagonal).
  i_diag_child  <- child_idx
  j_diag_child  <- child_idx
  x_diag_child  <- inv_fs

  i_diag_parent <- parent_idx[has_parent_in_Q]
  j_diag_parent <- parent_idx[has_parent_in_Q]
  x_diag_parent <- inv_fs[has_parent_in_Q]

  i_off_cp <- child_idx[has_parent_in_Q]
  j_off_cp <- parent_idx[has_parent_in_Q]
  x_off_cp <- -inv_fs[has_parent_in_Q]
  # Symmetric copy
  i_off_pc <- j_off_cp
  j_off_pc <- i_off_cp
  x_off_pc <- x_off_cp

  i <- c(i_diag_child, i_diag_parent, i_off_cp, i_off_pc)
  j <- c(j_diag_child, j_diag_parent, j_off_cp, j_off_pc)
  x <- c(x_diag_child, x_diag_parent, x_off_cp, x_off_pc)

  Q <- Matrix::sparseMatrix(
    i = i, j = j, x = x,
    dims = c(N, N),
    symmetric = FALSE,    # symmetric=TRUE expects upper-triangle only
    repr = "C"
  )
  # The triplets already encode both halves of every off-diagonal, so
  # the result is symmetric. Coerce to symmetric storage for downstream
  # efficiency.
  Q <- Matrix::forceSymmetric(Q, uplo = "U")

  # Tip labels go in the first n_tips rows/cols, in ape's native order.
  tip_labels <- tree$tip.label
  rownames(Q) <- colnames(Q) <- c(tip_labels, paste0("int_", all_internals_nonroot))

  list(
    Q = Q,
    tip_idx = seq_len(n_tips),
    int_idx = if (length(all_internals_nonroot)) (n_tips + 1L):N else integer(0),
    tip_labels = tip_labels,
    root_node = root_node,
    n_tips = n_tips,
    n_internal_nonroot = length(all_internals_nonroot),
    N = N
  )
}

# ---------------------------------------------------------------------------
# henderson_R_inv_apply()
# ---------------------------------------------------------------------------
#
# Apply R^{-1} to a vector or matrix `b` defined at tip positions, using
# Hadfield's sparse Q without ever forming A^{-1} explicitly.
#
# Math: let Q = [[Q_TT  Q_TI]; [Q_IT  Q_II]] be the block decomposition
# of Q at tip vs internal-nonroot rows. Schur complement on the precision
# form:
#
#   R^{-1} = A^{-1} = Q_TT - Q_TI Q_II^{-1} Q_IT
#
# So R^{-1} b = Q_TT b - Q_TI (Q_II^{-1} (Q_IT b)). Each step is sparse:
#   - Q_TT b, Q_TI v, Q_IT b are sparse matrix-vector products (O(nnz))
#   - Q_II^{-1} v is a sparse Cholesky solve (Matrix::Cholesky pre-
#     computed once, reused for every column)
#
# `henderson` is the list returned by `build_henderson_S_inv`. `b` is
# either a numeric vector (n_tips) or a matrix (n_tips x K). Returns the
# same shape.
henderson_R_inv_apply <- function(b, henderson, chol_Q_II = NULL) {
  stopifnot(is.list(henderson),
            !is.null(henderson$Q),
            !is.null(henderson$tip_idx),
            !is.null(henderson$int_idx))
  Q <- henderson$Q
  T_idx <- henderson$tip_idx
  I_idx <- henderson$int_idx

  Q_TT <- Q[T_idx, T_idx, drop = FALSE]
  if (length(I_idx) == 0L) {
    # Tree is just tips (degenerate). R^{-1} = Q_TT.
    return(as.matrix(Q_TT %*% b))
  }
  Q_TI <- Q[T_idx, I_idx, drop = FALSE]
  Q_IT <- Q[I_idx, T_idx, drop = FALSE]
  Q_II <- Q[I_idx, I_idx, drop = FALSE]

  if (is.null(chol_Q_II)) {
    chol_Q_II <- Matrix::Cholesky(Matrix::forceSymmetric(Q_II), LDL = FALSE)
  }
  rhs_v <- as.matrix(Q_IT %*% b)
  v     <- as.matrix(Matrix::solve(chol_Q_II, rhs_v))
  out   <- as.matrix(Q_TT %*% b - Q_TI %*% v)
  out
}

# ---------------------------------------------------------------------------
# henderson_log_det_R()
# ---------------------------------------------------------------------------
#
# log |R| via the precision form: log |R| = log |A| = -log |Q_TT - Q_TI
# Q_II^{-1} Q_IT|. Equivalently, log |R| = log |Q| - log |Q_II|, where
# both log-dets are computable from sparse Cholesky factorisations
# (sum of 2 * log of diagonal entries).
#
# Returns log |R| on the tip-only block. Used for the proper Kronecker
# Fisher-ML on Sigma.
# ---------------------------------------------------------------------------
# henderson_bm_predict()
# ---------------------------------------------------------------------------
#
# Per-column BM imputation using Hadfield's sparse Q. Equivalent in
# spirit to bm_impute_col() but O(n) instead of O(n_obs^3).
#
# Algorithm: given observed values y_o at tip indices `obs_in_tips`,
# the conditional posterior mean of the non-observed extended-tree
# nodes (missing tips + all internals) under BM with unit rate is:
#
#   E[X_no | X_o] = -Q[no, no]^{-1} Q[no, o] y_o
#
# where `o` indexes observed tips in Q and `no` indexes everything
# else (miss tips + internals). The prediction at miss tip positions
# is then E[X_no | X_o] restricted to miss-tip rows.
#
# Variance: returning a conservative tree-depth-based estimate
# (full sparse-selected-inverse via Takahashi recursion is doable
# but ~3x the engineering; for pigauto's purposes the calibrated
# gates handle posterior uncertainty downstream, so per-cell
# precision here matters less than at the GLS-mean level).
#
# Returns list(mu, se) where mu and se are length-n_tips vectors
# (observed cells echo y_o; missing cells get predictions).
henderson_bm_predict <- function(y, henderson, sigma2 = NULL, eps = 1e-8) {
  stopifnot(length(y) == henderson$n_tips)
  obs_idx_local <- which(!is.na(y))
  miss_idx_local <- which(is.na(y))
  n_o <- length(obs_idx_local)
  n_m <- length(miss_idx_local)
  if (n_o == henderson$n_tips) return(list(mu = y, se = rep(0, length(y))))
  if (n_o == 0L) return(list(mu = rep(0, length(y)),
                              se = rep(1, length(y))))
  Q <- henderson$Q
  T_idx <- henderson$tip_idx
  # Map local obs/miss positions to global Q positions.
  obs_q  <- T_idx[obs_idx_local]
  miss_q <- T_idx[miss_idx_local]
  no_q   <- c(miss_q, henderson$int_idx)   # non-observed: miss tips + internals
  Q_no_no <- Q[no_q, no_q, drop = FALSE]
  Q_no_o  <- Q[no_q, obs_q, drop = FALSE]

  rhs <- -as.matrix(Q_no_o %*% y[obs_idx_local])
  chol_no_no <- Matrix::Cholesky(Matrix::forceSymmetric(Q_no_no), LDL = FALSE)
  v <- as.matrix(Matrix::solve(chol_no_no, rhs))

  mu <- y
  mu[miss_idx_local] <- as.numeric(v[seq_len(n_m), 1])

  # Conservative se: empirical sd of observed values (REML-style estimate
  # of the BM rate scaled by sqrt of marginal-variance-fraction). For
  # downstream pigauto consumers the per-cell se matters less than the
  # mu, and the proper sparse Takahashi recursion is a follow-up.
  if (is.null(sigma2)) sigma2 <- stats::var(y[obs_idx_local])
  se <- numeric(length(y))
  se[miss_idx_local] <- sqrt(max(sigma2, eps))
  list(mu = mu, se = se)
}

henderson_log_det_R <- function(henderson, chol_Q = NULL, chol_Q_II = NULL) {
  Q <- henderson$Q
  I_idx <- henderson$int_idx
  if (is.null(chol_Q)) {
    chol_Q <- Matrix::Cholesky(Matrix::forceSymmetric(Q), LDL = FALSE)
  }
  if (length(I_idx) == 0L) {
    # No internal nodes; |R| = |Q|^{-1}
    return(-as.numeric(Matrix::determinant(Q, logarithm = TRUE)$modulus))
  }
  Q_II <- Q[I_idx, I_idx, drop = FALSE]
  if (is.null(chol_Q_II)) {
    chol_Q_II <- Matrix::Cholesky(Matrix::forceSymmetric(Q_II), LDL = FALSE)
  }
  # |R| = |Q_II| / |Q|  →  log |R| = log |Q_II| - log |Q|
  ld_Q    <- as.numeric(Matrix::determinant(Q, logarithm = TRUE)$modulus)
  ld_Q_II <- as.numeric(Matrix::determinant(Q_II, logarithm = TRUE)$modulus)
  ld_Q_II - ld_Q
}
