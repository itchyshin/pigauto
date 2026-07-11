# Harness-only positive controls for the v0.10.0 MI validation.

.numeric_outcome <- function(y) {
  if (is.factor(y)) return(as.numeric(y) - 1)
  as.numeric(y)
}

.draw_oracle_conditional_imputations <- function(dgp, m, seed) {
  m <- as.integer(m)
  if (length(m) != 1L || is.na(m) || m < 2L) {
    stop("m must be an integer >= 2.", call. = FALSE)
  }
  missing <- !is.finite(dgp$observed$x)
  if (!any(missing)) return(rep(list(dgp$observed), m))

  mu_prior <- dgp$x_prior_mean[missing]
  sd_prior <- rep(dgp$x_prior_sd, sum(missing))
  y <- .numeric_outcome(dgp$observed$y[missing])
  z <- dgp$observed$z[missing]
  species <- if (dgp$dgp == "lmer") dgp$observed$species[missing] else NULL

  draw_one <- function(i) {
    set.seed(as.integer(seed + i))
    if (dgp$dgp %in% c("lm", "lmer")) {
      beta_x <- 0.60
      beta_z <- -0.40
      random_intercept <- if (dgp$dgp == "lmer") {
        unname(dgp$random_intercept[species])
      } else {
        rep(0, length(y))
      }
      outcome_offset <- 1 + beta_z * z + random_intercept
      posterior_var <- 1 / (1 / sd_prior^2 + beta_x^2)
      posterior_mean <- posterior_var * (
        mu_prior / sd_prior^2 + beta_x * (y - outcome_offset)
      )
      x_draw <- stats::rnorm(
        length(y), mean = posterior_mean, sd = sqrt(posterior_var)
      )
    } else {
      x_draw <- numeric(length(y))
      for (j in seq_along(y)) {
        accepted <- FALSE
        while (!accepted) {
          proposal <- stats::rnorm(1L, mu_prior[j], sd_prior[j])
          probability <- stats::plogis(-0.20 + 0.80 * proposal - 0.50 * z[j])
          likelihood <- if (y[j] == 1) probability else 1 - probability
          accepted <- stats::runif(1L) <= likelihood
        }
        x_draw[j] <- proposal
      }
    }
    completed <- dgp$observed
    completed$x[missing] <- x_draw
    completed
  }

  lapply(seq_len(m), draw_one)
}

.draw_smcfcs_imputations <- function(dgp, m, seed, numit = 20L) {
  if (!requireNamespace("smcfcs", quietly = TRUE)) {
    stop("The standard single-level comparator requires package 'smcfcs'.",
         call. = FALSE)
  }
  data <- dgp$observed[, c("x", "y", "z"), drop = FALSE]
  data$z_sq <- data$z^2
  if (dgp$dgp == "glm") data$y <- .numeric_outcome(data$y)
  method <- c(x = "norm", y = "", z = "", z_sq = "")
  predictor_matrix <- matrix(
    0, nrow = ncol(data), ncol = ncol(data),
    dimnames = list(names(data), names(data))
  )
  predictor_matrix["x", c("z", "z_sq")] <- 1

  set.seed(as.integer(seed))
  fit <- smcfcs::smcfcs(
    originaldata = data,
    smtype = if (dgp$dgp == "glm") "logistic" else "lm",
    smformula = "y ~ x + z",
    method = method,
    predictorMatrix = predictor_matrix,
    m = as.integer(m),
    numit = as.integer(numit),
    rjlimit = 10000L,
    noisy = FALSE
  )
  lapply(fit$impDatasets, function(completed) {
    completed$z_sq <- NULL
    completed$y <- dgp$observed$y
    completed
  })
}

.draw_jomo_lmer_imputations <- function(dgp, m, seed, nburn = 1000L,
                                         nbetween = 100L) {
  if (!requireNamespace("jomo", quietly = TRUE)) {
    stop("The standard multilevel comparator requires package 'jomo'.",
         call. = FALSE)
  }
  data <- dgp$observed[, c("y", "x", "z", "species"), drop = FALSE]
  data$z_sq <- data$z^2
  data$species <- factor(data$species)
  level <- rep(1L, ncol(data))

  set.seed(as.integer(seed))
  long <- jomo::jomo.smc(
    formula = y ~ x + z + (1 | species),
    data = data,
    level = level,
    nburn = as.integer(nburn),
    nbetween = as.integer(nbetween),
    nimp = as.integer(m),
    output = 0,
    model = "lmer"
  )
  lapply(seq_len(m), function(i) {
    completed <- long[long$Imputation == i, , drop = FALSE]
    completed <- completed[order(completed$id), c("y", "x", "z", "z_sq"),
                           drop = FALSE]
    completed$z_sq <- NULL
    completed$species <- as.character(dgp$observed$species)
    rownames(completed) <- NULL
    completed
  })
}

.draw_standard_smc_imputations <- function(dgp, m, seed) {
  if (dgp$dgp == "lmer") {
    .draw_jomo_lmer_imputations(dgp, m, seed)
  } else {
    .draw_smcfcs_imputations(dgp, m, seed)
  }
}
