library(torch)
library(ape)
library(Matrix)
library(Rphylopars)
library(ggplot2)
library(patchwork)

# ── 1. Device Setup ──────────────────────────────────────────────────────────
device <- if (cuda_is_available()) torch_device("cuda") else torch_device("cpu")
cat("Using device:", as.character(device), "\n")

# ── 2. Feature Helpers ───────────────────────────────────────────────────────
get_spectral_features <- function(tree, k=8) {
  D <- cophenetic(tree); sigma <- median(D)*0.5
  A <- exp(-(D^2)/(2*sigma^2)); diag(A)<-0; L<-diag(rowSums(A))-A
  eig <- eigen(L, symmetric=TRUE); n<-nrow(A)
  torch_tensor(eig$vectors[,(n-k+1):n], dtype=torch_float(), device=device)
}

get_stable_adj <- function(tree) {
  D <- cophenetic(tree); sigma <- median(D)*0.5
  A <- exp(-(D^2)/(2*sigma^2)); diag(A)<-0; row_sums<-rowSums(A); row_sums[row_sums==0]<-1
  torch_tensor(diag(1/row_sums)%*%A, dtype=torch_float(), device=device)
}

rphylopars_impute_with_se <- function(tree, obs, env) {
  df <- data.frame(species=tree$tip.label, trait=obs, env=env)
  p <- suppressWarnings(phylopars(df, tree, model="BM", pheno_error=FALSE))
  trait_col <- which(colnames(p$anc_recon) %in% c("trait","Trait"))[1]
  list(mu=as.numeric(p$anc_recon[1:length(obs), trait_col]),
       se=sqrt(pmax(as.numeric(p$anc_var[1:length(obs), trait_col]),0)))
}

# ── 3. The Calibrated ResNet Class (Robust Definition) ───────────────────────
PhyloDAE_Calibrated <- nn_module("PhyloDAE_Calibrated",
                                 initialize = function(input_dim, hidden_dim, coord_dim, cov_dim) {
                                   total_input <- input_dim + coord_dim + cov_dim
                                   
                                   # Architecture
                                   self$linear_shortcut <- nn_linear(total_input, input_dim)
                                   self$enc1 <- nn_linear(total_input, 64)
                                   self$enc2 <- nn_linear(64, 64)
                                   self$graph_mix <- nn_linear(64, 64)
                                   self$dec1 <- nn_linear(64, 64)
                                   self$head_mu_deep <- nn_linear(64, input_dim)
                                   self$head_var <- nn_linear(64, input_dim)
                                   
                                   self$act <- nn_relu()
                                   self$drop <- nn_dropout(0.1)
                                 },
                                 
                                 forward = function(x, coords, covs, adj) {
                                   combined <- torch_cat(list(x, coords, covs), dim=2)
                                   lin <- self$linear_shortcut(combined)
                                   
                                   h <- self$act(self$enc1(combined)); h <- self$drop(h)
                                   h <- self$act(self$enc2(h))
                                   h <- h + self$graph_mix(torch_matmul(adj, h))
                                   h <- self$act(self$dec1(h))
                                   
                                   # Prediction = Linear + Deep Correction
                                   list(mu = lin + self$head_mu_deep(h), logvar = self$head_var(h))
                                 }
)

# Helper to Zero-Init the model safely (prevents the crash)
init_weights_safely <- function(model) {
  with_no_grad({
    # Access the underlying tensor data using $data()
    model$head_mu_deep$weight$fill_(0)
    model$head_mu_deep$bias$fill_(0)
    model$head_var$bias$fill_(0.5)
  })
}

# ── 4. Loss Function ─────────────────────────────────────────────────────────
gaussian_nll_loss <- function(pred_mu, pred_logvar, target, mask) {
  precision <- torch_exp(-pred_logvar)
  loss <- 0.5 * (precision * (pred_mu-target)^2 + pred_logvar)
  mask_f <- mask$to(dtype=torch_float())
  torch_sum(loss*mask_f)/torch_sum(mask_f+1e-6)
}

# ── 5. Ensemble Imputation (With Crash Fix) ──────────────────────────────────
impute_ensemble <- function(sim_data, n_models=5, epochs=1500, lr=0.01, k_eigen=8, corruption=0.1) {
  
  trait_vec <- sim_data$obs_trait; tree <- sim_data$tree; env_vec <- sim_data$env
  is_missing <- is.na(trait_vec)
  X_filled <- trait_vec; X_filled[is_missing] <- mean(trait_vec, na.rm=TRUE)
  x_mean <- mean(X_filled); x_sd <- sd(X_filled); if(x_sd==0) x_sd<-1
  X_scaled <- (X_filled - x_mean)/x_sd
  
  X_full <- torch_tensor(matrix(X_scaled,ncol=1), dtype=torch_float(), device=device)
  Env_t  <- torch_tensor(matrix(env_vec,ncol=1), dtype=torch_float(), device=device)
  M_true_obs <- torch_tensor(matrix(!is_missing,ncol=1), dtype=torch_bool(), device=device)
  adj_t <- get_stable_adj(tree); coords_t <- get_spectral_features(tree, k=k_eigen)
  
  ensemble_mu <- matrix(0, nrow=length(trait_vec), ncol=n_models)
  ensemble_se <- matrix(0, nrow=length(trait_vec), ncol=n_models)
  
  cat(sprintf("Training Ensemble (%d models): ", n_models))
  
  for(m in 1:n_models) {
    # 1. Instantiate on CPU first
    model <- PhyloDAE_Calibrated(1, 64, k_eigen, 1)
    
    # 2. Init weights safely
    init_weights_safely(model)
    
    # 3. Move to Device
    model$to(device = device)
    
    opt <- optim_adam(model$parameters, lr=lr)
    
    for(i in 1:epochs) {
      model$train(); opt$zero_grad()
      # Warmup logic
      curr_corr <- if(i < 300) 0.0 else corruption
      noise_mask <- torch_rand_like(X_full) > curr_corr
      out <- model(X_full*noise_mask, coords_t, Env_t, adj_t)
      
      target_mask <- M_true_obs & (!noise_mask)
      if(target_mask$sum()$item() > 0) {
        loss <- gaussian_nll_loss(out$mu, out$logvar, X_full, target_mask)
        loss$backward(); opt$step()
      }
    }
    
    model$eval()
    with_no_grad({
      out <- model(X_full, coords_t, Env_t, adj_t)
      ensemble_mu[,m] <- as.numeric(out$mu$cpu())
      ensemble_se[,m] <- sqrt(exp(as.numeric(out$logvar$cpu())))
    })
    cat(sprintf("%d..", m))
  }
  cat("Done.\n")
  
  avg_mu_scaled <- rowMeans(ensemble_mu)
  final_mu <- (avg_mu_scaled * x_sd) + x_mean
  
  # Total Uncertainty = Average Aleatoric + Variance of Epistemic
  avg_var_scaled <- rowMeans(ensemble_se^2) + apply(ensemble_mu, 1, var)
  final_se <- sqrt(avg_var_scaled) * x_sd
  
  list(imputed=final_mu, se=final_se)
}

# ── 6. Scaling Experiment ────────────────────────────────────────────────────
run_scaling_experiment <- function() {
  # Comparing Small (200) vs Large (1000) Data
  sample_sizes <- c(200, 1000)
  results <- data.frame()
  
  for(N in sample_sizes) {
    cat(sprintf("\n── RUNNING SAMPLE SIZE N=%d ──\n", N))
    
    # Data: Gentle slope (Linear-ish)
    tree <- rtree(N); tree$edge.length <- tree$edge.length + 1e-5
    env <- rTraitCont(tree, model="BM", sigma=0.5)
    env_std <- as.numeric(scale(env))
    phylo_noise <- rTraitCont(tree, model="BM", sigma=0.8)
    true_trait <- 10 / (1 + exp(-1 * env_std)) + phylo_noise
    
    obs_trait <- true_trait; missing_idx <- sample(1:N, size=0.3*N)
    obs_trait[missing_idx] <- NA
    dat <- list(tree=tree, env=env_std, obs_trait=obs_trait)
    
    # 1. Rphylopars
    cat("Rphylopars... ")
    phy <- tryCatch({ rphylopars_impute_with_se(tree, obs_trait, env_std) },
                    error=function(e) list(mu=rep(NA,N)))
    
    # 2. Phylo-DAE Ensemble
    # Note: corruption=0.1 is good for linear-ish data
    res <- impute_ensemble(dat, n_models=5, epochs=1500, lr=0.01, corruption=0.1)
    
    truth <- true_trait[missing_idx]
    rmse_phy <- sqrt(mean((truth - phy$mu[missing_idx])^2))
    rmse_dae <- sqrt(mean((truth - res$imputed[missing_idx])^2))
    
    cat(sprintf("Result N=%d | Phylo-DAE RMSE: %.3f | Rphylopars RMSE: %.3f\n", 
                N, rmse_dae, rmse_phy))
    
    results <- rbind(results, data.frame(N=N, Method="Rphylopars", RMSE=rmse_phy))
    results <- rbind(results, data.frame(N=N, Method="Phylo-DAE (Ens)", RMSE=rmse_dae))
  }
  
  ggplot(results, aes(x=factor(N), y=RMSE, fill=Method)) +
    geom_bar(stat="identity", position="dodge", alpha=0.8) +
    scale_fill_manual(values=c("Phylo-DAE (Ens)"="#1f78b4", "Rphylopars"="#e31a1c")) +
    theme_minimal() +
    labs(title="Impact of Sample Size & Ensembling", 
         subtitle="Neural Networks gain advantage at larger N",
         y="RMSE (Lower is Better)", x="Number of Species (N)")
}

# Run it
run_scaling_experiment()