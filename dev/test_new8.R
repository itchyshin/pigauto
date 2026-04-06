library(torch)
library(ape)
library(Matrix)
library(Rphylopars)
library(ggplot2)
library(patchwork)

# ── 1. Graph Helpers ─────────────────────────────────────────────────────────
get_gcn_operators <- function(tree, k=8) {
  D <- cophenetic(tree); sigma <- median(D)*0.5
  A <- exp(-(D^2)/(2*sigma^2)); diag(A)<-1
  deg <- rowSums(A); deg_inv_sqrt <- diag(1/sqrt(deg))
  A_norm <- deg_inv_sqrt %*% A %*% deg_inv_sqrt
  
  L <- diag(deg) - A; eig <- eigen(L, symmetric=TRUE)
  n <- nrow(A); coords <- eig$vectors[,(n-k+1):n]
  
  # FORCE CPU tensors for safety
  list(adj=torch_tensor(A_norm, dtype=torch_float()),
       coords=torch_tensor(coords, dtype=torch_float()))
}

# ── 2. The Residual Model (Stable) ───────────────────────────────────────────
PhyloResidualGCN <- nn_module("PhyloResidualGCN",
                              initialize = function(input_dim, hidden_dim, coord_dim, cov_dim) {
                                total_input <- input_dim + coord_dim + cov_dim
                                
                                # Standard Linear Layers (Stable on all devices)
                                self$proj1 <- nn_linear(total_input, hidden_dim)
                                self$proj2 <- nn_linear(hidden_dim, hidden_dim)
                                
                                self$head_mu <- nn_linear(hidden_dim, input_dim)
                                self$head_var <- nn_linear(hidden_dim, input_dim)
                                
                                self$act <- nn_relu()
                                self$drop <- nn_dropout(0.1)
                              },
                              
                              forward = function(x, coords, covs, adj) {
                                h <- torch_cat(list(x, coords, covs), dim=2)
                                
                                # GCN 1
                                h_trans <- self$proj1(h)
                                h <- torch_matmul(adj, h_trans)
                                h <- self$drop(self$act(h))
                                
                                # GCN 2
                                h_trans2 <- self$proj2(h)
                                h <- torch_matmul(adj, h_trans2)
                                h <- self$act(h)
                                
                                mu <- self$head_mu(h)
                                logvar <- self$head_var(h)
                                
                                # SAFETY CLAMP: Prevent logvar from exploding
                                logvar <- torch_clamp(logvar, min=-5, max=5)
                                
                                return(list(mu=mu, logvar=logvar))
                              }
)

# ── 3. Initialization Helper ─────────────────────────────────────────────────
init_weights_cpu <- function(model) {
  with_no_grad({
    nn_init_xavier_uniform_(model$proj1$weight)
    nn_init_xavier_uniform_(model$proj2$weight)
    model$head_mu$weight$fill_(0)
    model$head_mu$bias$fill_(0)
    model$head_var$bias$fill_(0.5)
  })
}

# ── 4. Loss Function ─────────────────────────────────────────────────────────
gaussian_nll_loss <- function(pred_mu, pred_logvar, target, mask) {
  precision <- torch_exp(-pred_logvar)
  loss <- 0.5 * (precision * (pred_mu - target)^2 + pred_logvar)
  
  mask_f <- mask$to(dtype=torch_float())
  # Add epsilon to denominator
  torch_sum(loss * mask_f) / (torch_sum(mask_f) + 1e-6)
}

# ── 5. The Hybrid Engine (Safe CPU Mode) ─────────────────────────────────────
impute_hybrid <- function(sim_data, epochs=1500, lr=0.005) {
  
  trait_vec <- sim_data$obs_trait; tree <- sim_data$tree; env_vec <- sim_data$env
  n_sp <- length(trait_vec)
  
  # A) Rphylopars Baseline
  df <- data.frame(species=tree$tip.label, trait=trait_vec, env=env_vec)
  phy_model <- suppressWarnings(tryCatch({
    phylopars(df, tree, model="BM", pheno_error=FALSE)
  }, error=function(e) NULL))
  
  if(is.null(phy_model)) return(list(mu=rep(NA, n_sp), se=rep(NA, n_sp)))
  
  trait_col <- which(colnames(phy_model$anc_recon) %in% c("trait","Trait"))[1]
  base_pred <- as.numeric(phy_model$anc_recon[1:n_sp, trait_col])
  base_se   <- sqrt(pmax(as.numeric(phy_model$anc_var[1:n_sp, trait_col]),0))
  
  # B) Residuals
  residuals <- trait_vec - base_pred
  res_mean <- mean(residuals, na.rm=TRUE)
  res_sd   <- sd(residuals, na.rm=TRUE); if(res_sd==0) res_sd<-1
  res_scaled <- (residuals - res_mean) / res_sd
  Res_Full <- res_scaled; Res_Full[is.na(Res_Full)] <- 0 
  
  # FORCE CPU TENSORS (100% Stability)
  X_res <- torch_tensor(matrix(Res_Full, ncol=1), dtype=torch_float())
  Base_Context <- torch_tensor(matrix(scale(base_pred), ncol=1), dtype=torch_float())
  Env_t <- torch_tensor(matrix(env_vec, ncol=1), dtype=torch_float())
  M_obs <- torch_tensor(matrix(!is.na(residuals), ncol=1), dtype=torch_bool())
  
  gcn_ops <- get_gcn_operators(tree) # Returns CPU tensors now
  
  # C) Train GCN
  model <- PhyloResidualGCN(1, 64, 8, 2)
  init_weights_cpu(model)
  # No model$to(device) call needed - defaults to CPU
  
  opt <- optim_adam(model$parameters, lr=lr, weight_decay=1e-3)
  
  for(i in 1:epochs) {
    model$train(); opt$zero_grad()
    
    noise_mask <- torch_rand_like(X_res) > 0.1
    X_corr <- X_res * noise_mask
    
    out <- model(X_corr, gcn_ops$coords, torch_cat(list(Env_t, Base_Context), dim=2), gcn_ops$adj)
    
    target_mask <- M_obs & (!noise_mask)
    if(target_mask$sum()$item() > 0) {
      loss <- gaussian_nll_loss(out$mu, out$logvar, X_res, target_mask)
      loss$backward()
      
      # CRITICAL: Gradient Clipping (Prevents Exploding NaNs)
      nn_utils_clip_grad_norm_(model$parameters, max_norm=1.0)
      
      opt$step()
    }
  }
  
  # D) Combine
  model$eval()
  with_no_grad({
    out <- model(X_res, gcn_ops$coords, torch_cat(list(Env_t, Base_Context), dim=2), gcn_ops$adj)
    pred_res <- as.numeric(out$mu)
    pred_var <- exp(as.numeric(out$logvar))
  })
  
  # Check for NaNs in output
  if(any(is.na(pred_res))) {
    warning("Hybrid GCN produced NaNs. Falling back to base prediction.")
    return(list(mu=base_pred, se=base_se))
  }
  
  final_mu <- base_pred + ((pred_res * res_sd) + res_mean)
  final_se <- sqrt(base_se^2 + (sqrt(pred_var) * res_sd)^2)
  
  return(list(mu=final_mu, se=final_se))
}

# ── 6. Benchmark ─────────────────────────────────────────────────────────────
run_hybrid_benchmark <- function(n_sims=5, n_sp=300) {
  results <- data.frame()
  cat(sprintf("Running Hybrid Benchmark (N=%d)...\n", n_sp))
  plot_data_list <- list()
  
  for(i in 1:n_sims) {
    cat(sprintf("[Sim %d/%d] ", i, n_sims))
    
    # Data: Step Function
    tree <- rtree(n_sp); tree$edge.length <- tree$edge.length + 1e-5
    env <- rTraitCont(tree, model="BM", sigma=0.5)
    env_std <- as.numeric(scale(env))
    phylo_noise <- rTraitCont(tree, model="BM", sigma=0.8)
    
    # TRAP: Step Function (Non-Linear)
    true_trait <- ifelse(env_std > 0, 2, -2) + phylo_noise
    
    obs_trait <- true_trait; missing_idx <- sample(1:n_sp, size=0.3*n_sp)
    obs_trait[missing_idx] <- NA
    dat <- list(tree=tree, env=env_std, obs_trait=obs_trait)
    
    # 1. Rphylopars
    cat("Rphy... ")
    phy_res <- suppressWarnings(tryCatch({
      df <- data.frame(species=tree$tip.label, trait=obs_trait, env=env_std)
      p <- phylopars(df, tree, model="BM", pheno_error=FALSE)
      trait_col <- which(colnames(p$anc_recon) %in% c("trait","Trait"))[1]
      list(mu=as.numeric(p$anc_recon[1:n_sp, trait_col]),
           se=sqrt(pmax(as.numeric(p$anc_var[1:n_sp, trait_col]),0)))
    }, error=function(e) list(mu=rep(NA,n_sp), se=rep(NA,n_sp))))
    
    # 2. Hybrid
    cat("Hybrid... ")
    hybrid_res <- impute_hybrid(dat, epochs=1500, lr=0.005) # Lower LR for safety
    
    # Metrics
    truth <- true_trait[missing_idx]
    
    rmse_phy <- sqrt(mean((truth - phy_res$mu[missing_idx])^2, na.rm=TRUE))
    rmse_hyb <- sqrt(mean((truth - hybrid_res$mu[missing_idx])^2, na.rm=TRUE))
    
    cor_phy <- cor(truth, phy_res$mu[missing_idx], use="complete.obs")
    cor_hyb <- cor(truth, hybrid_res$mu[missing_idx], use="complete.obs")
    
    results <- rbind(results, data.frame(Sim=i, Method="Rphylopars", RMSE=rmse_phy, Cor=cor_phy))
    results <- rbind(results, data.frame(Sim=i, Method="Phylo-Boost", RMSE=rmse_hyb, Cor=cor_hyb))
    
    if(i==n_sims) {
      env_m <- env_std[missing_idx]
      df_p <- rbind(
        data.frame(Method="Rphylopars", Env=env_m, Truth=truth, Pred=phy_res$mu[missing_idx], SE=phy_res$se[missing_idx]),
        data.frame(Method="Phylo-Boost", Env=env_m, Truth=truth, Pred=hybrid_res$mu[missing_idx], SE=hybrid_res$se[missing_idx])
      )
      plot_data_list[[1]] <- df_p
    }
    cat("Done.\n")
  }
  
  list(res=results, plots=plot_data_list[[1]])
}

# ── 7. Execute & Visualize ───────────────────────────────────────────────────
out <- run_hybrid_benchmark(n_sims=5, n_sp=300)
final_res <- out$res
plot_df <- out$plots

cat("\n--- FINAL SCOREBOARD ---\n")
print(aggregate(cbind(RMSE, Cor) ~ Method, data = final_res, mean))

# Visualization
p1 <- ggplot(final_res, aes(x=Method, y=RMSE, fill=Method)) +
  geom_boxplot(alpha=0.6) + geom_jitter(width=0.1) +
  scale_fill_manual(values=c("Phylo-Boost"="#1f78b4", "Rphylopars"="#e31a1c")) +
  theme_minimal() + theme(legend.position="none") + labs(title="RMSE (Lower is Better)")

# Add ribbons
plot_df$Lower <- plot_df$Pred - 1.96 * plot_df$SE
plot_df$Upper <- plot_df$Pred + 1.96 * plot_df$SE
plot_df$Inside <- (plot_df$Truth >= plot_df$Lower) & (plot_df$Truth <= plot_df$Upper)

p3 <- ggplot(plot_df, aes(x=Env)) +
  geom_ribbon(aes(ymin=Lower, ymax=Upper, fill=Method), alpha=0.2) +
  geom_point(aes(y=Truth, color=Inside)) +
  geom_line(aes(y=Pred, color=Method), linewidth=1) +
  facet_wrap(~Method) +
  scale_fill_manual(values=c("Phylo-Boost"="#1f78b4", "Rphylopars"="#e31a1c")) +
  scale_color_manual(values=c("Phylo-Boost"="#1f78b4", "Rphylopars"="#e31a1c", 
                              "TRUE"="grey", "FALSE"="black")) +
  theme_minimal() + labs(title="Calibration Check (Step Function)")

p1 / p3