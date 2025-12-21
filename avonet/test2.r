# ── 5. Clade-Priority Model Architecture (Dimension Fixed) ──────────────────
PhyloDAE_V3 <- nn_module(
  "PhyloDAE_V3",
  initialize = function(n_traits, n_phylo, hidden_dim = 512) {
    n_traits <- as.integer(n_traits)
    n_phylo  <- as.integer(n_phylo)
    hidden_dim <- as.integer(hidden_dim)
    
    # DYNAMIC DIMENSION: traits (4) + mask (4) + phylo (k)
    input_dim <- as.integer(n_traits + n_traits + n_phylo)
    
    self$fc1 <- nn_linear(input_dim, hidden_dim)
    self$bn1 <- nn_batch_norm1d(hidden_dim)
    
    self$gcn1 <- nn_linear(hidden_dim, hidden_dim)
    self$gate <- nn_linear(hidden_dim, hidden_dim)
    
    self$out <- nn_linear(hidden_dim, n_traits)
    self$act <- nn_gelu()
    self$drop <- nn_dropout(0.05)
    
    self$mask_token <- nn_parameter(torch_randn(1, n_traits, dtype = torch_float()))
  },
  
  forward = function(x, m, p, adj) {
    x_input <- x * m + (1 - m) * self$mask_token
    h_init  <- torch_cat(list(x_input, m, p), dim = as.integer(2))
    
    # 1. Phylogenetic Aggregation (First Pass)
    h_agg <- torch_matmul(adj, h_init)
    
    h <- self$fc1(h_agg) %>% self$bn1() %>% self$act()
    
    # 2. Residual Gating (Phylogenetic Smoothing)
    h_res <- h
    h_graph <- torch_matmul(adj, h)
    g <- torch_sigmoid(self$gate(h))
    h <- h_res + g * self$act(self$gcn1(h_graph))
    
    self$out(h)
  }
)

# ── 6. Training with Dynamic Dimension Check ──────────────────────────────
X_fill <- as.matrix(X_with_NA)
X_fill[is.na(X_fill)] <- 0

t_X            <- torch_tensor(X_fill, dtype = torch_float(), device = device)
t_mask         <- torch_tensor(as.matrix(mask_true), dtype = torch_float(), device = device)
t_truth_tensor <- torch_tensor(as.matrix(X_truth), dtype = torch_float(), device = device)
t_phylo        <- torch_tensor(as.matrix(phylo_mat), dtype = torch_float(), device = device)
t_adj          <- torch_tensor(as.matrix(adj_mat), dtype = torch_float(), device = device)

# GET CURRENT DIMENSIONS FROM DATA
curr_n_traits <- as.integer(ncol(X_truth))
curr_k_phylo  <- as.integer(ncol(phylo_mat)) # This will be 32 now

model <- PhyloDAE_V3(
  n_traits   = curr_n_traits, 
  n_phylo    = curr_k_phylo, 
  hidden_dim = as.integer(512)
)$to(device = device)

optimizer <- optim_adamw(model$parameters, lr = 5e-4, weight_decay = 1e-2)
scheduler <- lr_one_cycle(optimizer, max_lr = 0.002, epochs = 10000, steps_per_epoch = 1)

trait_weights <- torch_tensor(c(5.0, 1.0, 1.0, 1.0), device = device)$view(as.integer(c(1, 4)))

cat(sprintf("\n--- Training V3 (Clade Priority) | Input Dim: %d ---\n", 
            curr_n_traits * 2 + curr_k_phylo))

for (epoch in 1:10000) {
  model$train()
  optimizer$zero_grad()
  
  u <- torch_rand_like(t_X)
  probs <- torch_tensor(c(0.3, 0.1, 0.1, 0.1), device = device)
  training_mask <- (u < probs) & (t_mask == 1)
  
  curr_mask_tensor <- t_mask * (1 - training_mask$to(torch_float()))
  
  preds <- model(t_X, curr_mask_tensor, t_phylo, t_adj)
  
  diff_sq <- (preds - t_truth_tensor)^2
  weighted_diff <- diff_sq * trait_weights
  loss <- weighted_diff[training_mask]$mean()
  
  loss$backward()
  optimizer$step()
  scheduler$step()
  
  if (epoch %% 1000 == 0) {
    cat(sprintf("Epoch %d | Loss: %.6f\n", epoch, as.numeric(loss$cpu())))
  }
}