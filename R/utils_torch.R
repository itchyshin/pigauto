# Internal torch helpers — not exported

# Detect best available device (CUDA > MPS > CPU)
get_device <- function() {
  if (torch::cuda_is_available()) {
    torch::torch_device("cuda")
  } else if (torch::backends_mps_is_available()) {
    torch::torch_device("mps")
  } else {
    torch::torch_device("cpu")
  }
}

# RMSE computed entirely in torch (avoids per-epoch CPU copy)
rmse_torch <- function(pred, truth, mask_bool) {
  d <- pred - truth
  torch::torch_sqrt(torch::torch_mean(d[mask_bool]$pow(2)))
}
