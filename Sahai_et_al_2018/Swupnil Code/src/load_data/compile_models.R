require("rstan");
options(mc.cores = parallel::detectCores());

# ---------------------------------
# ---------- Stan Models ----------
# ---------------------------------

# model_kernel <- stan_model("src/stan/degree_kernel.stan")
model_kernel <- stan_model(file = paste0(code_path, 
                                         "src/stan/degree_kernel.stan"))
model_kernel_spline <- stan_model(file = paste0(code_path,
                                                "src/stan/degree_kernel_spline.stan"))
