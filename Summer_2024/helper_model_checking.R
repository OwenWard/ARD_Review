## helper function to compute the ppc check of zheng et al 2006



plot_ests <- function(ppc_y, true_y, prop_val = 0) {
  ## takes in two tibbles, ppc_y and true_y, of
  ## the draws of each entry of y and of the true y 
  ## respectively
  
  ppc_prop <- ppc_y |> 
    group_by(sub_pop_id, draw) |> 
    summarise(prop = sum(count == prop_val)/n()) |> 
    group_by(sub_pop_id) |> 
    summarise(avg = mean(prop), 
              lower = quantile(prop, 0.025),
              upper = quantile(prop, 0.975))
  
  true_prop <- true_y |> 
    group_by(sub_pop_id) |> 
    summarise(true_prop = sum(count == prop_val)/n()) 
  
  ## need to set this value in a better way I think
  max_prop <- min(1, max(true_prop$true_prop) + 0.15)
  
  final_plot <- ppc_prop |> 
    left_join(true_prop) |> 
    arrange(true_prop) |> 
    mutate(index = row_number()) |> 
    ggplot(aes(y = true_prop, x = avg)) +
    geom_point() +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.005, alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, alpha = 0.25, col = "red") +
    labs(x = "Simulated", y = "Data",
         subtitle = paste0("Prop = ", prop_val)) +
    # xlim(c(-0.05, max_prop)) +
    # ylim(c(-0.05, max_prop)) +
    # coord_flip() +
    coord_obs_pred() +
    NULL
  
  final_plot
}


## take in the stan fit and the true y matrix and convert to
## tidyverse equiv for ppc plot

construct_ppc <- function(stan_fit, y_sim){
  y_draws <- stan_fit$draws() |> 
    as_draws_df()
  
  ## first get the generated quantities of interest here
  
  ppc_y <- y_draws |> 
    dplyr::select(starts_with("y_sim")) |> 
    mutate(draw = row_number()) |> 
    pivot_longer(cols = starts_with("y_sim"), values_to = "count") |> 
    mutate(node_id = as.numeric(str_extract(name, pattern = "\\d+")),
           sub_pop_id = str_extract(name, pattern = "\\d+]"),
           sub_pop_id = as.numeric(str_replace(sub_pop_id, "\\]", ""))) 
  
  true_y_mat <- y_sim
  
  matrix_df <- as.data.frame(as.table(true_y_mat))
  colnames(matrix_df) <- c("node_id", "sub_pop_id", "count")
  matrix_df$node_id <- as.numeric(matrix_df$node_id)
  matrix_df$sub_pop_id <- as.numeric(matrix_df$sub_pop_id)
  true_y <- as_tibble(matrix_df)
  
  
  list(ppc_draws = ppc_y, y_tibble = true_y)
}
