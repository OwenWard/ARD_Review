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
    #coord_obs_pred() +
    theme_bw() +
    NULL
  
  final_plot
}

plot_ests_all <- function(ppc_y, true_y, prop_vals = 0:10) {
  ## takes in two tibbles, ppc_y and true_y, of
  ## the draws of each entry of y and of the true y 
  ## respectively
  
  count_labs <- paste0("P(y_ik =",prop_vals,")")
  
  ppc_prop <- NULL
  true_prop <- NULL
  
  for (i in 1:length(prop_vals)){
    ppc_prop_i <- ppc_y |> 
      group_by(sub_pop_id, draw) |> 
      summarise( prop = sum(count == prop_vals[i])/n(),
                 .groups = "drop_last") |> # removes warnings
      group_by(sub_pop_id) |> 
      summarise( avg = mean(prop), 
               lower = quantile(prop, 0.025),
               upper = quantile(prop, 0.975) ) %>%
      mutate( count = count_labs[i] )
    
    ppc_prop <- ppc_prop %>%
      bind_rows(ppc_prop_i)
    
    true_prop_i <- true_y |> 
      group_by(sub_pop_id) |> 
      summarise(true_prop = sum(count == prop_vals[i])/n()) %>%
      mutate( count = count_labs[i] )
    
    true_prop <- true_prop %>%
      bind_rows(true_prop_i)
  }
  
  ## need to set this value in a better way I think
  max_prop <- min(1, max(true_prop$true_prop) + 0.15)
  
  final_plot_data <- ppc_prop |> 
    left_join(true_prop) |> 
    arrange(true_prop) |> 
    mutate(index = row_number(),
           count = factor(count,
                          levels = count_labs))
  
  final_plot <- final_plot_data |>
    ggplot(aes(y = true_prop, x = avg)) +
    geom_point() +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.005, alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, alpha = 0.25, col = "red") +
    labs(x = "Simulated", y = "Data",
         #subtitle = paste0("Prop = ", prop_val)
         ) +
    # xlim(c(-0.05, max_prop)) +
    # ylim(c(-0.05, max_prop)) +
    # coord_flip() +
    # coord_obs_pred() +
    facet_wrap(vars(count),scales="free") +
    theme_bw() +
    NULL
  
  return(list( final_plot = final_plot,
               final_plot_data = final_plot_data ))
}


## take in the stan fit and the true y matrix and convert to
## tidyverse equiv for ppc plot

construct_ppc <- function(stan_fit, y_sim){
  if ("CmdStanModel" %in% class(stan_fit)){
    ysim_draws <- stan_fit$draws() |> 
      posterior::as_draws_df() |>
      dplyr::select(starts_with("y_sim"))
  } else if ("stanfit" %in% class(stan_fit)){
    ysim_draws <- rstan::As.mcmc.list( stan_fit,
                                pars = "y_sim" ) %>%
      posterior::as_draws_df() 
  } else if ("draws_df" %in% class(stan_fit)) ysim_draws <- stan_fit
  ## first get the generated quantities of interest here
  
  ppc_y <- ysim_draws |> 
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
