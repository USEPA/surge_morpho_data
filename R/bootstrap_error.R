
boot_error <- function(df,
                       param = c("volume", "mean_depth","max_depth"),
                       metric = c("perc_error", "rmse")){

  param <- match.arg(param)
  metric <- match.arg(metric)
  df_b <- df[sample(1:nrow(df), 100, replace = TRUE),]
  if(param == "volume"){
    perc_err <- mean(abs(df_b$b_volume - df_b$s_volume)/df_b$b_volume*100)
    rmse <- sqrt(mean((df_b$b_volume - df_b$s_volume)^2))
  } else if(param == "mean_depth"){
    perc_err <- mean(abs(df_b$b_mean_depth - df_b$s_mean_depth)/df_b$b_mean_depth*100)
    rmse <- sqrt(mean((df_b$b_mean_depth - df_b$s_mean_depth)^2))
  } else if(param == "max_depth"){
    perc_err <- mean(abs(df_b$b_max_depth - df_b$s_max_depth)/df_b$b_max_depth*100)
    rmse <- sqrt(mean((df_b$b_max_depth - df_b$s_max_depth)^2))
  }
  if(metric == "perc_error"){
    err <- perc_error
  } else if(metric == "rmse"){
    err <- rmse
  }
  err
}

formula <- vector("numeric", 100)
for(i in 1:100){
  formula[i] <- boot_error(compare_data_w, "volume", "rmse")
}
err_df_vol <- tibble(variable = "formula", value = formula) |>
  mutate(value = as.numeric(units::set_units(units::set_units(value, "m3"), "acre feet")))
meand <- vector("numeric", )
for(i in 1:100){
  meand[i] <- boot_error(compare_data_w, "mean_depth", "rmse")
}
err_df_meand <- tibble(variable = "mean_depth", value = meand)
maxd <- vector("numeric", 3)
for(i in 1:100){
  maxd[i] <- boot_error(compare_data_w, "max_depth", "rmse")
}
err_df_maxd <- tibble(variable = "max_depth", value = maxd)
err_df <- bind_rows(err_df_vol, err_df_meand, err_df_maxd)
boxplot(value ~ variable, data = err_df)
