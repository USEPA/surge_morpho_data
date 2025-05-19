# Load Packages
library(sf)
library(terra)
library(dplyr)
library(readr)
library(here)
library(tidyr)
library(purrr)
library(ggplot2)
library(patchwork)

# Calculate "actual" metrics from bathy
clipped_bathy_files <- list.files(here("data/bathymetry_clipped"),
                                  full.names = TRUE)
surge_bathy_res <- st_read(here("data/bathymetry_reservoirs.gpkg"))

calc_actual_metrics <- function(file, reserv){
  lake_id <- strsplit(basename(file), "_")[[1]][1]
  reserv <- reserv[reserv$lake_id == as.numeric(lake_id),]
  bathy_c <- rast(file)
  bathy_c_values <- values(bathy_c)[!is.na(values(bathy_c))]
  max_depth <- max(bathy_c_values, na.rm = TRUE)
  mean_depth <- mean(bathy_c_values, na.rm = TRUE)
  volume <- sum(prod(res(bathy_c)) * bathy_c_values, na.rm = TRUE)
  area <- st_area(reserv)
  mean_depth_va <- volume/area
  perc_less_2_5 <- sum(bathy_c_values <= 2.5)/length(bathy_c_values) * 100
  q_bathy_shape <- max_depth/mean_depth - 1
  perc_litt <- (1-(1-2.5/max_depth)^q_bathy_shape)*100
  data.frame(lake_id = lake_id, surface_area = as.numeric(area),
             max_depth = max_depth, mean_depth = mean_depth,
             volume = volume, perc_less_2_5 = perc_less_2_5,
             perc_litt = perc_litt)
}

actual_metrics_1 <- lapply(clipped_bathy_files, calc_actual_metrics,
                         reserv = surge_bathy_res)
actual_metrics <- bind_rows(actual_metrics_1) |>
  mutate(lake_id = as.numeric(lake_id), source = "bathymetry") |>
  pivot_longer(surface_area:perc_litt, names_to = "variables", values_to = "values") |>
  unique()


# Compare estimated data to bathymetry actual
surge_morpho <- readr::read_csv(here("data/surge_morpho.csv")) |>
  filter(lake_id %in% unique(actual_metrics$lake_id))
compare_data <- bind_rows(actual_metrics, surge_morpho) |>
  select(-lake_name) |>
  arrange(lake_id) |>
  filter(variables %in% c("max_depth", "mean_depth", "volume"))

compare_data_w <- compare_data |>
  mutate(source = case_when(source == "bathymetry" ~
                              "b",
                            source == "surge_sites" ~
                              "s",
                            TRUE ~ source)) |>
  pivot_wider(id_cols = lake_id, names_from = c("source", "variables"),
              names_sep = "_", values_from = "values") |>
  na.omit()

rmse_vol <- as.numeric(units::set_units(units::set_units(
  sqrt(mean((compare_data_w$b_volume - compare_data_w$s_volume)^2)), "m3"),
  "km3"))
rmse_mean_depth <- sqrt(mean((compare_data_w$b_mean_depth - compare_data_w$s_mean_depth)^2))
rmse_max_depth <- sqrt(mean((compare_data_w$b_max_depth - compare_data_w$s_max_depth)^2))

rmse_vol <-round(rmse_vol, 2)
rmse_mean_depth <- round(rmse_mean_depth, 2)
rmse_max_depth <- round(rmse_max_depth, 2)


# Figure - scatterplot of actual vs estimated
#' @param x vector of x values
#' @param y vector of y values
gg_scatter <- function(x, y, x_lab, y_lab, title, rmse){
  df <- data.frame(x,y)
  lim <- range(c(x, y))
  r2 <- round(summary(lm(y ~ x))$r.squared, 3)
  beta <- round(summary(lm(y ~ x))$coefficients[2,1], 3)
  int <- round(summary(lm(y ~ x))$coefficients[1,1], 3)
  lbl <- paste("R^2 == ", r2)
  rmse <- paste("RMSE:", rmse)
  ggplot(df, aes(x = x, y = y)) +
    geom_point() +
    #geom_smooth(method = "lm") +
    geom_abline(slope = 1, intercept = 0) +
    hrbrthemes::theme_ipsum_rc() +
    scale_x_continuous(limits = lim, n.breaks = 4) +
    scale_y_continuous(limits = lim, n.breaks = 4) +
    coord_fixed() +
    annotate("text", y = lim[2]*0.95, x = lim[2]*0.25, label = rmse, parse = TRUE) +
    xlab(x_lab) +
    ylab(y_lab) +
    labs(title = title) +
    theme(plot.title = element_text(size = 14))
}

max_depth_valid <- gg_scatter(compare_data_w$b_max_depth,
                              compare_data_w$s_max_depth,
                              x_lab = "Bathymetry Maximum Depth (m)",
                              y_lab = "Sampled Maximum Depth (m)",
                              title = "A",
                              rmse_max_depth)

mean_depth_valid <- gg_scatter(compare_data_w$b_mean_depth,
                              compare_data_w$s_mean_depth,
                              x_lab = "Bathymetry Mean Depth (m)",
                              y_lab = "Sampled Mean Depth (m)",
                              title = "B", rmse_mean_depth)

volume_valid <- gg_scatter(compare_data_w$b_volume/(1000*1000*1000),
                           compare_data_w$s_volume/(1000*1000*1000),
                              x_lab = "Bathymetry Volume (cubic km)",
                              y_lab = "Sampled Volume Depth (cubic km)",
                              title = "C", rmse_vol)

valid <- max_depth_valid + mean_depth_valid + volume_valid

ggsave("validate.jpg", valid, dpi = 600, width = 10, height = 7, units = "in")
