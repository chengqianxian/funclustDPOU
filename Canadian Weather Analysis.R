############################################################
## Real data analysis: Canadian Weather data
## fda::CanadianWeather
## Daily temperature curves
## 35 stations, 365 daily time points
## B-spline basis M = 6, truncation H = 8
############################################################

library(fda)
library(splines)
library(ggplot2)

set.seed(1)

############################################################
## Load Canadian Weather data
############################################################

data(CanadianWeather)

str(CanadianWeather)

############################################################
## Extract daily temperature data
############################################################

dailyAv <- CanadianWeather$dailyAv

# Usually dailyAv has dimensions: day x station x variable
dim(dailyAv)
dimnames(dailyAv)

# Robustly select temperature variable
var_names <- dimnames(dailyAv)[[3]]

if (!is.null(var_names)) {
  temp_idx <- grep("Temperature", var_names, ignore.case = TRUE)
  if (length(temp_idx) == 0) temp_idx <- 1
} else {
  temp_idx <- 1
}

Y_raw <- dailyAv[, , temp_idx]

# Rows = stations, columns = days
Y_mat <- t(Y_raw)

N <- nrow(Y_mat)
n_grid <- ncol(Y_mat)

station_names <- rownames(Y_mat)
if (is.null(station_names)) {
  station_names <- paste0("Station_", seq_len(N))
}

############################################################
## Time grid
############################################################

day_grid <- seq_len(n_grid)

# Rescale days to [0, 1] for OU covariance and B-spline basis
t_grid <- (day_grid - min(day_grid)) / (max(day_grid) - min(day_grid))

############################################################
## B-spline basis, M = 6
############################################################

M <- 6
H <- 8

B <- splines::bs(
  x = t_grid,
  df = M,
  degree = 3,
  intercept = TRUE
)

B <- as.matrix(B)

############################################################
## Construct model lists
############################################################

Y_list <- vector("list", N)
B_list <- vector("list", N)
t_list <- vector("list", N)

for (i in seq_len(N)) {
  Y_list[[i]] <- as.numeric(Y_mat[i, ])
  B_list[[i]] <- B
  t_list[[i]] <- t_grid
}

############################################################
## Quick visualization
############################################################

matplot(
  day_grid,
  t(Y_mat),
  type = "l",
  lty = 1,
  col = adjustcolor("grey50", alpha.f = 0.6),
  xlab = "Day",
  ylab = "Temperature",
  main = "Canadian Weather Daily Temperature Curves"
)


############################################################
## Hyperparameters
############################################################

m0 <- rep(0, M)

# Temperature scale is moderate; use diffuse prior
S0 <- diag(100, M)

a0 <- 2
b0 <- 1
alpha <- 1

############################################################
## Initialization by k-means
############################################################

set.seed(1)

km <- kmeans(Y_mat, centers = H, nstart = 20)

r_init <- matrix(1e-3, nrow = N, ncol = H)

for (i in seq_len(N)) {
  r_init[i, km$cluster[i]] <- 1
}

r_init <- r_init / rowSums(r_init)

############################################################
## Initialize stick-breaking parameters
############################################################

gamma_1_init <- vapply(seq_len(H), function(h) {
  update_gamma_h1(r_init[, h])
}, numeric(1))

gamma_2_init <- vapply(seq_len(H), function(h) {
  update_gamma_h2(alpha, r_init, h)
}, numeric(1))

############################################################
## Initialize tau parameters
############################################################

a_tilde_init <- rep(a0 + 0.5 * n_grid * N / H, H)
b_tilde_init <- rep(b0 + 100, H)

############################################################
## Fit VBEM model
############################################################

VBEMfit_weather <- variational_EM_DP_OU(
  Y_list = Y_list,
  B_list = B_list,
  t_list = t_list,
  H = H,
  m0 = m0,
  S0 = S0,
  a0 = a0,
  b0 = b0,
  alpha = alpha,
  init = list(
    r = r_init,
    a_tilde = a_tilde_init,
    b_tilde = b_tilde_init,
    gamma_1 = gamma_1_init,
    gamma_2 = gamma_2_init,
    delta = 8
  ),
  max_iter = 200,
  tol = 1e-6,
  delta_bounds = c(1e-4, 30),
  update_delta = TRUE,
  jitter = 1e-8,
  verbose = TRUE
)

############################################################
## Diagnostics
############################################################

delta_hat <- VBEMfit_weather$delta

exp(-delta_hat / 364)
exp(-delta_hat * 7 / 364)
exp(-delta_hat * 10 / 364)

plot(
  VBEMfit_weather$elbo_history,
  type = "l",
  xlab = "Iteration",
  ylab = "ELBO",
  main = "ELBO convergence: Canadian Weather data"
)

plot(
  VBEMfit_weather$delta_history,
  type = "l",
  xlab = "Iteration",
  ylab = expression(delta),
  main = "Delta estimates: Canadian Weather data"
)

cluster_size <- colSums(VBEMfit_weather$r)
cluster_size

active_clusters <- which(cluster_size > 5)
active_clusters

weather_cluster <- VBEMfit_weather$cluster

table(weather_cluster)


############################################################
## Prepare ggplot data
############################################################

active_clusters <- which(cluster_size > 5)

raw_df <- data.frame(
  day = rep(day_grid, times = N),
  temperature = as.vector(t(Y_mat)),
  station = rep(station_names, each = n_grid),
  cluster = factor(rep(weather_cluster, each = n_grid))
)

estimated_list <- vector("list", length(active_clusters))

for (j in seq_along(active_clusters)) {
  h <- active_clusters[j]
  
  mu_h <- VBEMfit_weather$mu_list[[h]]
  Sigma_h <- VBEMfit_weather$Sigma_list[[h]]
  
  mean_h <- as.vector(B %*% mu_h)
  se_h <- sqrt(rowSums((B %*% Sigma_h) * B))
  
  lower_h <- mean_h - 1.96 * se_h
  upper_h <- mean_h + 1.96 * se_h
  
  estimated_list[[j]] <- data.frame(
    day = day_grid,
    mean = mean_h,
    lower = lower_h,
    upper = upper_h,
    component = factor(h),
    cluster = factor(h),
    n_eff = round(cluster_size[h], 1)
  )
}

estimated_df <- do.call(rbind, estimated_list)

############################################################
## ggplot: raw curves + estimated cluster means + 95% CI
############################################################

p_weather <- ggplot() +
  geom_line(
    data = raw_df,
    aes(
      x = day,
      y = temperature,
      group = station
    ),
    color = "grey75",
    alpha = 0.45,
    linewidth = 0.35
  ) +
  geom_ribbon(
    data = estimated_df,
    aes(
      x = day,
      ymin = lower,
      ymax = upper,
      fill = cluster,
      group = component
    ),
    alpha = 0.18,
    color = NA
  ) +
  geom_line(
    data = estimated_df,
    aes(
      x = day,
      y = mean,
      color = cluster,
      group = component
    ),
    linewidth = 1.1
  ) +
  labs(
    x = "Day",
    y = "Daily average temperature",
    color = "Estimated cluster",
    fill = "Estimated cluster"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "right"
  )

p_weather

p_weather_clean <- ggplot() +
  geom_ribbon(
    data = estimated_df,
    aes(
      x = day,
      ymin = lower,
      ymax = upper,
      fill = cluster,
      group = component
    ),
    alpha = 0.20,
    color = NA
  ) +
  geom_line(
    data = estimated_df,
    aes(
      x = day,
      y = mean,
      color = cluster,
      group = component
    ),
    linewidth = 1.2
  ) +
  labs(
    x = "Day",
    y = "Daily average temperature",
    color = "Estimated cluster",
    fill = "95% CI"
  ) +
  theme_bw(base_size = 13)

p_weather_clean


library(fda)
library(ggplot2)
library(maps)

data(CanadianWeather)

############################################################
## Extract station metadata
############################################################

station_names <- CanadianWeather$place
province <- CanadianWeather$province
coords <- CanadianWeather$coordinates

head(coords)


############################################################
## Build station data frame
############################################################

weather_cluster <- VBEMfit_weather$cluster

cluster_size <- colSums(VBEMfit_weather$r)

station_df <- data.frame(
  station = station_names,
  province = province,
  latitude = coords[, "N.latitude"],
  longitude = -coords[, "W.longitude"],
  cluster = factor(weather_cluster)
)

head(station_df)

############################################################
## Canada map
############################################################

canada_map <- map_data("world", region = "Canada")

p_map <- ggplot() +
  geom_polygon(
    data = canada_map,
    aes(x = long, y = lat, group = group),
    fill = "grey90",
    color = "grey60",
    linewidth = 0.3
  ) +
  geom_point(
    data = station_df,
    aes(
      x = longitude,
      y = latitude,
      color = cluster
    ),
    size = 3,
    alpha = 0.9
  ) +
  coord_quickmap(
    xlim = c(-145, -50),
    ylim = c(40, 75)
  ) +
  labs(
    x = "Longitude",
    y = "Latitude",
    color = "Estimated cluster"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "right"
  )

p_map

############################################################
## Only show active clusters on map
############################################################

active_threshold <- 5   # 或者你之前用的阈值，比如 5
cluster_size <- colSums(VBEMfit_weather$r)
active_clusters <- which(cluster_size > active_threshold)

active_clusters
cluster_size[active_clusters]

############################################################
## Build station data frame
############################################################

weather_cluster <- VBEMfit_weather$cluster

lat_col <- grep("lat", colnames(CanadianWeather$coordinates), ignore.case = TRUE, value = TRUE)[1]
lon_col <- grep("long", colnames(CanadianWeather$coordinates), ignore.case = TRUE, value = TRUE)[1]

station_df <- data.frame(
  station = CanadianWeather$place,
  province = CanadianWeather$province,
  latitude = CanadianWeather$coordinates[, lat_col],
  longitude = -abs(CanadianWeather$coordinates[, lon_col]),
  cluster = weather_cluster
)

############################################################
## Keep only active clusters
############################################################

station_df_active <- subset(
  station_df,
  cluster %in% active_clusters
)

station_df_active$cluster <- factor(station_df_active$cluster)

############################################################
## Plot active clusters only
############################################################

library(ggplot2)
library(maps)

canada_map <- map_data("world", region = "Canada")

p_map_active <- ggplot() +
  geom_polygon(
    data = canada_map,
    aes(x = long, y = lat, group = group),
    fill = "grey90",
    color = "grey60",
    linewidth = 0.3
  ) +
  geom_point(
    data = station_df_active,
    aes(
      x = longitude,
      y = latitude,
      color = cluster
    ),
    size = 3,
    alpha = 0.9
  ) +
  coord_quickmap(
    xlim = c(-145, -50),
    ylim = c(40, 75)
  ) +
  labs(
    x = "Longitude",
    y = "Latitude",
    color = "Active cluster"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "right"
  )

p_map_active

station_df$status <- ifelse(
  station_df$cluster %in% active_clusters,
  "active",
  "inactive"
)

station_df$cluster_active <- ifelse(
  station_df$cluster %in% active_clusters,
  as.character(station_df$cluster),
  NA
)

station_df$plot_group <- ifelse(
  station_df$cluster %in% active_clusters,
  paste0("Cluster ", station_df$cluster),
  "Inactive"
)

station_df$plot_group <- factor(
  station_df$plot_group,
  levels = c("Inactive", paste0("Cluster ", active_clusters))
)

p_map_with_inactive_legend <- ggplot() +
  geom_polygon(
    data = canada_map,
    aes(x = long, y = lat, group = group),
    fill = "grey90",
    color = "grey60",
    linewidth = 0.3
  ) +
  geom_point(
    data = station_df,
    aes(
      x = longitude,
      y = latitude,
      color = plot_group,
      alpha = plot_group
    ),
    size = 3
  ) +
  scale_color_manual(
    values = c(
      "Inactive" = "grey60",
      setNames(
        scales::hue_pal()(length(active_clusters)),
        paste0("Cluster ", active_clusters)
      )
    )
  ) +
  scale_alpha_manual(
    values = c(
      "Inactive" = 0.45,
      setNames(rep(0.95, length(active_clusters)), paste0("Cluster ", active_clusters))
    ),
    guide = "none"
  ) +
  coord_quickmap(
    xlim = c(-145, -50),
    ylim = c(40, 75)
  ) +
  labs(
    x = "Longitude",
    y = "Latitude",
    color = "Station group"
  ) +
  theme_bw(base_size = 13)

p_map_with_inactive_legend
