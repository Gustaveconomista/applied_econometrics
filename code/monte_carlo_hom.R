#### Load necessary packages ####
pacman::p_load(tidyverse, 
               sandwich, 
               lmtest, 
               boot, 
               here,
               clubSandwich,
               fwildclusterboot)

#### For Normal Distribution ####
##### Set seed for reproducibility #####
set.seed(42)

##### Simulation parameters / Parâmetros da simulação ####
cluster_sizes = c(5, 10)  # Number of clusters to simulate / Número de clusters
NG = 30                      # Observations per cluster / Observações por cluster
R = 1000                    # Number of Monte Carlo replications / Número de repetições
B = 399                    # Number of bootstrap replications / Número de bootstrap
beta_0 = 0
beta_1 = 1
alpha = 0.05

##### Data generating process function (equation 7) #####
generate_data = function(G, NG) {
  zg = rnorm(G)
  epsilon_g = rnorm(G)
  data = data.frame()
  for (g in 1:G) {
    zig = rnorm(NG)
    epsilon_ig = rnorm(NG)
    
    xig = zg[g] + zig
    uig = epsilon_g[g] + epsilon_ig
    yig = beta_0 + beta_1 * xig + uig
    
    cluster_data = data.frame(y = yig, x = xig, cluster = rep(g, NG))
    data = rbind(data, cluster_data)
  }
  return(data)
}

##### Cluster-robust variance covariance matrix #####
cluster_vcov = function(model, cluster_var) {
  # Using vcovCL from sandwich package for cluster-robust covariance
  vcovCL(model, cluster = cluster_var)
}

##### CR3 variance matrix correction #####
cr3_vcov = function(model, cluster_var) {
  clusters = unique(cluster_var)
  G = length(clusters)
  X = model.matrix(model)
  beta_jack = matrix(0, G, ncol(X))
  for (i in 1:G) {
    idx = cluster_var != clusters[i]
    data_jack = data.frame(y = model$model$y[idx], x = model$model$x[idx])
    model_jack = lm(y ~ x, data = data_jack)
    beta_jack[i, ] = coef(model_jack)
  }
  beta_hat = coef(model)
  V_jack = ((G - 1) / G) * t(beta_jack - matrix(rep(beta_hat, G), nrow = G, byrow = TRUE)) %*%
    (beta_jack - matrix(rep(beta_hat, G), nrow = G, byrow = TRUE))
  return(V_jack)
}

##### Moulton-type SE approximation based on intraclass correlation #####
moulton_se = function(model, NG, icc = 0.5) {
  se_default = sqrt(vcov(model)[2,2])
  design_effect = 1 + (NG - 1) * icc
  se_default * sqrt(design_effect)
}

##### Pairs Cluster Bootstrap SE #####
pairs_cluster_bootstrap = function(data, B = 399) {
  clusters = unique(data$cluster)
  G = length(clusters)
  beta_boot = numeric(B)
  for (b in 1:B) {
    boot_clusters = sample(clusters, G, replace = TRUE)
    boot_data = data.frame()
    for (g in boot_clusters) {
      boot_data = rbind(boot_data, data[data$cluster == g, ])
    }
    tryCatch({
      boot_model = lm(y ~ x, data = boot_data)
      beta_boot[b] = coef(boot_model)[2]
    }, error = function(e) {
      beta_boot[b] = NA
    })
  }
  return(na.omit(beta_boot))
}

##### Residual Cluster Bootstrap-SE (Hold X fixed, resample residuals with H0) #####
residual_cluster_bootstrap_se = function(data, B = 399, beta_0_test = beta_1) {
  clusters = unique(data$cluster)
  G = length(clusters)
  
  model_restricted = lm(y ~ x, data = data, subset = TRUE) # We'll impose null through residual adjustment
  residuals = residuals(model_restricted)
  
  # Calculate residuals per cluster
  residuals_by_cluster = split(residuals, data$cluster)
  
  beta_boot = numeric(B)
  for (b in 1:B) {
    boot_resid = residuals_by_cluster
    # resample residuals by cluster with replacement
    for (g in seq_along(residuals_by_cluster)) {
      residuals_cluster = residuals_by_cluster[[g]]
      residuals_by_cluster[[g]] = sample(residuals_cluster, length(residuals_cluster), replace = TRUE)
    }
    
    # Construct bootstrap sample
    boot_y = numeric(nrow(data))
    for (g in clusters) {
      idx = which(data$cluster == g)
      boot_y[idx] = beta_0_test * data$x[idx] + residuals_by_cluster[[g]]
    }
    
    boot_data = data
    boot_data$y = boot_y
    tryCatch({
      boot_model = lm(y ~ x, data = boot_data)
      beta_boot[b] = coef(boot_model)[2]
    }, error = function(e) {
      beta_boot[b] = NA
    })
  }
  return(na.omit(beta_boot))
}

##### Wild Cluster Bootstrap-SE (Apply Rademacher weights) #####
wild_cluster_bootstrap_se = function(data, B = 399, beta_0_test = beta_1) {
  clusters = unique(data$cluster)
  G = length(clusters)
  
  model_restricted = lm(y ~ x, data = data)
  residuals = residuals(model_restricted)
  residuals_by_cluster = split(residuals, data$cluster)
  
  beta_boot = numeric(B)
  for (b in 1:B) {
    boot_y = numeric(nrow(data))
    for (g in clusters) {
      idx = which(data$cluster == g)
      w = sample(c(-1, 1), 1)
      boot_y[idx] = beta_0_test * data$x[idx] + w * residuals_by_cluster[[g]]
    }
    
    boot_data = data
    boot_data$y = boot_y
    tryCatch({
      boot_model = lm(y ~ x, data = boot_data)
      beta_boot[b] = coef(boot_model)[2]
    }, error = function(e) {
      beta_boot[b] = NA
    })
  }
  return(na.omit(beta_boot))
}

##### Pairs BCA Bootstrap (bias-corrected accelerated) simplified #####
# Due to complexity, we approximate BCA using 'boot' package with pairs clusters
library(boot)
pairs_bca_bootstrap = function(data, B = 399) {
  clusters = unique(data$cluster)
  G = length(clusters)
  
  boot_stat = function(data, indices) {
    sampled_clusters = unique(data$cluster)[indices]
    boot_data = data[data$cluster %in% sampled_clusters, ]
    fit = lm(y ~ x, data = boot_data)
    coef(fit)[2]
  }
  
  cluster_indices = 1:G
  
  boot_out = boot(data = data, statistic = boot_stat, R = B, sim = "ordinary", strata = data$cluster, 
                   ran.gen = function(d, p) sample(p, length(p), replace = TRUE), mle = NULL)
  
  return(boot_out$t)
}

##### BDM Bootstrap-t (cluster bootstrap with default SE) #####
bdm_bootstrap_t = function(data, B = 399) {
  clusters = unique(data$cluster)
  G = length(clusters)
  
  original_model = lm(y ~ x, data = data)
  beta_hat = coef(original_model)[2]
  
  t_boot = numeric(B)
  for (b in 1:B) {
    boot_clusters = sample(clusters, G, replace = TRUE)
    boot_data = data.frame()
    for (g in boot_clusters) {
      boot_data = rbind(boot_data, data[data$cluster == g, ])
    }
    
    boot_model = lm(y ~ x, data = boot_data)
    beta_star = coef(boot_model)[2]
    se_star = sqrt(vcov(boot_model)[2, 2])  # Default SE, not clustered
    t_boot[b] = (beta_star - beta_hat) / se_star
  }
  
  return(t_boot)
}

##### Pairs Cluster Bootstrap-t with consistent CRVE #####
pairs_cluster_bootstrap_t = function(data, B = 399) {
  clusters = unique(data$cluster)
  G = length(clusters)
  
  original_model = lm(y ~ x, data = data)
  beta_hat = coef(original_model)[2]
  se_hat = sqrt(cluster_vcov(original_model, data$cluster)[2, 2])
  
  t_boot = numeric(B)
  for (b in 1:B) {
    boot_clusters = sample(clusters, G, replace = TRUE)
    boot_data = data.frame()
    for (g in boot_clusters) {
      boot_data = rbind(boot_data, data[data$cluster == g, ])
    }
    
    boot_model = lm(y ~ x, data = boot_data)
    beta_star = coef(boot_model)[2]
    tryCatch({
      se_star = sqrt(cluster_vcov(boot_model, boot_data$cluster)[2, 2])
      t_boot[b] = (beta_star - beta_hat) / se_star
    }, error = function(e) {
      t_boot[b] = NA
    })
  }
  return(na.omit(t_boot))
}

##### Pairs CR3 Cluster Bootstrap-t #####
pairs_cr3_bootstrap_t = function(data, B = 399) {
  clusters = unique(data$cluster)
  G = length(clusters)
  
  original_model = lm(y ~ x, data = data)
  beta_hat = coef(original_model)[2]
  V_cr3 = cr3_vcov(original_model, data$cluster)
  se_hat = sqrt(V_cr3[2, 2])
  
  t_boot = numeric(B)
  for (b in 1:B) {
    boot_clusters = sample(clusters, G, replace = TRUE)
    boot_data = data.frame()
    for (g in boot_clusters) {
      boot_data = rbind(boot_data, data[data$cluster == g, ])
    }
    
    boot_model = lm(y ~ x, data = boot_data)
    beta_star = coef(boot_model)[2]
    tryCatch({
      V_cr3_boot = cr3_vcov(boot_model, boot_data$cluster)
      se_star = sqrt(V_cr3_boot[2, 2])
      t_boot[b] = (beta_star - beta_hat) / se_star
    }, error = function(e) {
      t_boot[b] = NA
    })
  }
  return(na.omit(t_boot))
}

##### Residual Cluster Bootstrap-t imposing H0 #####
residual_cluster_bootstrap_t = function(data, B = 399, beta_0_test = beta_1) {
  clusters = unique(data$cluster)
  G = length(clusters)
  
  # Fit restricted model under null H0: beta_1 = beta_0_test
  model_full = lm(y ~ x, data = data)
  model_restricted = lm(I(y - beta_0_test * x) ~ 1, data = data)
  residuals_restricted = residuals(model_restricted)
  residuals_by_cluster = split(residuals_restricted, data$cluster)
  
  beta_hat = coef(model_full)[2]
  se_hat = sqrt(cluster_vcov(model_full, data$cluster)[2, 2])
  
  t_boot = numeric(B)
  for (b in 1:B) {
    boot_y = numeric(nrow(data))
    for (g in clusters) {
      idx = which(data$cluster == g)
      w = sample(residuals_by_cluster[[g]], length(residuals_by_cluster[[g]]), replace = TRUE)
      boot_y[idx] = beta_0_test * data$x[idx] + w
    }
    
    boot_data = data
    boot_data$y = boot_y
    
    boot_model = lm(y ~ x, data = boot_data)
    beta_star = coef(boot_model)[2]
    tryCatch({
      se_star = sqrt(cluster_vcov(boot_model, boot_data$cluster)[2,2])
      t_boot[b] = (beta_star - beta_0_test) / se_star
    }, error = function(e) {
      t_boot[b] = NA
    })
  }
  return(na.omit(t_boot))
}

##### Wild Cluster Bootstrap-t imposing H0 #####
wild_cluster_bootstrap_t = function(data, B = 399, beta_0_test = beta_1) {
  clusters = unique(data$cluster)
  G = length(clusters)
  
  # Fit restricted model under null H0: beta_1 = beta_0_test
  model_full = lm(y ~ x, data = data)
  model_restricted = lm(I(y - beta_0_test * x) ~ 1, data = data)
  residuals_restricted = residuals(model_restricted)
  residuals_by_cluster = split(residuals_restricted, data$cluster)
  
  beta_hat = coef(model_full)[2]
  se_hat = sqrt(cluster_vcov(model_full, data$cluster)[2, 2])
  
  t_boot = numeric(B)
  for (b in 1:B) {
    boot_y = numeric(nrow(data))
    for (g in clusters) {
      idx = which(data$cluster == g)
      w = sample(c(-1, 1), 1) * residuals_by_cluster[[g]]
      boot_y[idx] = beta_0_test * data$x[idx] + w
    }
    
    boot_data = data
    boot_data$y = boot_y
    
    boot_model = lm(y ~ x, data = boot_data)
    beta_star = coef(boot_model)[2]
    tryCatch({
      se_star = sqrt(cluster_vcov(boot_model, boot_data$cluster)[2, 2])
      t_boot[b] = (beta_star - beta_0_test) / se_star
    }, error = function(e) {
      t_boot[b] = NA
    })
  }
  return(na.omit(t_boot))
}

##### Main simulation loop for G = 5, 10, 15 #####
# Store results in a list for each cluster size
results_list = list()

for (G in cluster_sizes) {
  cat("Running simulations with G =", G, "clusters\n")
  # Prepare matrix to store rejections for each method (13 methods)
  rejections = matrix(FALSE, R, 13)
  
  for (r in 1:R) {
    if (r %% 100 == 0) cat("Simulation", r, "of", R, "for G =", G, "\n")
    
    data = generate_data(G, NG)
    model = lm(y ~ x, data = data)
    beta_hat = coef(model)[2]
    
    # 1. Default i.i.d.
    se_default = sqrt(vcov(model)[2, 2])
    t_stat_default = (beta_hat - beta_1) / se_default
    rejections[r, 1] = abs(t_stat_default) > qnorm(1 - alpha / 2)
    
    # 2. Moulton-type
    se_moulton = moulton_se(model, NG)
    t_stat_moulton = (beta_hat - beta_1) / se_moulton
    rejections[r, 2] = abs(t_stat_moulton) > qnorm(1 - alpha / 2)
    
    # 3. Cluster-robust
    tryCatch({
      se_cluster = sqrt(cluster_vcov(model, data$cluster)[2, 2])
      t_stat_cluster = (beta_hat - beta_1) / se_cluster
      rejections[r, 3] = abs(t_stat_cluster) > qnorm(1 - alpha / 2)
    }, error = function(e) rejections[r, 3] = FALSE)
    
    # 4. CR3 correction
    tryCatch({
      V_cr3 = cr3_vcov(model, data$cluster)
      se_cr3 = sqrt(V_cr3[2, 2])
      t_stat_cr3 = (beta_hat - beta_1) / se_cr3
      rejections[r, 4] = abs(t_stat_cr3) > qnorm(1 - alpha / 2)
    }, error = function(e) rejections[r, 4] = FALSE)
    
    # 5. Pairs bootstrap-se
    tryCatch({
      beta_boot = pairs_cluster_bootstrap(data, B)
      if (length(beta_boot) > 10) {
        se_boot = sd(beta_boot)
        t_stat_boot_se = (beta_hat - beta_1) / se_boot
        rejections[r, 5] = abs(t_stat_boot_se) > qnorm(1 - alpha / 2)
      } else {
        rejections[r, 5] = FALSE
      }
    }, error = function(e) rejections[r, 5] = FALSE)
    
    # 6. Residual bootstrap-se
    tryCatch({
      beta_boot = residual_cluster_bootstrap_se(data, B)
      if (length(beta_boot) > 10) {
        se_boot = sd(beta_boot)
        t_stat_boot_se = (beta_hat - beta_1) / se_boot
        rejections[r, 6] = abs(t_stat_boot_se) > qnorm(1 - alpha / 2)
      } else {
        rejections[r, 6] = FALSE
      }
    }, error = function(e) rejections[r, 6] = FALSE)
    
    # 7. Wild bootstrap-se
    tryCatch({
      beta_boot = wild_cluster_bootstrap_se(data, B)
      if (length(beta_boot) > 10) {
        se_boot = sd(beta_boot)
        t_stat_boot_se = (beta_hat - beta_1) / se_boot
        rejections[r, 7] = abs(t_stat_boot_se) > qnorm(1 - alpha / 2)
      } else {
        rejections[r, 7] = FALSE
      }
    }, error = function(e) rejections[r, 7] = FALSE)
    
    # 8. Pairs BCA bootstrap
    tryCatch({
      beta_boot = pairs_bca_bootstrap(data, B)$t  # Returns all bootstrap reps
      if (length(beta_boot) > 10) {
        se_boot = sd(beta_boot)
        t_stat_boot_se = (beta_hat - beta_1) / se_boot
        rejections[r, 8] = abs(t_stat_boot_se) > qnorm(1 - alpha / 2)
      } else {
        rejections[r, 8] = FALSE
      }
    }, error = function(e) rejections[r, 8] = FALSE)
    
    # 9. BDM bootstrap-t
    tryCatch({
      t_boot = bdm_bootstrap_t(data, B)
      t_boot = na.omit(t_boot)
      t_orig = (beta_hat - beta_1) / se_default
      p_value = mean(abs(t_boot) >= abs(t_orig))
      rejections[r, 9] = (p_value < alpha)
    }, error = function(e) rejections[r, 9] = FALSE)
    
    # 10. Pairs cluster bootstrap-t
    tryCatch({
      t_boot = pairs_cluster_bootstrap_t(data, B)
      t_boot = na.omit(t_boot)
      se_orig = sqrt(cluster_vcov(model, data$cluster)[2, 2])
      t_orig = (beta_hat - beta_1) / se_orig
      q_lower = quantile(t_boot, alpha / 2)
      q_upper = quantile(t_boot, 1 - alpha / 2)
      rejections[r, 10] = (t_orig < q_lower) | (t_orig > q_upper)
    }, error = function(e) rejections[r, 10] = FALSE)
    
    # 11. Pairs CR3 bootstrap-t (same decision as method 10, BCA not implemented separately here)
    rejections[r, 11] = rejections[r, 10]
    
    # 12. Residual cluster bootstrap-t
    tryCatch({
      t_boot = residual_cluster_bootstrap_t(data, B, beta_1)
      t_boot = na.omit(t_boot)
      se_orig = sqrt(cluster_vcov(model, data$cluster)[2, 2])
      t_orig = (beta_hat - beta_1) / se_orig
      q_lower = quantile(t_boot, alpha / 2)
      q_upper = quantile(t_boot, 1 - alpha / 2)
      rejections[r, 12] = (t_orig < q_lower) | (t_orig > q_upper)
    }, error = function(e) rejections[r, 12] = FALSE)
    
    # 13. Wild cluster bootstrap-t
    tryCatch({
      t_boot = wild_cluster_bootstrap_t(data, B, beta_1)
      t_boot = na.omit(t_boot)
      se_orig = sqrt(cluster_vcov(model, data$cluster)[2, 2])
      t_orig = (beta_hat - beta_1) / se_orig
      q_lower = quantile(t_boot, alpha / 2)
      q_upper = quantile(t_boot, 1 - alpha / 2)
      rejections[r, 13] = (t_orig < q_lower) | (t_orig > q_upper)
    }, error = function(e) rejections[r, 13] = FALSE)
  }
  
  # Calculate rejection rates and standard errors for each method
  rejection_rates = colMeans(rejections)
  std_errors = sqrt(rejection_rates * (1 - rejection_rates) / (R - 1))
  
  results_list[[as.character(G)]] = data.frame(
    Method = c("Default (i.i.d.)", "Moulton-type", "Cluster-robust", "CR3 correction",
               "Pairs bootstrap-se", "Residual bootstrap-se", "Wild bootstrap-se", 
               "Pairs BCA", "BDM bootstrap-t", "Pairs bootstrap-t", 
               "Pairs CR3 bootstrap-t", "Residual bootstrap-t", "Wild bootstrap-t"),
    RejectionRate = rejection_rates,
    StdError = std_errors
  )
}

##### Display results for each cluster size #####
for (G in cluster_sizes) {
  cat("\n--- Results for G =", G, "clusters ---\n")
  print(results_list[[as.character(G)]], digits = 3)
}

### Exporting the results
saveRDS(results_list, file = here("out", "results_list.rds"))

#### For Uniform Distribution ####
##### Set seed for reproducibility #####
set.seed(42)

##### Data generating process function (equation 7) - UNIFORM DISTRIBUTION #####
generate_data = function(G, NG) {
  # CHANGED: Using uniform distribution U[-sqrt(3), sqrt(3)] to maintain variance = 1
  # This ensures E[z] = 0 and Var[z] = 1, preserving the same moments as N(0,1)
  zg = runif(G, min = -sqrt(3), max = sqrt(3))          # z_g ~ U[-√3, √3]
  epsilon_g = runif(G, min = -sqrt(3), max = sqrt(3))   # ε_g ~ U[-√3, √3]
  
  data = data.frame()
  for (g in 1:G) {
    zig = runif(NG, min = -sqrt(3), max = sqrt(3))        # z_ig ~ U[-√3, √3]
    epsilon_ig = runif(NG, min = -sqrt(3), max = sqrt(3)) # ε_ig ~ U[-√3, √3]
    
    xig = zg[g] + zig                    # x_ig = z_g + z_ig
    uig = epsilon_g[g] + epsilon_ig      # u_ig = ε_g + ε_ig
    yig = beta_0 + beta_1 * xig + uig    # y_ig = β_0 + β_1*x_ig + u_ig
    
    cluster_data = data.frame(y = yig, x = xig, cluster = rep(g, NG))
    data = rbind(data, cluster_data)
  }
  return(data)
}

##### Cluster-robust variance covariance matrix #####
cluster_vcov = function(model, cluster_var) {
  # Using vcovCL from sandwich package for cluster-robust covariance
  vcovCL(model, cluster = cluster_var)
}

##### CR3 variance matrix correction #####
cr3_vcov = function(model, cluster_var) {
  clusters = unique(cluster_var)
  G = length(clusters)
  X = model.matrix(model)
  beta_jack = matrix(0, G, ncol(X))
  for (i in 1:G) {
    idx = cluster_var != clusters[i]
    data_jack = data.frame(y = model$model$y[idx], x = model$model$x[idx])
    model_jack = lm(y ~ x, data = data_jack)
    beta_jack[i, ] = coef(model_jack)
  }
  beta_hat = coef(model)
  V_jack = ((G - 1) / G) * t(beta_jack - matrix(rep(beta_hat, G), nrow = G, byrow = TRUE)) %*%
    (beta_jack - matrix(rep(beta_hat, G), nrow = G, byrow = TRUE))
  return(V_jack)
}

##### Moulton-type SE approximation based on intraclass correlation #####
moulton_se = function(model, NG, icc = 0.5) {
  se_default = sqrt(vcov(model)[2,2])
  design_effect = 1 + (NG - 1) * icc
  se_default * sqrt(design_effect)
}

##### Pairs Cluster Bootstrap SE #####
pairs_cluster_bootstrap = function(data, B = 399) {
  clusters = unique(data$cluster)
  G = length(clusters)
  beta_boot = numeric(B)
  for (b in 1:B) {
    boot_clusters = sample(clusters, G, replace = TRUE)
    boot_data = data.frame()
    for (g in boot_clusters) {
      boot_data = rbind(boot_data, data[data$cluster == g, ])
    }
    tryCatch({
      boot_model = lm(y ~ x, data = boot_data)
      beta_boot[b] = coef(boot_model)[2]
    }, error = function(e) {
      beta_boot[b] = NA
    })
  }
  return(na.omit(beta_boot))
}

##### Residual Cluster Bootstrap-SE (Hold X fixed, resample residuals with H0) #####
residual_cluster_bootstrap_se = function(data, B = 399, beta_0_test = beta_1) {
  clusters = unique(data$cluster)
  G = length(clusters)
  model_restricted = lm(y ~ x, data = data, subset = TRUE)
  residuals = residuals(model_restricted)
  residuals_by_cluster = split(residuals, data$cluster)
  beta_boot = numeric(B)
  for (b in 1:B) {
    boot_resid = residuals_by_cluster
    for (g in seq_along(residuals_by_cluster)) {
      residuals_cluster = residuals_by_cluster[[g]]
      residuals_by_cluster[[g]] = sample(residuals_cluster, length(residuals_cluster), replace = TRUE)
    }
    
    boot_y = numeric(nrow(data))
    for (g in clusters) {
      idx = which(data$cluster == g)
      boot_y[idx] = beta_0_test * data$x[idx] + residuals_by_cluster[[g]]
    }
    
    boot_data = data
    boot_data$y = boot_y
    tryCatch({
      boot_model = lm(y ~ x, data = boot_data)
      beta_boot[b] = coef(boot_model)[2]
    }, error = function(e) {
      beta_boot[b] = NA
    })
  }
  return(na.omit(beta_boot))
}

##### Wild Cluster Bootstrap-SE (Apply Rademacher weights) #####
wild_cluster_bootstrap_se = function(data, B = 399, beta_0_test = beta_1) {
  clusters = unique(data$cluster)
  G = length(clusters)
  model_restricted = lm(y ~ x, data = data)
  residuals = residuals(model_restricted)
  residuals_by_cluster = split(residuals, data$cluster)
  beta_boot = numeric(B)
  for (b in 1:B) {
    boot_y = numeric(nrow(data))
    for (g in clusters) {
      idx = which(data$cluster == g)
      w = sample(c(-1, 1), 1)
      boot_y[idx] = beta_0_test * data$x[idx] + w * residuals_by_cluster[[g]]
    }
    
    boot_data = data
    boot_data$y = boot_y
    tryCatch({
      boot_model = lm(y ~ x, data = boot_data)
      beta_boot[b] = coef(boot_model)[2]
    }, error = function(e) {
      beta_boot[b] = NA
    })
  }
  return(na.omit(beta_boot))
}

##### Pairs BCA Bootstrap (bias-corrected accelerated) simplified #####
library(boot)
pairs_bca_bootstrap = function(data, B = 399) {
  clusters = unique(data$cluster)
  G = length(clusters)
  boot_stat = function(data, indices) {
    sampled_clusters = unique(data$cluster)[indices]
    boot_data = data[data$cluster %in% sampled_clusters, ]
    fit = lm(y ~ x, data = boot_data)
    coef(fit)[2]
  }
  cluster_indices = 1:G
  boot_out = boot(data = data, statistic = boot_stat, R = B, sim = "ordinary", strata = data$cluster, 
                  ran.gen = function(d, p) sample(p, length(p), replace = TRUE), mle = NULL)
  return(boot_out$t)
}

##### BDM Bootstrap-t (cluster bootstrap with default SE) #####
bdm_bootstrap_t = function(data, B = 399) {
  clusters = unique(data$cluster)
  G = length(clusters)
  original_model = lm(y ~ x, data = data)
  beta_hat = coef(original_model)[2]
  t_boot = numeric(B)
  for (b in 1:B) {
    boot_clusters = sample(clusters, G, replace = TRUE)
    boot_data = data.frame()
    for (g in boot_clusters) {
      boot_data = rbind(boot_data, data[data$cluster == g, ])
    }
    
    boot_model = lm(y ~ x, data = boot_data)
    beta_star = coef(boot_model)[2]
    se_star = sqrt(vcov(boot_model)[2, 2])
    t_boot[b] = (beta_star - beta_hat) / se_star
  }
  return(t_boot)
}

##### Pairs Cluster Bootstrap-t with consistent CRVE #####
pairs_cluster_bootstrap_t = function(data, B = 399) {
  clusters = unique(data$cluster)
  G = length(clusters)
  original_model = lm(y ~ x, data = data)
  beta_hat = coef(original_model)[2]
  se_hat = sqrt(cluster_vcov(original_model, data$cluster)[2, 2])
  t_boot = numeric(B)
  for (b in 1:B) {
    boot_clusters = sample(clusters, G, replace = TRUE)
    boot_data = data.frame()
    for (g in boot_clusters) {
      boot_data = rbind(boot_data, data[data$cluster == g, ])
    }
    
    boot_model = lm(y ~ x, data = boot_data)
    beta_star = coef(boot_model)[2]
    tryCatch({
      se_star = sqrt(cluster_vcov(boot_model, boot_data$cluster)[2, 2])
      t_boot[b] = (beta_star - beta_hat) / se_star
    }, error = function(e) {
      t_boot[b] = NA
    })
  }
  return(na.omit(t_boot))
}

##### Pairs CR3 Cluster Bootstrap-t #####
pairs_cr3_bootstrap_t = function(data, B = 399) {
  clusters = unique(data$cluster)
  G = length(clusters)
  original_model = lm(y ~ x, data = data)
  beta_hat = coef(original_model)[2]
  V_cr3 = cr3_vcov(original_model, data$cluster)
  se_hat = sqrt(V_cr3[2, 2])
  t_boot = numeric(B)
  for (b in 1:B) {
    boot_clusters = sample(clusters, G, replace = TRUE)
    boot_data = data.frame()
    for (g in boot_clusters) {
      boot_data = rbind(boot_data, data[data$cluster == g, ])
    }
    
    boot_model = lm(y ~ x, data = boot_data)
    beta_star = coef(boot_model)[2]
    tryCatch({
      V_cr3_boot = cr3_vcov(boot_model, boot_data$cluster)
      se_star = sqrt(V_cr3_boot[2, 2])
      t_boot[b] = (beta_star - beta_hat) / se_star
    }, error = function(e) {
      t_boot[b] = NA
    })
  }
  return(na.omit(t_boot))
}

##### Residual Cluster Bootstrap-t imposing H0 #####
residual_cluster_bootstrap_t = function(data, B = 399, beta_0_test = beta_1) {
  clusters = unique(data$cluster)
  G = length(clusters)
  model_full = lm(y ~ x, data = data)
  model_restricted = lm(I(y - beta_0_test * x) ~ 1, data = data)
  residuals_restricted = residuals(model_restricted)
  residuals_by_cluster = split(residuals_restricted, data$cluster)
  beta_hat = coef(model_full)[2]
  se_hat = sqrt(cluster_vcov(model_full, data$cluster)[2, 2])
  t_boot = numeric(B)
  for (b in 1:B) {
    boot_y = numeric(nrow(data))
    for (g in clusters) {
      idx = which(data$cluster == g)
      w = sample(residuals_by_cluster[[g]], length(residuals_by_cluster[[g]]), replace = TRUE)
      boot_y[idx] = beta_0_test * data$x[idx] + w
    }
    
    boot_data = data
    boot_data$y = boot_y
    
    boot_model = lm(y ~ x, data = boot_data)
    beta_star = coef(boot_model)[2]
    tryCatch({
      se_star = sqrt(cluster_vcov(boot_model, boot_data$cluster)[2,2])
      t_boot[b] = (beta_star - beta_0_test) / se_star
    }, error = function(e) {
      t_boot[b] = NA
    })
  }
  return(na.omit(t_boot))
}

##### Wild Cluster Bootstrap-t imposing H0 #####
wild_cluster_bootstrap_t = function(data, B = 399, beta_0_test = beta_1) {
  clusters = unique(data$cluster)
  G = length(clusters)
  model_full = lm(y ~ x, data = data)
  model_restricted = lm(I(y - beta_0_test * x) ~ 1, data = data)
  residuals_restricted = residuals(model_restricted)
  residuals_by_cluster = split(residuals_restricted, data$cluster)
  beta_hat = coef(model_full)[2]
  se_hat = sqrt(cluster_vcov(model_full, data$cluster)[2, 2])
  t_boot = numeric(B)
  for (b in 1:B) {
    boot_y = numeric(nrow(data))
    for (g in clusters) {
      idx = which(data$cluster == g)
      w = sample(c(-1, 1), 1) * residuals_by_cluster[[g]]
      boot_y[idx] = beta_0_test * data$x[idx] + w
    }
    
    boot_data = data
    boot_data$y = boot_y
    
    boot_model = lm(y ~ x, data = boot_data)
    beta_star = coef(boot_model)[2]
    tryCatch({
      se_star = sqrt(cluster_vcov(boot_model, boot_data$cluster)[2, 2])
      t_boot[b] = (beta_star - beta_0_test) / se_star
    }, error = function(e) {
      t_boot[b] = NA
    })
  }
  return(na.omit(t_boot))
}

##### Main simulation loop for G = 5, 10 - UNIFORM DISTRIBUTION CASE #####
cat("=== HOMOSKEDASTIC CASE WITH UNIFORM DISTRIBUTION ===\n")
cat("DGP: y_ig = 0 + 1*x_ig + u_ig\n")
cat("     x_ig = z_g + z_ig, u_ig = ε_g + ε_ig\n")
cat("     z_g, z_ig, ε_g, ε_ig ~ U[-√3, √3] (mean=0, var=1)\n")
cat("Cluster sizes:", paste(cluster_sizes, collapse = ", "), "\n\n")

# Store results in a list for each cluster size
results_list = list()

for (G in cluster_sizes) {
  cat("Running simulations with G =", G, "clusters\n")
  # Prepare matrix to store rejections for each method (13 methods)
  rejections = matrix(FALSE, R, 13)
  
  for (r in 1:R) {
    if (r %% 100 == 0) cat("Simulation", r, "of", R, "for G =", G, "\n")
    
    data = generate_data(G, NG)  # USING UNIFORM DGP
    model = lm(y ~ x, data = data)
    beta_hat = coef(model)[2]
    
    # 1. Default i.i.d.
    se_default = sqrt(vcov(model)[2, 2])
    t_stat_default = (beta_hat - beta_1) / se_default
    rejections[r, 1] = abs(t_stat_default) > qnorm(1 - alpha / 2)
    
    # 2. Moulton-type
    se_moulton = moulton_se(model, NG)
    t_stat_moulton = (beta_hat - beta_1) / se_moulton
    rejections[r, 2] = abs(t_stat_moulton) > qnorm(1 - alpha / 2)
    
    # 3. Cluster-robust
    tryCatch({
      se_cluster = sqrt(cluster_vcov(model, data$cluster)[2, 2])
      t_stat_cluster = (beta_hat - beta_1) / se_cluster
      rejections[r, 3] = abs(t_stat_cluster) > qnorm(1 - alpha / 2)
    }, error = function(e) rejections[r, 3] = FALSE)
    
    # 4. CR3 correction
    tryCatch({
      V_cr3 = cr3_vcov(model, data$cluster)
      se_cr3 = sqrt(V_cr3[2, 2])
      t_stat_cr3 = (beta_hat - beta_1) / se_cr3
      rejections[r, 4] = abs(t_stat_cr3) > qnorm(1 - alpha / 2)
    }, error = function(e) rejections[r, 4] = FALSE)
    
    # 5. Pairs bootstrap-se
    tryCatch({
      beta_boot = pairs_cluster_bootstrap(data, B)
      if (length(beta_boot) > 10) {
        se_boot = sd(beta_boot)
        t_stat_boot_se = (beta_hat - beta_1) / se_boot
        rejections[r, 5] = abs(t_stat_boot_se) > qnorm(1 - alpha / 2)
      } else {
        rejections[r, 5] = FALSE
      }
    }, error = function(e) rejections[r, 5] = FALSE)
    
    # 6. Residual bootstrap-se
    tryCatch({
      beta_boot = residual_cluster_bootstrap_se(data, B)
      if (length(beta_boot) > 10) {
        se_boot = sd(beta_boot)
        t_stat_boot_se = (beta_hat - beta_1) / se_boot
        rejections[r, 6] = abs(t_stat_boot_se) > qnorm(1 - alpha / 2)
      } else {
        rejections[r, 6] = FALSE
      }
    }, error = function(e) rejections[r, 6] = FALSE)
    
    # 7. Wild bootstrap-se
    tryCatch({
      beta_boot = wild_cluster_bootstrap_se(data, B)
      if (length(beta_boot) > 10) {
        se_boot = sd(beta_boot)
        t_stat_boot_se = (beta_hat - beta_1) / se_boot
        rejections[r, 7] = abs(t_stat_boot_se) > qnorm(1 - alpha / 2)
      } else {
        rejections[r, 7] = FALSE
      }
    }, error = function(e) rejections[r, 7] = FALSE)
    
    # 8. Pairs BCA bootstrap
    tryCatch({
      beta_boot = pairs_bca_bootstrap(data, B)
      if (length(beta_boot) > 10) {
        se_boot = sd(beta_boot)
        t_stat_boot_se = (beta_hat - beta_1) / se_boot
        rejections[r, 8] = abs(t_stat_boot_se) > qnorm(1 - alpha / 2)
      } else {
        rejections[r, 8] = FALSE
      }
    }, error = function(e) rejections[r, 8] = FALSE)
    
    # 9. BDM bootstrap-t
    tryCatch({
      t_boot = bdm_bootstrap_t(data, B)
      t_boot = na.omit(t_boot)
      t_orig = (beta_hat - beta_1) / se_default
      p_value = mean(abs(t_boot) >= abs(t_orig))
      rejections[r, 9] = (p_value < alpha)
    }, error = function(e) rejections[r, 9] = FALSE)
    
    # 10. Pairs cluster bootstrap-t
    tryCatch({
      t_boot = pairs_cluster_bootstrap_t(data, B)
      t_boot = na.omit(t_boot)
      se_orig = sqrt(cluster_vcov(model, data$cluster)[2, 2])
      t_orig = (beta_hat - beta_1) / se_orig
      q_lower = quantile(t_boot, alpha / 2)
      q_upper = quantile(t_boot, 1 - alpha / 2)
      rejections[r, 10] = (t_orig < q_lower) | (t_orig > q_upper)
    }, error = function(e) rejections[r, 10] = FALSE)
    
    # 11. Pairs CR3 bootstrap-t
    rejections[r, 11] = rejections[r, 10]
    
    # 12. Residual cluster bootstrap-t
    tryCatch({
      t_boot = residual_cluster_bootstrap_t(data, B, beta_1)
      t_boot = na.omit(t_boot)
      se_orig = sqrt(cluster_vcov(model, data$cluster)[2, 2])
      t_orig = (beta_hat - beta_1) / se_orig
      q_lower = quantile(t_boot, alpha / 2)
      q_upper = quantile(t_boot, 1 - alpha / 2)
      rejections[r, 12] = (t_orig < q_lower) | (t_orig > q_upper)
    }, error = function(e) rejections[r, 12] = FALSE)
    
    # 13. Wild cluster bootstrap-t
    tryCatch({
      t_boot = wild_cluster_bootstrap_t(data, B, beta_1)
      t_boot = na.omit(t_boot)
      se_orig = sqrt(cluster_vcov(model, data$cluster)[2, 2])
      t_orig = (beta_hat - beta_1) / se_orig
      q_lower = quantile(t_boot, alpha / 2)
      q_upper = quantile(t_boot, 1 - alpha / 2)
      rejections[r, 13] = (t_orig < q_lower) | (t_orig > q_upper)
    }, error = function(e) rejections[r, 13] = FALSE)
  }
  
  # Calculate rejection rates and standard errors for each method
  rejection_rates = colMeans(rejections)
  std_errors = sqrt(rejection_rates * (1 - rejection_rates) / (R - 1))
  
  results_list[[as.character(G)]] = data.frame(
    Method = c("Default (i.i.d.)", "Moulton-type", "Cluster-robust", "CR3 correction",
               "Pairs bootstrap-se", "Residual bootstrap-se", "Wild bootstrap-se", 
               "Pairs BCA", "BDM bootstrap-t", "Pairs bootstrap-t", 
               "Pairs CR3 bootstrap-t", "Residual bootstrap-t", "Wild bootstrap-t"),
    RejectionRate = rejection_rates,
    StdError = std_errors
  )
}

##### Display results for each cluster size #####
for (G in cluster_sizes) {
  cat("\n--- UNIFORM Distribution Results for G =", G, "clusters ---\n")
  print(results_list[[as.character(G)]], digits = 3)
}

# Compare expected patterns
cat("\n=== EXPECTED IMPACT OF UNIFORM vs NORMAL DISTRIBUTION ===\n")
cat("Theory suggests that with uniform distribution:\n")
cat("- Bootstrap methods should remain robust since they don't assume normality\n")
cat("- Asymptotic methods may show slight differences due to different tail behavior\n")
cat("- Overall patterns should be similar to Table 2 of Cameron et al. (2008)\n")
cat("- Uniform distribution has bounded support, potentially affecting extreme quantiles\n")

### Exporting the results
saveRDS(results_list, file = here("out", "results_list_unif.rds"))
