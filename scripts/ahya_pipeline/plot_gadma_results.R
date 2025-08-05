# Title: Plotting gadma results
# Author: Katharine Prata
# Description: taking in combined and formatted groups and plotting the parameter
# values

library(tidyverse)
# tibble 3.2.1
library(scales)

# Converting dadi units to physical units constants
mu = 1.86*10^(-8)
#gen.small = 2
#gen.big = 5
L.genom = 390206788

# Functions
calc_snps <- function(gadma_results, proj) {
  skip_groups = character(length = 0L)
  gadma_results$unfil.snps = NA
  gadma_results$sfs.snps = NA
  for (i in 1:nrow(gadma_results)) {
    # effective sequence params
    if (gadma_results[i, 1] == "group1-group2") {
      cat("\nunfiltered snps for: out_QCbasic_noSites_group1_dadi_group2_dadi.pos.gz")
      if (proj == "proj") {
        gadma_results$sfs.snps[i] = 285376.87
      } else {
        gadma_results$sfs.snps[i] = 323336.46
      }
      gadma_results$unfil.snps[i] = 33523695
    } else if (gadma_results[i, 1] == "group1-Amil") {
      cat("\nunfiltered snps for: out_QCbasic_noSites_group1_dadi_Amil_dadi.pos.gz")
      if (proj == "proj") {
        gadma_results$sfs.snps[i] = 111991.19
      } else {
        gadma_results$sfs.snps[i] = 444001.5
      }
      gadma_results$unfil.snps[i]  = 32965107
    } else if (gadma_results[i, 1] == "group1-group3") {
      cat("\nunfiltered snps for: out_QCbasic_noSites_group1_dadi_group3_dadi.pos.gz")
      if (proj == "proj") {
        gadma_results$sfs.snps[i] = 70052.17
      } else {
        gadma_results$sfs.snps[i] = 339902.54
      }
      gadma_results$unfil.snps[i] = 33812349
    } else if (gadma_results[i, 1] == "group1-group4") {
      cat("\nunfiltered snps for: out_QCbasic_noSites_group1_dadi_group4_dadi.pos.gz")
      if (proj == "proj") {
        gadma_results$sfs.snps[i] = 76247.09
      } else {
        gadma_results$sfs.snps[i] = 380214.39
      }
      gadma_results$unfil.snps[i] = 34517876
    } else if (gadma_results[i, 1] == "group2-group3") {
      cat("\nunfiltered snps for: out_QCbasic_noSites_group2_dadi_group3_dadi.pos.gz")
      if (proj == "proj") {
        gadma_results$sfs.snps[i] = 80497.72
      } else {
        gadma_results$sfs.snps[i] = 371442.47
      }
      gadma_results$unfil.snps[i] = 34492471
    } else if (gadma_results[i, 1] == "group2-group4") {
      cat("\nunfiltered snps for: out_QCbasic_noSites_group2_dadi_group4_dadi.pos.gz")
      if (proj == "proj") {
        gadma_results$sfs.snps[i] = 84843.52
      } else {
        gadma_results$sfs.snps[i] = 413972.89
      }
      gadma_results$unfil.snps[i] = 36046531
    } else if (gadma_results[i, 1] == "group3-group4") {
      cat("\nunfiltered snps for: out_QCbasic_noSites_group3_dadi_group4_dadi.pos.gz")
      if (proj == "proj") {
        gadma_results$sfs.snps[i] = 92922.75
      } else {
        gadma_results$sfs.snps[i] = 421890.54
      }
      gadma_results$unfil.snps[i] = 35072420
    } else if (gadma_results[i, 1] == "group2-Amil") {
      cat("\nunfiltered snps for: out_QCbasic_noSites_group2_dadi_Amil_dadi.pos.gz")
      if (proj == "proj") {
        gadma_results$sfs.snps[i] = 118359.4
      } else {
        gadma_results$sfs.snps[i] = 472699.58
      }
      gadma_results$unfil.snps[i] = 34184556
    } else if (gadma_results[i, 1] == "group3-Amil") {
      cat("\nunfiltered snps for: out_QCbasic_noSites_group3_dadi_Amil_dadi.pos.gz")
      if (proj == "proj") {
        gadma_results$sfs.snps[i] = 127781.61
      } else {
        gadma_results$sfs.snps[i] = 500811.35
      }
      gadma_results$unfil.snps[i] = 34598003
    } else if (gadma_results[i, 1] == "group4-Amil") {
      cat("\nunfiltered snps for: out_QCbasic_noSites_group4_dadi_Amil_dadi.pos.gz")
      if (proj == "proj") {
        gadma_results$sfs.snps[i] = 139754.46
      } else {
        gadma_results$sfs.snps[i] = 556220.74
      }
      gadma_results$unfil.snps[i] = 35256168
    } else {
      skip_groups = append(skip_groups, paste(gadma_results[i, 1]))
      cat("\nDon't have params for:", paste(gadma_results[i, 1]), "skipping")
      next
    }
  }
  skip_groups = unique(skip_groups)
  cat("\nNo values for:", paste(skip_groups), "removing from dataframe")
  gadma_results <- gadma_results[!is.na(gadma_results$unfil.snps), ]
  return(gadma_results)
}
calc_L.eff <- function(L.genom, snps) {
  L.eff <- L.genom * snps
  return(L.eff)
}
calc_nref <- function(theta, mu, L.eff) {
  nref <- theta / (4*mu*L.eff)
  return(nref)
}
calc_Tgen <- function(nref, t) {
  Tgen <-   ifelse(t == 0, 0,
                   2 * nref * t)
  return(Tgen)
}
calc_N <- function(nref, n) {
  N <- nref * n
  return(N)
}
calc_rM <- function(nref, m) {
  rM <- ifelse(m == 0, 0,
         m/(2*nref))
  return(rM)
}
calc_M <- function(N, rM) {
  M <- ifelse(rM == 0, 0,
               rM * N)
  return(M)
}
convert_params <- function(gadma_results, model, proj, L.genom, mu) {
 
  gadma_results <- calc_snps(gadma_results, proj)

  # Step 1. Calculate effective L
  snps <- (gadma_results$sfs.snps/gadma_results$unfil.snps)
  gadma_results$L.effective <- calc_L.eff(L.genom, snps)
  
  # Step 2. Calculate Nref (theta = 4Nref*mu*L)
  gadma_results$nref <- calc_nref(gadma_results$theta, mu,
                                  gadma_results$L.effective)
  
  # Divergence time periods
  gadma_results$T1gen <- calc_Tgen(gadma_results$nref, gadma_results$t1)
  # T2
  gadma_results$T2gen <- calc_Tgen(gadma_results$nref, gadma_results$t2)
  gadma_results$T1gen <- gadma_results$T1gen + gadma_results$T2gen
  
  # Population sizes after divergence
  gadma_results$N1.0 <- calc_N(gadma_results$nref, gadma_results$nu_1)
  gadma_results$N2.0 <- calc_N(gadma_results$nref, gadma_results$nu_2)
  
  # T1 population sizes
  gadma_results$N1.T1 <- calc_N(gadma_results$nref, gadma_results$nu11)
  gadma_results$N2.T1 <- calc_N(gadma_results$nref, gadma_results$nu12)
  
  # T2 population sizes
  gadma_results$N1.T2 <- calc_N(gadma_results$nref, gadma_results$nu21)
  gadma_results$N2.T2 <- calc_N(gadma_results$nref, gadma_results$nu22)
  
  if (model != "1het") {
    # T1 migration rates
    gadma_results$rM12.T1 <- calc_rM(gadma_results$nref, gadma_results$m1_12)
    gadma_results$rM21.T1 <- calc_rM(gadma_results$nref, gadma_results$m1_21)
    # T1 migrants
    gadma_results$M12.T1 <- calc_M(gadma_results$N1.T1, gadma_results$rM12.T1) 
    gadma_results$M21.T1 <- calc_M(gadma_results$N2.T1, gadma_results$rM21.T1) 
    # T2 migration rates
    gadma_results$rM12.T2 <- calc_rM(gadma_results$nref, gadma_results$m2_12)
    gadma_results$rM21.T2 <- calc_rM(gadma_results$nref, gadma_results$m2_21)
    # T2 migrants
    gadma_results$M12.T2 <- calc_M(gadma_results$N1.T2, gadma_results$rM12.T2)
    gadma_results$M21.T2 <- calc_M(gadma_results$N2.T1, gadma_results$rM21.T2)
  }
  
  if (model == "2het" | model == "1het") {
    # T1 migration rates
    gadma_results$rMe12.T1 <- calc_rM(gadma_results$nref, gadma_results$me1_12)
    gadma_results$rMe21.T1 <- calc_rM(gadma_results$nref, gadma_results$me1_21)
    # T1 migrants
    gadma_results$Me12.T1 <- calc_M(gadma_results$N1.T1, gadma_results$rMe12.T1)
    gadma_results$Me21.T1 <- calc_M(gadma_results$N2.T1, gadma_results$rMe12.T1)
    # T2 migration rates
    gadma_results$rMe12.T2 <- calc_rM(gadma_results$nref, gadma_results$me2_12)
    gadma_results$rMe21.T2 <- calc_rM(gadma_results$nref, gadma_results$me2_12)
    # T2 migrants
    gadma_results$Me12.T2 <- calc_M(gadma_results$N1.T2, gadma_results$rMe12.T2)
    gadma_results$Me21.T2 <- calc_M(gadma_results$N2.T2, gadma_results$rMe21.T2)
  }
  return(gadma_results)
}

setwd("~/git/kp_dadi/scripts/")

# Import results ####
gadma_results <- read.delim("../results/gadma_proj_1het_results_combined_noerror.txt")
proj <- "proj"
model <- "1het"
optim <- paste(proj, model, sep = "_")

# Load confidence interval files ####
conf_groups <- c(
  "group1-group2", "group1-group3", "group1-group4",
  "group2-group3", "group2-group4", "group3-group4",
  "group1-Amil", "group2-Amil", "group3-Amil", "group4-Amil"
)
conf_files <- paste0("../results/", conf_groups, "_projected0.8_1het_confidence_intervals.txt")
all_conf_list <- lapply(seq_along(conf_files), function(i) {
  df <- read.delim(conf_files[i])
  df$Group <- conf_groups[i]
  df
})
all_conf <- do.call(rbind, all_conf_list)
all_conf$Optimised_thres_max <- all_conf$Optimised/10000
all_conf$Optimised_thres_min <- all_conf$Optimised/1000
all_conf$chosen_eps <- (all_conf$Optimised < 10 & all_conf$eps > all_conf$Optimised_thres_max & all_conf$eps < all_conf$Optimised_thres_min) |
      (all_conf$Optimised >= 10 & all_conf$eps == 0.01)
all_conf <- all_conf %>%
  group_by(Group, Parameter) %>%
  filter(if (any(chosen_eps)) chosen_eps else eps == max(eps)) %>%
  ungroup()
# Parameter order
order_params <- c("nu_1", "nu_2", "t1", "nu11", "nu12", "me1_12", "me1_21", "t2", "nu21", "nu22", "me2_12", "me2_21", "P1", "P2", "theta")

all_conf$Parameter <- factor(all_conf$Parameter, levels = order_params, ordered = TRUE)
all_conf <- all_conf[order(all_conf$Parameter), ]
all_conf$Optimised[all_conf$Group == "group2-group3"]

# Reformat gadma results so each parameter is a column ####
param_names <- gadma_results$ModelParameters[1] |>
  str_remove_all("[()]") |>
  str_split(",\\s*") |>
  unlist()

model_values_matrix <- gadma_results$ModelValues |>
  str_remove_all("[()]") |>
  str_split(",\\s*") |>
  map(~ as_tibble_row(set_names(as.numeric(.x), param_names))) |>
  list_rbind()

model_values_df <- as_tibble(model_values_matrix)
colnames(model_values_df) <- param_names

gadma_results <- bind_cols(gadma_results |> select(-ModelParameters, -ModelValues), model_values_df)

rm(model_values_df, model_values_matrix, param_names)

if (model == "") {
  num1 <- 4
  num2 <- 3
} else {
  num1 <- 5
  num2 <- 4
}
# Fix group comparison names
gadma_results <- gadma_results %>%
  mutate(Groups = str_split_fixed(Directory, "_", num1)[, num2]) %>%
  mutate(Groups = str_replace_all(Groups, "grp", "group")) %>%
  mutate(Groups = str_replace_all(Groups, "het", "group1-group2")) %>%
  mutate(Groups = str_replace_all(Groups, "-process.*", ""))

gadma_results$Groups <- as.factor(gadma_results$Groups)

gadma_results <- gadma_results %>% 
  select(Groups, everything(), -Directory)

str(gadma_results)
colnames(gadma_results)[colnames(gadma_results) == "Theta"] <- "theta"
# Remove results where theta could not be calculated or model messed up optimisation (e.g., resulted in unplottable model)
gadma_results$theta[gadma_results$theta == "--"] = NA
gadma_results <- na.omit(gadma_results)
gadma_results$theta <- as.numeric(gadma_results$theta)
gadma_results <- gadma_results[gadma_results$theta > 100 & gadma_results$theta < 1e6,]
if (optim == "full_1het" | optim == "full_2het") {
  gadma_results <- gadma_results[gadma_results$LogLikelihood < -2500,]
}

# Reshape to long format for plotting
gadma_long <- gadma_results %>%
  pivot_longer(cols = -c(Groups, Run, LogLikelihood),
               names_to = "Parameter",
               values_to = "Value")

# Normalise LogLikelihood for color scale (lower loglik = worse)
gadma_long <- gadma_long %>%
  mutate(
    loglike_scaled = rescale(LogLikelihood, to = c(0, 1)),  # 0 = low (red), 1 = high (green)
    Run = factor(Run)  # So lines group correctly
  )

# plots
for (group in levels(gadma_long$Groups)) {
  df_group = gadma_long[gadma_long$Groups == group,]
  p <- ggplot(df_group, aes(x = Parameter, y = Value)) +
    geom_line(aes(color = loglike_scaled), alpha = 0.4, linewidth = 0.6) +
    geom_point(aes(color = loglike_scaled), size = 2) +
    facet_wrap(~Parameter, scales = "free_y") +
    scale_color_gradient(low = "red", high = "green", name = "Log-Likelihood") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.text.y = element_text(angle = 0),
      strip.text.x = element_text(angle = 45, hjust = 0)
    ) +
    labs(
      title = group,
      subtitle = "Parameter Values Across Runs"
    )
  print(p)
}

# best models
for (group in unique(gadma_results$Groups)) {
  gadma_group <- gadma_results[gadma_results$Groups == group,]
  max_idx <- which.max(gadma_group$LogLikelihood)
  print(gadma_group[max_idx,])
}

# Convert params ####
gadma_results_pu <- convert_params(gadma_results,
                                   model = model,
                                   proj = proj,
                                   L.genom = L.genom,
                                   mu = mu)

# Confidence intervals ####
gadma_conf <- all_conf[,1:6]
conf_long <- gadma_conf %>%
  pivot_longer(cols = c("Optimised", "Lower_CI", "Upper_CI", "eps"),
               names_to = "stat",
               values_to = "value")
conf_long <- conf_long %>%
  mutate(param_stat = paste(Parameter, stat, sep = "_"))

conf_wide <- conf_long %>%
  select(Group, param_stat, value) %>%
  pivot_wider(names_from = param_stat, values_from = value)

conf_wide <- calc_snps(conf_wide, proj)
conf_wide$snps <- (conf_wide$sfs.snps/conf_wide$unfil.snps)
conf_wide$L.eff <- calc_L.eff(L.genom, conf_wide$snps)
conf_wide$nref <- calc_nref(conf_wide$theta_Optimised, mu, conf_wide$L.eff)
# Pop sizes
pop_cols <- c("nu_1_Optimised", "nu_1_Lower_CI", "nu_1_Upper_CI",
              "nu_2_Optimised", "nu_2_Lower_CI", "nu_2_Upper_CI",
              "nu11_Optimised", "nu11_Lower_CI", "nu11_Upper_CI",
              "nu12_Optimised", "nu12_Lower_CI", "nu12_Upper_CI",
              "nu21_Optimised", "nu21_Lower_CI", "nu21_Upper_CI",
              "nu22_Optimised", "nu22_Lower_CI", "nu22_Upper_CI")
for (col in pop_cols) {
  conf_wide[[col]] <- calc_N(conf_wide$nref, conf_wide[[col]])
}
# Tgen
t_cols <- c("t1_Optimised", "t1_Lower_CI", "t1_Upper_CI",
            "t2_Optimised", "t2_Lower_CI", "t2_Upper_CI")
for (col in t_cols) {
  conf_wide[[col]] <- calc_Tgen(conf_wide$nref, conf_wide[[col]])
}
# rM
m_cols <- c("me1_12_Optimised", "me1_12_Lower_CI","me1_12_Upper_CI",
            "me1_21_Optimised", "me1_21_Lower_CI","me1_21_Upper_CI",
            "me2_12_Optimised", "me2_12_Lower_CI","me2_12_Upper_CI",
            "me2_21_Optimised", "me2_21_Lower_CI","me2_21_Upper_CI")
for (col in m_cols) {
  conf_wide[[col]] <- calc_rM(conf_wide$nref, conf_wide[[col]])
}

me1_12_cols <- c("me1_12_Optimised", "me1_12_Lower_CI","me1_12_Upper_CI")
for (col in me1_12_cols) {
  conf_wide[[col]] <- calc_M(conf_wide$nu11_Optimised, conf_wide[[col]])
}
me1_21_cols <- c("me1_21_Optimised", "me1_21_Lower_CI","me1_21_Upper_CI")
for (col in me1_21_cols) {
  conf_wide[[col]] <- calc_M(conf_wide$nu12_Optimised, conf_wide[[col]])
}
me2_12_cols <- c("me2_12_Optimised", "me2_12_Lower_CI","me2_12_Upper_CI")
for (col in me2_12_cols) {
  conf_wide[[col]] <- calc_M(conf_wide$nu21_Optimised, conf_wide[[col]])
}
me2_21_cols <- c("me2_21_Optimised", "me2_21_Lower_CI","me2_21_Upper_CI")
for (col in me2_21_cols) {
  conf_wide[[col]] <- calc_M(conf_wide$nu22_Optimised, conf_wide[[col]])
}

# Get unique parameter names by removing the suffixes
param_names <- unique(str_replace(names(conf_wide), "_Optimised|_Lower_CI|_Upper_CI|_eps", ""))
param_names <- param_names[param_names != "Group"] # Remove Group if present
param_names <- param_names[param_names != ""]      # Remove any empty strings
other_cols <- c("unfil.snps", "sfs.snps", "snps", "L.eff", "nref")
param_names <- param_names[!param_names %in% other_cols]      # Remove any empty strings


# High number parameters
high_params <- c("nu_1", "nu_2", "t1", "t2", "nu11", "nu12", "nu21", "nu22")

# Build the full column names to modify
high_params_all <- unlist(lapply(high_params, function(p) {
  c(paste0(p, "_Optimised"), paste0(p, "_Lower_CI"), paste0(p, "_Upper_CI"))
}))

# Divide each column by 1000
conf_wide[high_params_all] <- lapply(conf_wide[high_params_all], function(x) x / 1000)

# Start here for formatted confidence intervals ####
combined_cols <- data.frame(Group = conf_wide$Group)

format_sig <- function(x) {
  ifelse(abs(x) < 1, signif(x, 2), signif(x, 4))
}

for (p in param_names) {
  opt_col <- paste0(p, "_Optimised")
  lower_col <- paste0(p, "_Lower_CI")
  upper_col <- paste0(p, "_Upper_CI")
  combined_col <- paste0(p)
  
  combined_cols[[combined_col]] <- paste0(
    format_sig(conf_wide[[opt_col]]),
    " (", format_sig(conf_wide[[lower_col]]),
    ", ", format_sig(conf_wide[[upper_col]]), ")"
  )
}

str(combined_cols)

# Change order here as well ####
if (model == "2het") {
  param_cols <- c(26,27,28,29,36:39,30,31,21,
                  40,47:50,41,42,22)
} else if (model == "1het") {
  param_cols <- c(23:25,28:31,26:27,17,32,35:38,33:34,18)
} else if (model == "") {
  param_cols <- c(21:23,28,29,26:27,24:25,30,33:36,31:32)
}

# Reshape to long format
gadma_pu_long <- gadma_results_pu %>%
  pivot_longer(
    cols = all_of(param_cols),
    names_to = "Physical.Parameter",
    values_to = "Physical.Value"
  )

parameter_order <- colnames(gadma_results_pu[,param_cols])

gadma_pu_long <- gadma_pu_long |>
  mutate(Physical.Parameter = factor(Physical.Parameter, levels = parameter_order))

# Plot
for (group in unique(gadma_pu_long$Groups)) {
  df_group = gadma_pu_long[gadma_pu_long$Groups == group,]
  # Filter to top 5 rows by LogLikelihood (highest = least negative)
  top5_runs <- df_group |>
    #dplyr::group_by(Run) |>
    #dplyr::summarise(LogLikelihood = mean(LogLikelihood), .groups = "drop") |>
    dplyr::arrange(desc(LogLikelihood)) |>
    dplyr::slice_head(n = 5)
  df_group <- df_group[df_group$Run %in% top5_runs$Run, ]
  p <- ggplot(df_group, aes(x = Physical.Parameter, y = Physical.Value)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  facet_wrap(~ Physical.Parameter, scales = "free_y") +
  scale_y_continuous(labels = function(x) {
      ifelse(abs(x) > 1000, scales::scientific_format()(x), x)
    }) +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  ylab("Parameter Value") +
  xlab("Parameter") +
  labs(
      title = group,
      subtitle = "Confidence intervals"
    ) +
    theme(
      strip.text = element_text(face = "bold"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_text(margin = margin(t = 10))
    )
  print(p)
  ggsave(paste0("../plots/gadma_params/", group, "_", optim,
                "_physical-param-confi.pdf"),
         units = "cm", width = 20, height = 20)
}

summarise_data <- function(gadma_results_pu, param_cols) {
  summary_list <- list()
  for (group in unique(gadma_results_pu$Groups)) {
    df_group <- gadma_results_pu[gadma_results_pu$Groups == group,]
    
    top_runs <- df_group |>
      arrange(desc(LogLikelihood)) |>
      slice_head(n = 3)
    
    # Calculate min and max for columns in thousands
    min_vals <- top_runs[, param_cols] |> summarise(across(everything(), min))
    max_vals <- top_runs[, param_cols] |> summarise(across(everything(), max))
    # dived by 1000 for generations and pop sizes
    desired_cols <- c("nref", "T1gen", "N1.0", "N2.0", "N1.T1", "N2.T1", "T2gen", "N1.T2", "N2.T2")
    existing_cols <- desired_cols[desired_cols %in% names(min_vals)]
    min_vals[existing_cols] <- lapply(min_vals[existing_cols], function(x) x / 1000)
    max_vals[existing_cols] <- lapply(max_vals[existing_cols], function(x) x / 1000)
    #if (min_vals$P1 > 0.5 & max_vals$P1 >0.5) {
      # Combine min and max as "min-max"
      #min_max <- Map(function(min, max) sprintf("%.1f–%.1f", min, max), min_vals, max_vals)
      #min_max_df <- as.data.frame(min_max)
      #names(min_max_df) <- names(min_vals)
    #}
    min_max_df <- rbind(round(min_vals,1), round(max_vals,1))
    
    min_max_df$Group <- group
    summary_list[[group]] <- min_max_df
  }
  
  summary_table <- bind_rows(summary_list) |> relocate(Group)
  
  custom_order <- c(
    "group1-group2", "group1-group3", "group1-group4",
    "group2-group3", "group2-group4", "group3-group4",
    "group1-Amil", "group2-Amil", "group3-Amil", "group4-Amil"
  )
  
  summary_table$Group <- factor(summary_table$Group, levels = custom_order)
  summary_table <- summary_table[order(summary_table$Group),]
  return(summary_table)
}

summary_table <- summarise_data(gadma_results_pu, param_cols)

# Alternative way because of proportion columns ####
gadma_top <- gadma_results_pu[0,] 
for (group in unique(gadma_results_pu$Groups)) {
  df_group <- gadma_results_pu[gadma_results_pu$Groups == group,]
  
  top_run <- df_group |>
    arrange(desc(LogLikelihood)) |>
    slice_head(n = 3) |>
    slice(c(1, n()))
  
  gadma_top <- rbind(gadma_top, top_run)
}

gadma_top <- gadma_top[,c(1,3,4,param_cols)]

desired_cols <- c("nref", "T1gen", "N1.0", "N2.0", "N1.T1", "N2.T1", "T2gen", "N1.T2", "N2.T2")
existing_cols <- desired_cols[desired_cols %in% names(gadma_top)]
gadma_top[existing_cols] <- lapply(gadma_top[existing_cols], function(x) x / 1000)
custom_order <- c(
  "group1-group2", "group1-group3", "group1-group4",
  "group2-group3", "group2-group4", "group3-group4",
  "group1-Amil", "group2-Amil", "group3-Amil", "group4-Amil"
)

gadma_top$Groups <- factor(gadma_top$Groups, levels = custom_order)
gadma_top <- gadma_top[order(gadma_top$Groups),]
# get rid of rM for now
gadma_top <- gadma_top %>%
  select(-any_of(c("rMe12.T1", "rMe21.T1", "rM12.T1", "rM21.T1",
                   "rMe12.T2", "rMe21.T2", "rM.12.T2", "rM21.T2")))
gadma_top[,2:ncol(gadma_top)] <- round(gadma_top[,2:ncol(gadma_top)], 2)


# Find optimal run
#1
gadma_results[gadma_results$Groups == "group1-group2" & gadma_results$LogLikelihood == -5683.32,]
#2
gadma_results[gadma_results$Groups == "group1-group2" & gadma_results$LogLikelihood == -5993.41,]

gadma_results[gadma_results$Groups == "group1-group3" & gadma_results$LogLikelihood == -6067.67,]

gadma_results[gadma_results$Groups == "group1-group4" & gadma_results$LogLikelihood == -6184.24,]

gadma_results[gadma_results$Groups == "group2-group3" & gadma_results$LogLikelihood == -4913.41,]

gadma_results[gadma_results$Groups == "group2-group4" & gadma_results$LogLikelihood == -4968.81,]

gadma_results[gadma_results$Groups == "group3-group4" & gadma_results$LogLikelihood == -4632.39,]

gadma_results[gadma_results$Groups == "group1-Amil" & gadma_results$LogLikelihood == -7230.64,]

gadma_results[gadma_results$Groups == "group2-Amil" & gadma_results$LogLikelihood == -7114.40,]

gadma_results[gadma_results$Groups == "group3-Amil" & gadma_results$LogLikelihood == -7319.05,]

gadma_results[gadma_results$Groups == "group4-Amil" & gadma_results$LogLikelihood == -6813.63,]

