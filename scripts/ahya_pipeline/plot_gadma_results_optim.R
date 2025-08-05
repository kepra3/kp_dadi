# Title: Plotting gadma results
# Author: Katharine Prata
# Description: taking in combined and formatted groups and plotting the parameter
# values

# USAGE: 
# Change the 'model' variable in the MODEL CONFIGURATION section below to:
# - "1het" for 1 heterogeneous gene flow model
# - "2het" for 2 heterogeneous gene flow model  
# - "hetsym" for symmetric heterogeneous gene flow model (single time period)
# - "hetsc" for secondary contact model (two time periods, no T1 gene flow)
# - "" for single migration rate model (two time periods)

library(tidyverse)
# tibble 3.2.1
library(scales)

# Constants ####
# Converting dadi units to physical units constants
mu = 1.86*10^(-8)
#gen.small = 2
#gen.big = 5
L.genom = 390206788

# Functions ####
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
  
  # Population sizes after divergence
  gadma_results$N1.0 <- calc_N(gadma_results$nref, gadma_results$nu_1)
  gadma_results$N2.0 <- calc_N(gadma_results$nref, gadma_results$nu_2)
  
  # T1 population sizes (all models have these)
  gadma_results$N1.T1 <- calc_N(gadma_results$nref, gadma_results$nu11)
  gadma_results$N2.T1 <- calc_N(gadma_results$nref, gadma_results$nu12)
  
  # Divergence time T1 (all models have t1)
  gadma_results$T1gen <- calc_Tgen(gadma_results$nref, gadma_results$t1)
  
  # Model-specific parameters
  if (model == "hetsym") {
    # Single time period model with symmetric heterogeneous gene flow
    # Migration rates for T1 period only
    gadma_results$rMe12.T1 <- calc_rM(gadma_results$nref, gadma_results$me1_12)
    gadma_results$rMe21.T1 <- calc_rM(gadma_results$nref, gadma_results$me1_21)
    # T1 migrants
    gadma_results$Me12.T1 <- calc_M(gadma_results$N1.T1, gadma_results$rMe12.T1)
    gadma_results$Me21.T1 <- calc_M(gadma_results$N2.T1, gadma_results$rMe21.T1)
    
  } else {
    # Two time period models (1het, 2het, hetsc, "") - have both t1 and t2
    # T2 time period
    gadma_results$T2gen <- calc_Tgen(gadma_results$nref, gadma_results$t2)
    gadma_results$T1gen <- gadma_results$T1gen + gadma_results$T2gen
    
    # T2 population sizes
    gadma_results$N1.T2 <- calc_N(gadma_results$nref, gadma_results$nu21)
    gadma_results$N2.T2 <- calc_N(gadma_results$nref, gadma_results$nu22)
    
    if (model == "hetsc") {
      # Secondary contact: NO gene flow in T1, only in T2
      # T2 migration rates only
      gadma_results$rMe12.T2 <- calc_rM(gadma_results$nref, gadma_results$me2_12)
      gadma_results$rMe21.T2 <- calc_rM(gadma_results$nref, gadma_results$me2_21)
      # T2 migrants
      gadma_results$Me12.T2 <- calc_M(gadma_results$N1.T2, gadma_results$rMe12.T2)
      gadma_results$Me21.T2 <- calc_M(gadma_results$N2.T2, gadma_results$rMe21.T2)
      
    } else if (model == "") {
      # Single migration rate model
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
      gadma_results$M21.T2 <- calc_M(gadma_results$N2.T2, gadma_results$rM21.T2)
      
    } else if (model %in% c("2het", "1het")) {
      # Heterogeneous gene flow models (1het and 2het)
      # T1 migration rates
      gadma_results$rMe12.T1 <- calc_rM(gadma_results$nref, gadma_results$me1_12)
      gadma_results$rMe21.T1 <- calc_rM(gadma_results$nref, gadma_results$me1_21)
      # T1 migrants
      gadma_results$Me12.T1 <- calc_M(gadma_results$N1.T1, gadma_results$rMe12.T1)
      gadma_results$Me21.T1 <- calc_M(gadma_results$N2.T1, gadma_results$rMe21.T1)
      # T2 migration rates
      gadma_results$rMe12.T2 <- calc_rM(gadma_results$nref, gadma_results$me2_12)
      gadma_results$rMe21.T2 <- calc_rM(gadma_results$nref, gadma_results$me2_21)
      # T2 migrants
      gadma_results$Me12.T2 <- calc_M(gadma_results$N1.T2, gadma_results$rMe12.T2)
      gadma_results$Me21.T2 <- calc_M(gadma_results$N2.T2, gadma_results$rMe21.T2)
    }
  }
  return(gadma_results)
}
summarise_data <- function(gadma_results_pu, physical_params, model) {
  summary_list <- list()
  for (group in unique(gadma_results_pu$Groups)) {
    df_group <- gadma_results_pu[gadma_results_pu$Groups == group,]
    
    top_runs <- df_group |>
      arrange(desc(LogLikelihood)) |>
      slice_head(n = 3)
    
    # Calculate min and max for physical parameters
    min_vals <- top_runs[, physical_params] |> summarise(across(everything(), min))
    max_vals <- top_runs[, physical_params] |> summarise(across(everything(), max))
    # Divide by 1000 for generations and pop sizes (adjust for different models)
    if (model == "hetsym") {
      # Single time period
      scale_params <- c("nref", "T1gen", "N1.0", "N2.0", "N1.T1", "N2.T1")
    } else {
      # Two time periods (1het, 2het, hetsc, "")
      scale_params <- c("nref", "T1gen", "N1.0", "N2.0", "N1.T1", "N2.T1", "T2gen", "N1.T2", "N2.T2")
    }
    existing_cols <- scale_params[scale_params %in% names(min_vals)]
    min_vals[existing_cols] <- lapply(min_vals[existing_cols], function(x) x / 1000)
    max_vals[existing_cols] <- lapply(max_vals[existing_cols], function(x) x / 1000)
    
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

setwd("~/git/kp_dadi/scripts/")

# MAIN ANALYSIS ####
# MODEL CONFIGURATION ####
# Set model type: "1het", "2het", "hetsym", "hetsc", or ""
model <- "hetsc"  # Change this to switch between models
proj <- "proj"
optim <- paste(proj, model, sep = "_")

# Load appropriate results file based on model
if (model == "1het") {
  results_file <- "../results/gadma_proj_1het_results_combined_noerror.txt"
} else if (model == "2het") {
  results_file <- "../results/gadma_proj_2het_results_combined_noerror.txt"
} else if (model == "hetsym") {
  results_file <- "../results/gadma_proj_hetsym_results_combined_noerror.txt"
} else if (model == "hetsc") {
  results_file <- "../results/gadma_proj_hetsc_results_combined_noerror.txt"
} else if (model == "") {
  results_file <- "../results/gadma_proj_results_combined.txt"
} else {
  stop("Unknown model type. Please use: '1het', '2het', 'hetsym', 'hetsc', or ''")
}

# Import results ####
cat("Loading results for model:", model, "\n")
cat("File:", results_file, "\n")
gadma_results <- read.delim(results_file)

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

# Clean up column names and data
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

# Plotting ####
# plots
for (group in levels(gadma_long$Groups)) {
  df_group <- gadma_long[gadma_long$Groups == group,]
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

# Find best models by group
for (group in unique(gadma_results$Groups)) {
  gadma_group <- gadma_results[gadma_results$Groups == group,]
  max_idx <- which.max(gadma_group$LogLikelihood)
  cat("Best model for", group, ":\n")
  print(gadma_group[max_idx,])
  cat("\n")
}

# Physical unit conversions ####
gadma_results_pu <- convert_params(gadma_results,
                                   model = model,
                                   proj = proj,
                                   L.genom = L.genom,
                                   mu = mu)

# Parameter selection ####
# Define explicit parameter columns for physical unit plotting based on model
if (model == "2het") {
  physical_params <- c("nref", "T1gen", "N1.0", "N2.0", "N1.T1", "N2.T1", 
                       "Me12.T1", "Me21.T1", "M12.T1", "M21.T1", "T2gen", 
                       "N1.T2", "N2.T2", "Me12.T2", "Me21.T2", "M12.T2", "M21.T2")
} else if (model == "1het") {
  physical_params <- c("nref", "T1gen", "N1.0", "N2.0", "N1.T1", "N2.T1", 
                       "Me12.T1", "Me21.T1", "T2gen", "N1.T2", "N2.T2", 
                       "Me12.T2", "Me21.T2", "P1", "P2")
} else if (model == "hetsym") {
  # Single time period with symmetric heterogeneous gene flow
  physical_params <- c("nref", "T1gen", "N1.0", "N2.0", "N1.T1", "N2.T1", 
                       "Me12.T1", "Me21.T1", "P1")
} else if (model == "hetsc") {
  # Secondary contact: two time periods, no T1 gene flow, only T2 gene flow
  physical_params <- c("nref", "T1gen", "N1.0", "N2.0", "N1.T1", "N2.T1", 
                       "T2gen", "N1.T2", "N2.T2", "Me12.T2", "Me21.T2", "P2")
} else if (model == "") {
  # Single migration rate model (two time periods)
  physical_params <- c("nref", "T1gen", "N1.0", "N2.0", "N1.T1", "N2.T1", 
                       "M12.T1", "M21.T1", "T2gen", "N1.T2", "N2.T2", 
                       "M12.T2", "M21.T2")
} else {
  stop("Unknown model type for parameter selection")
}

cat("Physical parameters for model", model, ":\n")
cat(paste(physical_params, collapse = ", "), "\n\n")

# Reshape to long format
gadma_pu_long <- gadma_results_pu %>%
  pivot_longer(
    cols = all_of(physical_params),
    names_to = "Physical.Parameter",
    values_to = "Physical.Value"
  )

gadma_pu_long <- gadma_pu_long |>
  mutate(Physical.Parameter = factor(Physical.Parameter, levels = physical_params))

# Plot physical parameters
for (group in unique(gadma_pu_long$Groups)) {
  df_group <- gadma_pu_long[gadma_pu_long$Groups == group,]
  # Filter to top 5 runs by LogLikelihood (highest = least negative)
  top5_runs <- df_group |>
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
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_text(margin = margin(t = 10))
    ) +
    ylab("Parameter Value") +
    xlab("Parameter") +
    labs(
        title = group,
        subtitle = "Confidence intervals"
      )
  print(p)
  ggsave(paste0("../plots/gadma_params/", group, "_", optim,
                "_physical-param-confi.pdf"),
         units = "cm", width = 20, height = 20)
}

# SUMMARY TABLES ####
summary_table <- summarise_data(gadma_results_pu, physical_params, model)
# Alternative approach for summary tables ####
gadma_top <- gadma_results_pu[0,] 
for (group in unique(gadma_results_pu$Groups)) {
  df_group <- gadma_results_pu[gadma_results_pu$Groups == group,]
  
  top_run <- df_group |>
    arrange(desc(LogLikelihood)) |>
    slice_head(n = 3) |>
    slice(c(1, n()))
  
  gadma_top <- rbind(gadma_top, top_run)
}

gadma_top <- gadma_top[,c("Groups", "LogLikelihood", physical_params)]

# Scale parameters based on model type
if (model == "hetsym") {
  # Single time period
  scale_params <- c("nref", "T1gen", "N1.0", "N2.0", "N1.T1", "N2.T1")
} else {
  # Two time periods (1het, 2het, hetsc, "")
  scale_params <- c("nref", "T1gen", "N1.0", "N2.0", "N1.T1", "N2.T1", "T2gen", "N1.T2", "N2.T2")
}
existing_cols <- scale_params[scale_params %in% names(gadma_top)]
gadma_top[existing_cols] <- lapply(gadma_top[existing_cols], function(x) x / 1000)
custom_order <- c(
  "group1-group2", "group1-group3", "group1-group4",
  "group2-group3", "group2-group4", "group3-group4",
  "group1-Amil", "group2-Amil", "group3-Amil", "group4-Amil"
)

gadma_top$Groups <- factor(gadma_top$Groups, levels = custom_order)
gadma_top <- gadma_top[order(gadma_top$Groups),]

gadma_top[,2:ncol(gadma_top)] <- round(gadma_top[,2:ncol(gadma_top)], 2)

# EXPORT RESULTS ####
# 1. Export top models with original dadi parameters
cat("\n=== SAVING RESULTS ===\n")

# Get best model for each group with original dadi parameters
best_models_dadi <- data.frame()
for (group in unique(gadma_results$Groups)) {
  gadma_group <- gadma_results[gadma_results$Groups == group,]
  best_idx <- which.max(gadma_group$LogLikelihood)
  best_model <- gadma_group[best_idx,]
  best_models_dadi <- rbind(best_models_dadi, best_model)
}

# Order by group
best_models_dadi$Groups <- factor(best_models_dadi$Groups, levels = custom_order)
best_models_dadi <- best_models_dadi[order(best_models_dadi$Groups),]

# Write dadi parameters
output_file_dadi <- paste0("../results/", optim, "_best_models_dadi_params.csv")
write.csv(best_models_dadi, output_file_dadi, row.names = FALSE)
cat("Top models with dadi parameters saved to:", output_file_dadi, "\n")

# 2. Export summary with physical parameters
output_file_physical <- paste0("../results/", optim, "_summary_physical_params.csv")
write.csv(gadma_top, output_file_physical, row.names = FALSE)
cat("Summary with physical parameters saved to:", output_file_physical, "\n")

# 3. Export summary statistics table
output_file_summary <- paste0("../results/", optim, "_parameter_summary_table.csv")
write.csv(summary_table, output_file_summary, row.names = FALSE)
cat("Parameter summary table saved to:", output_file_summary, "\n")

# Display quick summary
cat("\n=== QUICK SUMMARY ===\n")
cat("Model type:", model, "\n")
cat("Projection:", proj, "\n")
cat("Number of group comparisons:", nrow(best_models_dadi), "\n")
cat("Groups analysed:", paste(levels(best_models_dadi$Groups), collapse = ", "), "\n")
cat("\nFiles saved:\n")
cat("- Dadi parameters:", output_file_dadi, "\n")
cat("- Physical parameters:", output_file_physical, "\n") 
cat("- Summary table:", output_file_summary, "\n")


