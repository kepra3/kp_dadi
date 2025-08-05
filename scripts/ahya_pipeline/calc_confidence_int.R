# Title: Calculate and format confidence intervals for GADMA results
# Author: Katharine Prata
# Description: Process confidence interval files and convert to physical units

library(tidyverse)
library(scales)

# CONSTANTS ####
mu = 1.86*10^(-8)
L.genom = 390206788

# SHARED FUNCTIONS ####
# Note: These functions are shared with plot_gadma_results_optim.R
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

# CONSTANTS ####
mu = 1.86*10^(-8)
L.genom = 390206788

# SETTINGS ####
proj <- "proj"
model <- "1het"

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

# CONFIDENCE INTERVALS PROCESSING ####
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

# Convert population sizes to physical units
pop_cols <- c("nu_1_Optimised", "nu_1_Lower_CI", "nu_1_Upper_CI",
              "nu_2_Optimised", "nu_2_Lower_CI", "nu_2_Upper_CI",
              "nu11_Optimised", "nu11_Lower_CI", "nu11_Upper_CI",
              "nu12_Optimised", "nu12_Lower_CI", "nu12_Upper_CI",
              "nu21_Optimised", "nu21_Lower_CI", "nu21_Upper_CI",
              "nu22_Optimised", "nu22_Lower_CI", "nu22_Upper_CI")
for (col in pop_cols) {
  conf_wide[[col]] <- calc_N(conf_wide$nref, conf_wide[[col]])
}

# Convert time parameters to generations
t_cols <- c("t1_Optimised", "t1_Lower_CI", "t1_Upper_CI",
            "t2_Optimised", "t2_Lower_CI", "t2_Upper_CI")
for (col in t_cols) {
  conf_wide[[col]] <- calc_Tgen(conf_wide$nref, conf_wide[[col]])
}

# Convert migration rates
m_cols <- c("me1_12_Optimised", "me1_12_Lower_CI","me1_12_Upper_CI",
            "me1_21_Optimised", "me1_21_Lower_CI","me1_21_Upper_CI",
            "me2_12_Optimised", "me2_12_Lower_CI","me2_12_Upper_CI",
            "me2_21_Optimised", "me2_21_Lower_CI","me2_21_Upper_CI")
for (col in m_cols) {
  conf_wide[[col]] <- calc_rM(conf_wide$nref, conf_wide[[col]])
}

# Convert to number of migrants
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
param_names <- param_names[!param_names %in% other_cols]

# Scale high number parameters (divide by 1000)
high_params <- c("nu_1", "nu_2", "t1", "t2", "nu11", "nu12", "nu21", "nu22")
high_params_all <- unlist(lapply(high_params, function(p) {
  c(paste0(p, "_Optimised"), paste0(p, "_Lower_CI"), paste0(p, "_Upper_CI"))
}))
conf_wide[high_params_all] <- lapply(conf_wide[high_params_all], function(x) x / 1000)

# FORMAT CONFIDENCE INTERVALS ####
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

# Save results
write.csv(conf_wide, "../results/confidence_intervals_physical_units.csv", row.names = FALSE)
write.csv(combined_cols, "../results/confidence_intervals_formatted.csv", row.names = FALSE)

cat("Confidence intervals processed and saved to:\n")
cat("- ../results/confidence_intervals_physical_units.csv\n")
cat("- ../results/confidence_intervals_formatted.csv\n")
