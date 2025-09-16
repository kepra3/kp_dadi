# Title: Calculate and format confidence intervals for GADMA results
# Author: Katharine Prata
# Description: Process confidence interval files and convert to physical units
#
# USAGE: 
# Change the 'model' variable in the MODEL CONFIGURATION section below to:
# - "1het" for 1 heterogeneous gene flow model
# - "2het" for 2 heterogeneous gene flow model  
# - "hetsym" for continuous heterogeneous gene flow model (single time period)
# - "hetsc" for secondary contact model (two time periods, no T1 gene flow)
# - "" for single migration rate model (two time periods)

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

# SETTINGS ####
proj <- "proj"
# MODEL CONFIGURATION ####
# Set model type: "1het", "2het", "hetsym", "hetsc", or ""
model <- "hetsym"  # Change this to switch between models

# Load confidence interval files ####
conf_groups <- c(
  "group1-group2"#, "gadma/group1-group3", "gadma/group1-group4",
  #"group2-group3", "group2-group4", "group3-group4",
  #"gadma/group1-Amil", "gadma/group2-Amil", "gadma/group3-Amil", "gadma/group4-Amil"
)
conf_files <- paste0("../results/", conf_groups, "_projected0.8_", model, "_confidence_intervals.txt")

# Filter to only existing files
existing_files <- file.exists(conf_files)
conf_groups <- conf_groups[existing_files]
conf_files <- conf_files[existing_files]

if (length(conf_files) == 0) {
  stop("No confidence interval files found for model: ", model)
}

cat("Loading confidence interval files for model:", model, "\n")
cat("Files found:", length(conf_files), "\n")
for (i in seq_along(conf_files)) {
  cat("  -", conf_files[i], "\n")
}

all_conf_list <- lapply(seq_along(conf_files), function(i) {
  df <- read.delim(conf_files[i])
  df$Group <- conf_groups[i]
  df
})
all_conf <- do.call(rbind, all_conf_list)
all_conf$Optimised_thres_max <- all_conf$Optimised/10000
all_conf$Optimised_thres_min <- all_conf$Optimised/1000
all_conf$chosen_eps <- (all_conf$Optimised < 100 & all_conf$eps > all_conf$Optimised_thres_max & all_conf$eps < all_conf$Optimised_thres_min) |
      (all_conf$Optimised >= 100 & all_conf$eps == 1)

# Remove gadma prefix
all_conf$Group <- str_remove(all_conf$Group, "^gadma/")

# MODEL-SPECIFIC EPSILON ADJUSTMENTS ####
# For some groups/parameters, low eps settings increase the upper confidence interval
# This is probably because the lower bound is unable to be estimated
# Apply model-specific adjustments for problematic parameter/group combinations

if (model == "1het") {
  # For 1het model: adjust t2 parameter for specific groups
  groups_t2_adjust <- c("group2-group3", "group2-group4", "group3-group4")
  if (any(all_conf$Group %in% groups_t2_adjust & all_conf$Parameter == "t2")) {
    cat("Applying t2 epsilon adjustments for 1het model\n")
    all_conf <- all_conf %>%
      group_by(Group, Parameter) %>%
      mutate(chosen_eps = ifelse(Parameter == "t2" & Group %in% groups_t2_adjust, FALSE, chosen_eps),
             chosen_eps = ifelse(Parameter == "t2" & Group %in% groups_t2_adjust & Upper_CI == min(Upper_CI), TRUE, chosen_eps)) %>%
      ungroup()
    groups_t2_adjust2 <- c("group3-Amil")
    all_conf <- all_conf %>%
      group_by(Group, Parameter) %>%
      mutate(chosen_eps = ifelse(Parameter == "t2" & Group %in% groups_t2_adjust2, FALSE, chosen_eps),
             chosen_eps = ifelse(Parameter == "t2" & Group %in% groups_t2_adjust2 & eps == 1e-04, TRUE, chosen_eps)) %>%
      ungroup()
  }
  all_conf$Upper_CI[(all_conf$Parameter == "P1" | all_conf$Parameter == "P2") & all_conf$Upper_CI > 1] <- 1
} else if (model == "hetsc") {
  # For hetsc model: adjust specific parameters for group1-group2
  if (any(all_conf$Group == "group1-group2")) {
    cat("Applying parameter adjustments for hetsc model\n")
    # Adjust t2 parameter
    if (any(all_conf$Group == "group1-group2" & all_conf$Parameter == "t2")) {
      all_conf <- all_conf %>%
        group_by(Group, Parameter) %>%
        mutate(chosen_eps = ifelse(Parameter == "t2" & Group == "group1-group2", FALSE, chosen_eps),
               chosen_eps = ifelse(Parameter == "t2" & Group == "group1-group2" & eps == 1e-03, TRUE, chosen_eps)) %>%
        ungroup()
    }
    # Adjust P2 parameter
    if (any(all_conf$Group == "group1-group2" & all_conf$Parameter == "P2")) {
      all_conf <- all_conf %>%
        group_by(Group, Parameter) %>%
        mutate(chosen_eps = ifelse(Parameter == "P2" & Group == "group1-group2", FALSE, chosen_eps),
               chosen_eps = ifelse(Parameter == "P2" & Group == "group1-group2" & eps == 1, TRUE, chosen_eps)) %>%
        ungroup()
    }
  }
  
} else if (model == "2het") {
  # Add 2het-specific adjustments if needed
  cat("No specific epsilon adjustments defined for 2het model\n")
  
} else if (model == "hetsym") {
  # Add hetsym-specific adjustments if needed
  cat("No specific epsilon adjustments defined for hetsym model\n")
  
} else if (model == "") {
  # Add single migration rate model adjustments if needed
  cat("No specific epsilon adjustments defined for single migration rate model\n")
}


all_conf <- all_conf %>%
  group_by(Group, Parameter) %>%
  filter(if (any(chosen_eps)) chosen_eps else eps == max(eps)) %>%
  ungroup()

# MODEL-SPECIFIC PARAMETER ORDERING ####
# Define parameter order based on model type
if (model == "1het") {
  order_params <- c("t1", "nu_1", "nu_2", "me1_12", "me1_21", "P1", "nu11", "nu12", "t2", "me2_12", "me2_21", "P2", "nu21", "nu22", "theta")
} else if (model == "2het") {
  order_params <- c("t1", "nu_1", "nu_2", "me1_12", "me1_21", "M1_12", "M1_21", "nu11", "nu12", "t2", "nu21", "nu22", "me2_12", "me2_21", "M2_12", "M2_21", "theta")
} else if (model == "hetsym") {
  # Single time period with symmetric heterogeneous gene flow
  order_params <- c("t1", "nu_1", "nu_2", "me1_12", "me1_21", "P1", "nu11", "nu12", "theta")
} else if (model == "hetsc") {
  # Secondary contact: two time periods, no T1 gene flow, only T2 gene flow
  order_params <- c("t1", "nu_1", "nu_2", "nu11", "nu12", "t2", "me2_12", "me2_21", "P2", "nu21", "nu22", "theta")
} else if (model == "") {
  # Single migration rate model (two time periods)
  order_params <- c("t1", "nu_1", "nu_2", "m1_12", "m1_21", "nu11", "nu12", "t2", "m2_12", "m2_21", "nu21", "nu22", "theta")
} else {
  stop("Unknown model type: ", model, ". Please use: '1het', '2het', 'hetsym', 'hetsc', or ''")
}
# Apply parameter ordering (only for parameters that exist in the data)
existing_params <- intersect(order_params, unique(all_conf$Parameter))
all_conf$Parameter <- factor(all_conf$Parameter, levels = existing_params, ordered = TRUE)
all_conf <- all_conf[order(all_conf$Parameter), ]

cat("Parameters found for model", model, ":\n")
cat(paste(existing_params, collapse = ", "), "\n\n")

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

cat("Columns available in conf_wide:\n")
cat(paste(names(conf_wide), collapse = ", "), "\n\n")

# Convert population sizes to physical units
pop_cols <- c("nu_1_Optimised", "nu_1_Lower_CI", "nu_1_Upper_CI",
              "nu_2_Optimised", "nu_2_Lower_CI", "nu_2_Upper_CI",
              "nu11_Optimised", "nu11_Lower_CI", "nu11_Upper_CI",
              "nu12_Optimised", "nu12_Lower_CI", "nu12_Upper_CI",
              "nu21_Optimised", "nu21_Lower_CI", "nu21_Upper_CI",
              "nu22_Optimised", "nu22_Lower_CI", "nu22_Upper_CI")

# Only process existing columns
existing_pop_cols <- pop_cols[pop_cols %in% names(conf_wide)]
for (col in existing_pop_cols) {
  conf_wide[[col]] <- calc_N(conf_wide$nref, conf_wide[[col]])
}

# Convert time parameters to generations
t_cols <- c("t1_Optimised", "t1_Lower_CI", "t1_Upper_CI",
            "t2_Optimised", "t2_Lower_CI", "t2_Upper_CI")

# Only process existing columns
existing_t_cols <- t_cols[t_cols %in% names(conf_wide)]
for (col in existing_t_cols) {
  conf_wide[[col]] <- calc_Tgen(conf_wide$nref, conf_wide[[col]])
}

# Convert migration rates (model-specific)
if (model %in% c("1het", "2het", "hetsym", "hetsc")) {
  # Heterogeneous migration models use "me" parameters
  m_cols <- c("me1_12_Optimised", "me1_12_Lower_CI", "me1_12_Upper_CI",
              "me1_21_Optimised", "me1_21_Lower_CI", "me1_21_Upper_CI",
              "me2_12_Optimised", "me2_12_Lower_CI", "me2_12_Upper_CI",
              "me2_21_Optimised", "me2_21_Lower_CI", "me2_21_Upper_CI")
} else if (model == "") {
  # Single migration rate model uses "m" parameters
  m_cols <- c("m1_12_Optimised", "m1_12_Lower_CI", "m1_12_Upper_CI",
              "m1_21_Optimised", "m1_21_Lower_CI", "m1_21_Upper_CI",
              "m2_12_Optimised", "m2_12_Lower_CI", "m2_12_Upper_CI",
              "m2_21_Optimised", "m2_21_Lower_CI", "m2_21_Upper_CI")
} else {
  m_cols <- character(0)  # No migration columns
}

# Only process existing columns
existing_m_cols <- m_cols[m_cols %in% names(conf_wide)]
for (col in existing_m_cols) {
  conf_wide[[col]] <- calc_rM(conf_wide$nref, conf_wide[[col]])
}

# Convert to number of migrants (model-specific)
if (model %in% c("1het", "2het", "hetsym", "hetsc")) {
  # Heterogeneous migration models
  # Only process if the required columns exist
  
  # me1_12 migration (T1 period: pop1 -> pop2)
  me1_12_cols <- c("me1_12_Optimised", "me1_12_Lower_CI", "me1_12_Upper_CI")
  if (all(me1_12_cols %in% names(conf_wide)) && "nu11_Optimised" %in% names(conf_wide)) {
    for (col in me1_12_cols) {
      conf_wide[[col]] <- calc_M(conf_wide$nu11_Optimised, conf_wide[[col]])
    }
  }
  
  # me1_21 migration (T1 period: pop2 -> pop1)
  me1_21_cols <- c("me1_21_Optimised", "me1_21_Lower_CI", "me1_21_Upper_CI")
  if (all(me1_21_cols %in% names(conf_wide)) && "nu12_Optimised" %in% names(conf_wide)) {
    for (col in me1_21_cols) {
      conf_wide[[col]] <- calc_M(conf_wide$nu12_Optimised, conf_wide[[col]])
    }
  }
  
  # me2_12 migration (T2 period: pop1 -> pop2) - only for models with T2
  if (model %in% c("1het", "2het", "hetsc")) {
    me2_12_cols <- c("me2_12_Optimised", "me2_12_Lower_CI", "me2_12_Upper_CI")
    if (all(me2_12_cols %in% names(conf_wide)) && "nu21_Optimised" %in% names(conf_wide)) {
      for (col in me2_12_cols) {
        conf_wide[[col]] <- calc_M(conf_wide$nu21_Optimised, conf_wide[[col]])
      }
    }
    
    # me2_21 migration (T2 period: pop2 -> pop1)
    me2_21_cols <- c("me2_21_Optimised", "me2_21_Lower_CI", "me2_21_Upper_CI")
    if (all(me2_21_cols %in% names(conf_wide)) && "nu22_Optimised" %in% names(conf_wide)) {
      for (col in me2_21_cols) {
        conf_wide[[col]] <- calc_M(conf_wide$nu22_Optimised, conf_wide[[col]])
      }
    }
  }
  
} else if (model == "") {
  # Single migration rate model - similar structure but with "m" instead of "me"
  # m1_12 migration (T1 period)
  m1_12_cols <- c("m1_12_Optimised", "m1_12_Lower_CI", "m1_12_Upper_CI")
  if (all(m1_12_cols %in% names(conf_wide)) && "nu11_Optimised" %in% names(conf_wide)) {
    for (col in m1_12_cols) {
      conf_wide[[col]] <- calc_M(conf_wide$nu11_Optimised, conf_wide[[col]])
    }
  }
  
  # m1_21 migration (T1 period)
  m1_21_cols <- c("m1_21_Optimised", "m1_21_Lower_CI", "m1_21_Upper_CI")
  if (all(m1_21_cols %in% names(conf_wide)) && "nu12_Optimised" %in% names(conf_wide)) {
    for (col in m1_21_cols) {
      conf_wide[[col]] <- calc_M(conf_wide$nu12_Optimised, conf_wide[[col]])
    }
  }
  
  # m2_12 migration (T2 period)
  m2_12_cols <- c("m2_12_Optimised", "m2_12_Lower_CI", "m2_12_Upper_CI")
  if (all(m2_12_cols %in% names(conf_wide)) && "nu21_Optimised" %in% names(conf_wide)) {
    for (col in m2_12_cols) {
      conf_wide[[col]] <- calc_M(conf_wide$nu21_Optimised, conf_wide[[col]])
    }
  }
  
  # m2_21 migration (T2 period)
  m2_21_cols <- c("m2_21_Optimised", "m2_21_Lower_CI", "m2_21_Upper_CI")
  if (all(m2_21_cols %in% names(conf_wide)) && "nu22_Optimised" %in% names(conf_wide)) {
    for (col in m2_21_cols) {
      conf_wide[[col]] <- calc_M(conf_wide$nu22_Optimised, conf_wide[[col]])
    }
  }
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

# Only process existing columns
existing_high_params <- high_params_all[high_params_all %in% names(conf_wide)]
conf_wide[existing_high_params] <- lapply(conf_wide[existing_high_params], function(x) x / 1000)

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

# Save results with model-specific filenames
output_suffix <- ifelse(model == "", "single_migration", model)
physical_units_file <- paste0("../results/confidence_intervals_", output_suffix, "_physical_units.csv")
formatted_file <- paste0("../results/confidence_intervals_", output_suffix, "_formatted.csv")

write.csv(conf_wide, physical_units_file, row.names = FALSE)
write.csv(combined_cols, formatted_file, row.names = FALSE)

cat("Confidence intervals processed and saved to:\n")
cat("- ", physical_units_file, "\n")
cat("- ", formatted_file, "\n")

