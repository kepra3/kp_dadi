# Title: Plotting gadma results
# Author: Katharine Prata
# Description: taking in combined and formatted groups and plotting the parameter
# values

library(tidyverse)
# tibble 3.2.1
library(scales)

# Functions
convert_params <- function(gadma_results) {
  skip_groups = character(length = 0L)
  gadma_results$unfil.snps = NA
  gadma_results$sfs.snps = NA
  for (i in 1:nrow(gadma_results)) {
    # effective sequence params
    if (gadma_results[i, 1] == "group1-group2") {
      cat("\nunfiltered snps for: out_QCbasic_noSites_group1_dadi_group2_dadi.pos.gz")
      gadma_results$sfs.snps[i] = 285376.87
      #323336.46
      gadma_results$unfil.snps[i] = 33523695
      } else if (gadma_results[i, 1] == "group1-Amil") {
        cat("\nunfiltered snps for: out_QCbasic_noSites_group1_dadi_Amil_dadi.pos.gz")
        gadma_results$sfs.snps[i] = 111991.19
          #444001.5
        gadma_results$unfil.snps[i]  = 32965107
      } else if (gadma_results[i, 1] == "group1-group3") {
        cat("\nunfiltered snps for: out_QCbasic_noSites_group1_dadi_group3_dadi.pos.gz")
        gadma_results$sfs.snps[i] = 70052.17
        #339902.54
        gadma_results$unfil.snps[i] = 33812349
      } else if (gadma_results[i, 1] == "group1-group4") {
        cat("\nunfiltered snps for: out_QCbasic_noSites_group1_dadi_group4_dadi.pos.gz")
        gadma_results$sfs.snps[i] = 76247.09
        #380214.39
        gadma_results$unfil.snps[i] = 34517876
      } else if (gadma_results[i, 1] == "group2-group3") {
        cat("\nunfiltered snps for: out_QCbasic_noSites_group2_dadi_group3_dadi.pos.gz")
        gadma_results$sfs.snps[i] = 80497.72
        #371442.47
        gadma_results$unfil.snps[i] = 34492471
      } else if (gadma_results[i, 1] == "group2-group4") {
        cat("\nunfiltered snps for: out_QCbasic_noSites_group2_dadi_group4_dadi.pos.gz")
        gadma_results$sfs.snps[i] = 84843.52
          #413972.89
        gadma_results$unfil.snps[i] = 36046531
      } else if (gadma_results[i, 1] == "group3-group4") {
        cat("\nunfiltered snps for: out_QCbasic_noSites_group3_dadi_group4_dadi.pos.gz")
        gadma_results$sfs.snps[i] = 92922.75
          #421890.54
        gadma_results$unfil.snps[i] = 35072420
      } else if (gadma_results[i, 1] == "group2-Amil") {
        cat("\nunfiltered snps for: out_QCbasic_noSites_group2_dadi_Amil_dadi.pos.gz")
        gadma_results$sfs.snps[i] = 118359.4
          #472699.58
        gadma_results$unfil.snps[i] = 34184556
      } else if (gadma_results[i, 1] == "group3-Amil") {
        cat("\nunfiltered snps for: out_QCbasic_noSites_group3_dadi_Amil_dadi.pos.gz")
        gadma_results$sfs.snps[i] = 127781.61
          #500811.35
        gadma_results$unfil.snps[i] = 34598003
      } else if (gadma_results[i, 1] == "group4-Amil") {
        cat("\nunfiltered snps for: out_QCbasic_noSites_group4_dadi_Amil_dadi.pos.gz")
        gadma_results$sfs.snps[i] = 139754.46
        # snps before projection 556220.74
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
  # Converting dadi units to physical units
  # constants
  mu = 1.86*10^(-8)
  gen.small = 2
  gen.big = 5
  L.genom = 390206788
  
  # Step 1. Calculate effective L
  gadma_results$L.effective = L.genom * (gadma_results$sfs.snps/gadma_results$unfil.snps)
  
  # Step 2. Calculate Nref (theta = 4Nref*mu*L)
  gadma_results$nref = gadma_results$theta / (4*mu*gadma_results$L.effective)
  
  # Divergence time
  gadma_results$T1gen = 2 * gadma_results$nref * gadma_results$t1
  gadma_results$T1years.1 = gadma_results$T1gen * gen.small
  gadma_results$T1years.2 = gadma_results$T1gen * gen.big
  
  # Population sizes after divergence
  gadma_results$N1.0 = gadma_results$nu_1 * gadma_results$nref
  gadma_results$N2.0 = gadma_results$nu_2 * gadma_results$nref
  
  # T1 population sizes
  gadma_results$N1.T1 = gadma_results$nu11 * gadma_results$nref
  gadma_results$N2.T1 = gadma_results$nu21 * gadma_results$nref
  
  # T1 migration rate
  gadma_results$rM12.T1 = ifelse(gadma_results$m1_12 == 0, 0,
                                gadma_results$m1_12/(2*gadma_results$nref))
  gadma_results$rM21.T1 = ifelse(gadma_results$m1_21 == 0, 0,
                                gadma_results$m1_21/(2*gadma_results$nref))
  
  # T1 migration rate
  # number of migrants
  gadma_results$M12.T1 = ifelse(gadma_results$rM12.T1 == 0, 0,
                                 gadma_results$rM12.T1 * gadma_results$N1.T1)
  gadma_results$M21.T1 = ifelse(gadma_results$rM21.T1 == 0, 0,
                                 gadma_results$rM21.T1 * gadma_results$N2.T1)
  
  # T2
  gadma_results$T2gen = 2 * gadma_results$nref * gadma_results$t2
  gadma_results$T2years.1 = gadma_results$T2gen * gen.small
  gadma_results$T2years.2 = gadma_results$T2gen * gen.big
  
  # T2 population sizes
  gadma_results$N1.T2 = gadma_results$nu12 * gadma_results$nref
  gadma_results$N2.T2 = gadma_results$nu22 * gadma_results$nref
  
  # T2 migration rate
  gadma_results$rM12.T2 = gadma_results$m2_12/(2*gadma_results$nref)
  gadma_results$rM21.T2 = gadma_results$m2_21/(2*gadma_results$nref)
  
  # T2 migrants
  gadma_results$M12.T2 = ifelse(gadma_results$rM12.T2 == 0, 0,
                                 gadma_results$rM12.T2*gadma_results$N1.T2)
  gadma_results$M21.T2 = ifelse(gadma_results$rM21.T2 == 0, 0,
                                 gadma_results$rM21.T2*gadma_results$N2.T2)
  
  return(gadma_results)
}

setwd("~/git/kp_dadi/scripts/")
# Import results
gadma_results <- read.delim("../results/gadma_results_combined.txt")

# Reformat so each parameter is a column
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

# Fix group comparison names
gadma_results <- gadma_results %>%
  mutate(Groups = str_split_fixed(Directory, "_", 4)[, 3]) %>%
  mutate(Groups = str_replace_all(Groups, "grp", "group")) %>%
  mutate(Groups = str_replace_all(Groups, "-process.*", ""))

gadma_results$Groups <- as.factor(gadma_results$Groups)

for (group in unique(gadma_results$Groups)) {
  gadma_group <- gadma_results[gadma_results$Groups == group,]
  max_idx <- which.max(gadma_group$LogLikelihood)
  print(gadma_group[max_idx,])
}

gadma_results <- gadma_results %>% 
  select(Groups, everything(), -Directory)

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

gadma_results_pu <- convert_params(gadma_results)

# Reshape to long format
gadma_pu_long <- gadma_results_pu %>%
  pivot_longer(
    cols = c(nref:M21.T2),
    names_to = "Physical.Parameter",
    values_to = "Physical.Value"
  )

parameter_order <- c("nref", "T1gen", "T1years.1", "T1years.2", "N1.0", "N2.0",
                     "M12.T1", "M21.T1", "rM12.T1", "rM21.T1", "N1.T1", "N2.T1",
                     "T2gen", "T2years.1", "T2years.2",
                     "M12.T2", "M21.T2", "rM12.T2", "rM21.T2", "N1.T2", "N2.T2")

gadma_pu_long <- gadma_pu_long |>
  mutate(Physical.Parameter = factor(Physical.Parameter, levels = parameter_order))

# Plot
for (group in unique(gadma_pu_long$Groups)) {
  df_group = gadma_pu_long[gadma_pu_long$Groups == group,]
  # Filter to top 5 rows by LogLikelihood (highest = least negative)
  top5_runs <- df_group |>
    dplyr::group_by(Run) |>
    dplyr::summarise(LogLikelihood = mean(LogLikelihood), .groups = "drop") |>
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
  ggsave(paste0("../plots/gadma_params/", group, "_physical-param-confi.pdf"),
         units = "cm", width = 20, height = 20)
}
gadma_top <- gadma_results_pu[0,]
for (group in unique(gadma_results_pu$Groups)) {
  df_group = gadma_results_pu[gadma_results_pu$Groups == group,]
  # Filter to top 5 rows by LogLikelihood (highest = least negative)
  top_run <- df_group |>
    dplyr::arrange(desc(LogLikelihood)) |>
    dplyr::slice_head(n = 1)
  row <- df_group[df_group$LogLikelihood %in% top_run$LogLikelihood, ]
  gadma_top <- rbind(gadma_top, row)
}
div <- gadma_top[order(gadma_top$T1gen), c(1,21)]
div
#div b/w group1-Amil
#220,000 generations
#450,000 - 1,00,000 years

#div b/w group2-Amil
#130,000 generations
#260,000 - 640,000 years

#div b/w group3-Amil
#130,000 generations
#260,000 - 640,000 years

#div b/w group4-Amil
#132,000 generations
#265,000 - 660,000 years

#div b/w group1-group2
#90,000 generations
#170,000 - 450,000 years

#div b/w group1-group3
#200,000 generations
#400,000 - 1,000,000 years

#div b/w group1-group4
#200,000 generations
#400,000 - 1,00,000 years

#div b/w group2-group3
#97,000 generations
#193,000 - 485,000 years

#div b/w group2-group4
#40,000 generations
#170,000 - 410,000 years

#div b/w group3-group4
#84,000 generations
#170,000 - 420,000 years