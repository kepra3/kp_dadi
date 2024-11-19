# Title: summarise dadi results
# Author: Katharine Prata
# Date created: 18/03/24

# Packages
library(ggplot2)
library(tidyr)
library(dplyr)

# Functions
plot.results <- function(results) {
  for (pop in levels(results$Pop)) {
    print(pop)
    group <- results[results$Pop == pop,]
    p <- ggplot(group, aes(Model, AIC)) + geom_point() + theme_bw() + ggtitle(pop)
    print(p)
    ggsave(paste0("../plots/model_summaries/", pop, ".model.summary.pdf"), p, height = 5, width = 5,
           units = "cm", dpi = 400)
    p1 <- ggplot(group[group$AIC < quantile(group$AIC, 0.9),], aes(Model, AIC)) + 
      geom_point() + theme_bw() + ggtitle(paste0(pop, " - Top 90%")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    print(p1)
    ggsave(paste0("../plots/model_summaries/", pop, ".AIC.top90.summary.pdf"), p1, height = 15, width = 15,
           units = "cm", dpi = 400)
    p2 <- ggplot(group[group$log.likelihood > quantile(group$log.likelihood, 0.9),],
                 aes(Model, log.likelihood)) + 
      geom_point() + theme_bw() + ggtitle(paste0(pop, " - Top 90%")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    print(p2)
    ggsave(paste0("../plots/model_summaries/", pop, ".log.likelihood.top90.summary.pdf"), p1, height = 15, width = 15,
           units = "cm", dpi = 400)
    print(as.data.frame(group[group$AIC == min(group$AIC),])[1,])
    #for (model in levels(results$Model)) {
    #  print(as.data.frame(group[group$AIC == min(group$AIC[group$Model == model]),])[1,])
    #}
  }
}

setwd("~/git/kp_dadi/scripts/")

# Args

# Load data
results <- read.table("../results/dadi_optimisation_combined.txt", sep = "\t",
                      header = TRUE)
num_runs <- 3

#results <- results[331:length(results[,1]),]

results$Pop <- as.factor(results$Pop)
results$Model <- as.factor(results$Model)

table(results$Pop, results$Model)

str(results)
results$log.likelihood <- as.numeric(results$log.likelihood)
str(results)

results <- results[results$log.likelihood < -15,]
results[results$log.likelihood == -1297.27,] <- NA
results <- na.omit(results)

# Remove complicated models (including het and test between no migration vs. migration models)
#remove_models <- c("het_asym_mig", "anc_het_asym_mig", "sec_het_asym_mig",
#                   "split_bottle_het_asym_mig", "split_bottle_anc_het_asym_mig",
#                   "split_bottle_sec_het_asym_mig", "split_sizechange_het_asym_mig",
#                   "split_sizechange_anc_het_asym_mig",
#                   "split_sizechange_sec_het_asym_mig")

#results2 <- results[!results$Model %in% remove_models,]

plot.results(results)
#plot.results(results2)
#results <- results2

dat <- data.frame(Pop = numeric(0), Model = numeric(0), Opt = numeric(0))
for (group in sort(unique(results$Pop))) {
    subset <- results[results$Pop == group,]
    row <- head(subset[order(subset$AIC),],num_runs)
    row <- row[,c(1,3,5,6,10,11)]
    dat <- rbind(dat, row)
  }
dat2 <- data.frame(Pop = numeric(0), Model = numeric(0), Opt = numeric(0))
for (group in sort(unique(results$Pop))) {
  for (model in sort(unique(results$Model))) {
  subset <- results[results$Pop == group & results$Model == model,]
  row <- head(subset[order(subset$AIC),],num_runs)
  row <- row[,c(1,3,5,6,10,11)]
  dat2 <- rbind(dat2, row)
  }}

dat_nomig_mig <- data.frame(Pop = numeric(0), Model = numeric(0), Opt = numeric(0))
no.mig_models <- c("no_mig", "split_bottle", "split_sizechange")
mig_models <- levels(dat2$Model)
mig_models <- mig_models[!mig_models %in% no.mig_models]
for (group in sort(unique(results$Pop))) {
  subset <- results[results$Pop == group & results$Model %in% no.mig_models,]
  row <- head(subset[order(subset$AIC),],num_runs)
  row <- row[,c(1,3,5,6,10,11)]
  dat_nomig_mig <- rbind(dat_nomig_mig, row)
  subset <- results[results$Pop == group & results$Model %in% mig_models,]
  row <- head(subset[order(subset$AIC),],num_runs)
  row <- row[,c(1,3,5,6,10,11)]
  dat_nomig_mig <- rbind(dat_nomig_mig, row)
  }

# Shorten pop names
dat <- dat %>% separate(Pop, into = c("Pop1","2","Pop2"), sep = "_", remove = FALSE) %>% 
  separate(Pop1, into = c( "x", "y", "Pop1"), sep = "\\.") %>% 
  unite(Pop_short, c("Pop1", "Pop2"), sep = "-")
dat <- dat[,c(-2,-3,-5)]
dat2 <- dat2 %>% separate(Pop, into = c("Pop1","2","Pop2"), sep = "_", remove = FALSE) %>% 
  separate(Pop1, into = c( "x", "y", "Pop1"), sep = "\\.") %>% 
  unite(Pop_short, c("Pop1", "Pop2"), sep = "-")
dat2 <- dat2[,c(-2,-3,-5)]
dat_nomig_mig <- dat_nomig_mig %>% separate(Pop, into = c("Pop1","2","Pop2"), sep = "_", remove = FALSE) %>% 
  separate(Pop1, into = c( "x", "y", "Pop1"), sep = "\\.") %>% 
  unite(Pop_short, c("Pop1", "Pop2"), sep = "-")
dat_nomig_mig <- dat_nomig_mig[,c(-2,-3,-5)]


pop_order <- c("group1-group2", "group1-group3", "group1-group4",
               "group2-group3", "group2-group4", "group3-group4",
               "group1-Amil", "group2-Amil", "group3-Amil", "group4-Amil")

dat$Pop_short <- factor(dat$Pop_short, levels = pop_order)
dat2$Pop_short <- factor(dat2$Pop_short, levels = pop_order)
dat_nomig_mig$Pop_short <- factor(dat_nomig_mig$Pop_short, levels = pop_order)

# Best models ####

ggplot(dat, aes(AIC, Model, color = AIC)) + 
  geom_point() + 
  scale_color_gradient(low = "green", high = "red") +
  theme_bw() +
  facet_wrap(~Pop_short, nrow=4, ncol=3) +
  theme(axis.text.x = element_text(angle = 90))
ggsave(paste0("../plots/model_summaries/Ahya_", num_runs, "_13-11-24.pdf"), units = "cm",
       width = 25, height = 30)

ggplot(dat[dat$Pop_short %in% pop_order[1:6],], aes(AIC, Model, color = AIC)) + 
  geom_point() + 
  scale_color_gradient(low = "green", high = "red") +
  theme_bw() +
  facet_wrap(~Pop_short, nrow=2, ncol=3) +
  theme(axis.text.x = element_text(angle = 90))
ggsave(paste0("../plots/model_summaries/Ahya-within_", num_runs, "_13-11-24.pdf"), units = "cm",
       width = 30, height = 25)

ggplot(dat[dat$Pop_short %in% pop_order[7:10],], aes(AIC, Model, color = AIC)) + 
  geom_point() + 
  scale_color_gradient(low = "green", high = "red") +
  theme_bw() +
  facet_wrap(~Pop_short, nrow=1, ncol=4) +
  theme(axis.text.x = element_text(angle = 90))
ggsave(paste0("../plots/model_summaries/Ahya-Amil_", num_runs, "_13-11-24.pdf"), units = "cm",
       width = 30, height = 15)

# Per group pair
for (pop in levels(dat$Pop_short)) {
  ggplot(dat[dat$Pop_short == pop,], aes(AIC, Model, color = AIC)) + 
    geom_point() + 
    scale_color_gradient(low = "green", high = "red") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90))
  ggsave(paste0("../plots/model_summaries/",
                pop, "_", num_runs, "_11-11-24.pdf"),
         units = "cm", width = 15, height = 15)
}

# All models ####

ggplot(dat2[dat2$Pop_short %in% pop_order[1:6],], aes(AIC, Model, color = AIC)) + 
  geom_point() + 
  scale_color_gradient(low = "green", high = "red") +
  theme_bw() +
  facet_wrap(~Pop_short, nrow=2, ncol=3) +
  theme(axis.text.x = element_text(angle = 90))
ggsave(paste0("../plots/model_summaries/Ahya-within_all-models-", num_runs, "_13-11-24.pdf"), units = "cm",
       width = 30, height = 25)

ggplot(dat2[dat2$Pop_short %in% pop_order[7:10],], aes(AIC, Model, color = AIC)) + 
  geom_point() + 
  scale_color_gradient(low = "green", high = "red") +
  theme_bw() +
  facet_wrap(~Pop_short, nrow=1, ncol=4) +
  theme(axis.text.x = element_text(angle = 90))
ggsave(paste0("../plots/model_summaries/Ahya-Amil_all-models-", num_runs, "_13-11-24.pdf"), units = "cm",
       width = 30, height = 15)

# No migration vs. migration ####
ggplot(dat_nomig_mig[dat_nomig_mig$Pop_short %in% pop_order[1:6],], aes(AIC, Model, color = AIC)) + 
  geom_point() + 
  scale_color_gradient(low = "green", high = "red") +
  theme_bw() +
  facet_wrap(~Pop_short, nrow=2, ncol=3) +
  theme(axis.text.x = element_text(angle = 90))
ggsave(paste0("../plots/model_summaries/Ahya-within_nomig-mig-", num_runs, "_13-11-24.pdf"), units = "cm",
       width = 30, height = 25)

ggplot(dat_nomig_mig[dat_nomig_mig$Pop_short %in% pop_order[7:10],], aes(AIC, Model, color = AIC)) + 
  geom_point() + 
  scale_color_gradient(low = "green", high = "red") +
  theme_bw() +
  facet_wrap(~Pop_short, nrow=1, ncol=4) +
  theme(axis.text.x = element_text(angle = 90))
ggsave(paste0("../plots/model_summaries/Ahya-Amil_nomig-mig-", num_runs, "_13-11-24.pdf"), units = "cm",
       width = 30, height = 15)

# Per group pair
for (pop in levels(dat$Pop_short)) {
  ggplot(dat_nomig_mig[dat_nomig_mig$Pop_short == pop,], aes(AIC, Model, color = AIC)) + 
    geom_point() + 
    scale_color_gradient(low = "green", high = "red") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90))
  ggsave(paste0("../plots/model_summaries/",
                pop, "_nomig-mig_", num_runs, "_11-11-24.pdf"),
         units = "cm", width = 15, height = 15)
}


# Calculate the number of comma-separated values in each row
comma_counts <- sapply(results$optimised_params_labels, function(x) length(strsplit(x, ",")[[1]]))

# Find the maximum count and the row(s) where it occurs
max_count <- max(comma_counts)
max_rows <- which(comma_counts == max_count)

# Swap param labels and params
dat <- dat[,c(1:5,7,6)]

dat_params <- separate(dat, optimised_params,
                       into=unlist(strsplit(gsub("\\s+", "", names(max_rows[1])), ",")),
                       sep=",", remove = FALSE)
dat_reshape <- dat_params
for (i in 1:length(dat_params[,1])) {
  vars <- dat_params[i,8:22]
  if (dat_params$Model[i] == "asym_mig" | 
            dat_params$Model[i] == "anc_asym_mig" | 
            dat_params$Model[i] == "sec_cont_asym_mig") {
   dat_reshape[i,c(10:11)] <- NA
   dat_reshape[i,c(14:15)] <- vars[3:4]
 } else if (dat_params$Model[i] == "het_asym_mig") {
   dat_reshape[i,c(10:11,13,18)] <- NA
   dat_reshape[i,12] <- vars[7]
   dat_reshape[i, 14:15] <- vars[3:4]
   dat_reshape[i, 16:17] <- vars[5:6]
   dat_reshape[i, 22] <- vars[8]
 } else if (dat_params$Model[i] == "sec_het_asym_mig" |
            dat_params$Model[i] == "anc_het_asym_mig") {
   dat_reshape[i,c(10:11,13,18)] <- NA
   dat_reshape[i,12:13] <- vars[7:8]
   dat_reshape[i, 14:15] <- vars[3:4]
   dat_reshape[i, 16:17] <- vars[5:6]
   dat_reshape[i, 22] <- vars[9]
 } else if(dat_params$Model[i] == "split_bottle_anc_het_asym_mig" | 
           dat_params$Model[i] == "split_bottle_sec_het_asym_mig" | 
           dat_params$Model[i] == "split_sizechange_anc_het_asym_mig" | 
           dat_params$Model[i] == "split_sizechange_sec_het_asym_mig") {
   dat_reshape[i, 22] <- vars[11]
   dat_reshape[i, 18] <- NA
 }}

df <- dat_reshape[,c(2,3,5,8:22)]
df[,4:18] <- lapply(df[,4:18], as.numeric)

# Convert to long format
order_cols <- c("nu1", "nu2", "nu1F", "nu2F", "T1", "T2", "m12T1", "m21T1",
                "me12T1", "me21T1", "m12T2", "m21T2", "me12T2", "me21T2", "P")
df_long <- df %>%
  pivot_longer(cols = all_of(order_cols),
               names_to = "variable",
               values_to = "value")


# Create the plot
df_long$Pop_short <- as.factor(df_long$Pop_short)
df_long$variable <-factor(df_long$variable, levels = order_cols)
df_clean <- df_long %>% filter(!is.na(value))
for(pop in levels(df_clean$Pop_short)) {
  print(ggplot(df_clean[df_clean$Pop_short == pop,], aes(x = AIC, y = value, color = AIC)) +
    geom_point() +
    #geom_line() +
    facet_grid(variable ~ Model, drop =TRUE, scales = "free") +
    labs(x = "AIC", y = "Value") +
      scale_color_gradient2(low = "forestgreen",
                            mid = "gold",
                            high = "red",
                            midpoint = median(df_clean$AIC[df_clean$Pop_short == pop], na.rm = TRUE),
                            limits = c(min(df_clean$AIC[df_clean$Pop_short == pop], na.rm = TRUE),
                                       max(df_clean$AIC[df_clean$Pop_short == pop], na.rm = TRUE))) +
    theme_bw() +
    ggtitle(pop) +
    theme(text = element_text(size =8),
          strip.background = element_blank(),
          strip.text = element_text(size = 8, face = "bold"),
      strip.text.x = element_text(angle = 70),
      axis.text.x = element_text(angle = 90)
    ))
  ggsave(paste0("../plots/model_summaries/", pop, "_top", num_runs, ".pdf"),
         height = 23, width = 15, units= "cm")
}


