# Title: summarise dadi results
# Author: Katharine Prata
# Date created: 18/03/24

# Packages
library(ggplot2)

# Functions
plot.results <- function(results) {
  for (pop in levels(results$Pop)) {
    print(pop)
    group <- results[results$Pop == pop,]
    p <- ggplot(group, aes(Model, AIC)) + geom_point() + theme_bw() + ggtitle(pop)
    print(p)
    ggsave(paste0("../plots/", pop, ".model.summary.pdf"), p, height = 5, width = 5,
           units = "cm", dpi = 400)
    p1 <- ggplot(group[group$AIC < quantile(group$AIC, 0.9),], aes(Model, AIC)) + 
      geom_point() + theme_bw() + ggtitle(paste0(pop, " - Top 90%")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    print(p1)
    ggsave(paste0("../plots/", pop, ".AIC.top90.summary.pdf"), p1, height = 15, width = 15,
           units = "cm", dpi = 400)
    p2 <- ggplot(group[group$log.likelihood > quantile(group$log.likelihood, 0.9),],
                 aes(Model, log.likelihood)) + 
      geom_point() + theme_bw() + ggtitle(paste0(pop, " - Top 90%")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    print(p2)
    ggsave(paste0("../plots/", pop, ".log.likelihood.top90.summary.pdf"), p1, height = 15, width = 15,
           units = "cm", dpi = 400)
    print(as.data.frame(group[group$AIC == min(group$AIC),])[1,])
    #for (model in levels(results$Model)) {
    #  print(as.data.frame(group[group$AIC == min(group$AIC[group$Model == model]),])[1,])
    #}
  }
}

setwd("~/git/kp_dadi/scripts/")

# Load data
results <- read.table("../results/dadi_optimisation_combined.txt", sep = "\t",
                      header = TRUE)

results$Pop <- as.factor(results$Pop)
results$Model <- as.factor(results$Model)

table(results$Pop, results$Model)

results <- results[results$AIC > 10000,]

plot.results(results)

dat <- data.frame(Pop = numeric(0), Model = numeric(0), Opt = numeric(0))
for (group in sort(unique(results$Pop))) {
  subset <- results[results$Pop == group,]
  row <- head(subset[order(subset$AIC),],3)
  row <- row[,c(1,3,5,6,10, 11)]
  dat <- rbind(dat, row)
}
