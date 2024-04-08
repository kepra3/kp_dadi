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
    p1 <- ggplot(group[group$AIC < quantile(group$AIC, 0.5),], aes(Model, AIC)) + 
      geom_point() + theme_bw() + ggtitle(paste0(pop, " - Top 50%"))
    print(p1)
    ggsave(paste0("../plots/", pop, ".model-top50.summary.pdf"), p1, height = 5, width = 5,
           units = "cm", dpi = 400)
    print(as.data.frame(group[group$AIC == min(group$AIC[group$Model == "bottle"]),])[1,])
    print(as.data.frame(group[group$AIC == min(group$AIC[group$Model == "bottle_neck"]),])[1,])
    print(as.data.frame(group[group$AIC == min(group$AIC[group$Model == "size_change"]),])[1,])
  }
}

setwd("~/git/kp_dadi/scripts/")

# Load data
results <- read.table("../results/dadi_optimisation_ahya2.txt", sep = "\t",
                      header = TRUE)

results$Pop <- as.factor(results$Pop)
results$Model <- as.factor(results$Model)

table(results$Pop, results$Model)

plot.results(results)

dat <- data.frame(Pop = numeric(0), Model = numeric(0), Opt = numeric(0))
for (group in c("group1", "group2", "group3", "group4", "Amil")) {
  subset <- results[results$Pop == paste0("ahya.fold.het05.", group, ".major"),]
  row <- head(subset[order(subset$AIC),],1)
  row <- row[,c(1,2,9)]
  dat <- rbind(dat, row)
}
