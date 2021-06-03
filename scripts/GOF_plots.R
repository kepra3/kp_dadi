# Title     : Goodness-of-fit plots
# Objective : Plot the bootstrap likelihoods alongside analysis dataset
# Created by: Katharine Prata
# Created on: 26/11/20
#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
sfs <- args[1]
model <- args[2]
emp_ll <- as.numeric(args[3])
emp_chi <- as.numeric(args[4])

# example use of inputs
#sfs <- "AG1-AG2"
#model <- "mig_inbred"
#emp_ll <- -499.09
#emp_chi <- 591.89

########################## Make plots ######################################

bootstraps <- read.delim(paste0("../results/bootstraps/", sfs, "_", model, "_100_nonparametric_bootstraps_vcf.txt"))

ll_seq <- seq(250, 850, 2)
chi_seq <- seq(5, 10, 0.1)

# uncomment if running in R
#hist(-bootstraps$Likelihood, breaks = ll_seq, main = "Bootstrapping Results - Log-likelihood distribution",
#     xlab = "-logL", col = "blue",ylab = "Freq", lty = "blank")
#abline(v = abs(emp_ll), lwd = 3, col = 'red')

#hist(log(bootstraps$chi.squared), breaks = chi_seq, main = "Bootstrapping Results - Chi-squared distribution",
#     xlab = "log(chi2)",col = "blue", ylab = "Freq", lty = "blank")
#abline(v = log(emp_chi), lwd = 3, col = 'red')

########################## Save plots ######################################

pdf(file = paste0("../plots/", sfs, "_", model, "_log-likelihood_100.pdf"), height = 2.25, width = 2.25)
par(mfrow = c(1, 1))
hist(-bootstraps$Likelihood, breaks = ll_seq, xlab = NULL, col = "blue", ylab = NULL, main = NULL, lty = "blank")
abline(v = abs(emp_ll), lwd = 3, col = 'red')
dev.off()

pdf(file = paste0("../plots/", sfs, "_", model, "_chisq_100.pdf"), height = 2.25, width = 2.25)
par(mfrow = c(1, 1))
hist(log(bootstraps$chi.squared), breaks = chi_seq, xlab = NULL, col = "blue", ylab = NULL, main = NULL, lty = "blank")
abline(v = log(emp_chi), lwd = 3, col = 'red')
dev.off()
