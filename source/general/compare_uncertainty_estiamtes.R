# compare_uncertainty_estimates.R

load("data/simulation_study/empirical_uncertainty_estimates.RData")
load("data/simulation_study/mle_uncertainty_estimates.RData")

colnames(mle_cv) <- colnames(emp_cv)
colnames(mle_sd) <- colnames(emp_sd)

rel_diff_cv <- (emp_cv[, 2:9] - mle_cv[, 2:9])  / mle_cv[, 2:9] * 100
rel_diff_cv <- cbind(mle_cv[, 1], rel_diff_cv)

## Now extract stuff
