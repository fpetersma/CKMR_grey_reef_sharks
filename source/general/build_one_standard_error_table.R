load("~/University of St Andrews/PhD/Sharks/Code/CKMR_grey_reef_sharks/source/result_summaries/sd_cv_both_species.RData")
load("~/University of St Andrews/PhD/Sharks/Code/CKMR_grey_reef_sharks/data/simulation_study/simple/simulation_1000_schemes_all_scenarios_fit_results_sim=all_no_growth.RData")

## Extract failed fit scenarios
conv <- sapply(scenario_fits, function(scen) {
  all(sapply(scen, function(fit) fit$message) == "relative convergence (4)")
})

scenarios_to_keep <- (1:25)[conv] # convert TRUE/FALSE to indices

scen_names <- paste0(rep(paste0(rep("ME", 5), c("-67", "-33", "+0", "+33", "+67")), 
                         each = 5), ":",
                     rep(paste0(rep("GC", 5), c("-10", "-5", "+0", "+5", "+10")), 
                         times = 5))

## Tables for LaTeX
caption <- ""

kable(cbind(ci_coverage_burnham_simple, ci_coverage_burnham_complex[, -1]), 
      booktabs = T, 
      format = "latex", 
      digits = 2,
      col.names = c("Scenario", 
                    "Adult males", "Adult females", 
                    "Adult males", "Adult females"), 
      row.names = F, linesep = "", caption = caption) %>%
  add_header_above( c(" " = 1, "Simple species" = 2, "Complex species" = 2))