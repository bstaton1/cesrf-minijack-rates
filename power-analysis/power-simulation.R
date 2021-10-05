# :::::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT TO RUN STOCHASTIC SIMULATION POWER ANALYSES #
# :::::::::::::::::::::::::::::::::::::::::::::::::: #

# THIS SCRIPT CONTAINS EVERYTHING TO RUN THE POWER ANALYSIS CALCULATIONS
# IT IS INTENDED TO BE CALLED VIA COMMAND LINE SO THAT MULTIPLE PROCESSOR CORES CAN SIMULTANEOUSLY RUN CALCULATIONS

rm(list = ls(all = TRUE))

##### SPECIFY CORE FUNCTION #####

# a custom function to simulate data from one experiment and fit the GLMM
# "experiment" means one year of crosses in which all 4 sire ages are available

# features of the study design that can be altered:
# n_total_sires: the total number of sires in the experiment
# n_total_dams: the total number of dams in the experiment
# mean_progeny_per_cross: the expected number of progeny from a cross (male smolts only)

# features of the true biology that can be altered:
# minijack_probs: the true average probability that a smolt from a sire of a given age will be a minijack
# sigma_sire: SD of sire-specific random effects (logit-normal)
# sigma_sire: SD of dam-specific random effects (logit-normal)

# if fit == TRUE:
# function returns a data.frame with estimation output: estimates, CIs, p-values, etc.

# if fit == FALSE:
# returns the simulated data set

# example argument values
# uncomment and run these to allow running the function line-by-line to see how each part works
# n_total_sires = 34
# n_total_dams = 23
# mean_progeny_per_cross = 22
# minijack_probs = c("1" = 0.50, "3" = 0.40, "4" = 0.30, "5" = 0.20)
# sigma_sire = 0.85
# sigma_dam = 1

simulate_experiment = function(n_total_sires, n_total_dams, mean_progeny_per_cross, minijack_probs, sigma_sire = 0.85, sigma_dam = 1, fit = TRUE) {
  
  ### DEFINE FIXED FEATURES OF EXPERIMENTAL DESIGN ###
  
  # the expected age distribution of sires
  p_sire_age = rep(0.25, 4)
  
  # the expected number of crosses each male is used in
  p_n_cross = c("1" = 0.20, "2" = 0.55, "3" = 0.22, "4" = 0.02, "5" = 0.01)
  
  ### DEFINE FIXED FEATURES OF THE BIOLOGY/FIXED EFFECTS MODEL ###
  
  # mean and variability of progeny weight
  meanLog_progeny_wt = 3.42  # mean(log(progeny_wt))
  sdLog_progeny_wt = 0.28    # sd(log(progeny_wt))
  
  # the log odds ratio of progeny weight on minijack rate
  beta_progeny_wt = 0.27
  
  # determine the coefficients of the logit model
  # based on average progeny wt and provided average minijack rates by sire age
  beta_age1 = qlogis(minijack_probs["1"]) - beta_progeny_wt * exp(meanLog_progeny_wt)
  beta_age3 = qlogis(minijack_probs["3"]) - beta_age1 - beta_progeny_wt * exp(meanLog_progeny_wt)
  beta_age4 = qlogis(minijack_probs["4"]) - beta_age1 - beta_progeny_wt * exp(meanLog_progeny_wt)
  beta_age5 = qlogis(minijack_probs["5"]) - beta_age1 - beta_progeny_wt * exp(meanLog_progeny_wt)
  betas = c(beta_age1, beta_progeny_wt, beta_age3, beta_age4, beta_age5)
  
  # build a general design matrix. one row represents the correct settings for a given sire age
  # which will be used to build the full design matrix 
  # with the progeny_wt column replaced with real values.
  # see below.
  DM_age = cbind("age1" = c(1,1,1,1), 
                 "progeny_wt" = rep(NA, 4),
                 "age3" = c(0,1,0,0),
                 "age4" = c(0,0,1,0),
                 "age5" = c(0,0,0,1))
  rownames(DM_age) = c("1", "3", "4", "5")
  
  ### GENERATE RANDOM EXPERIMENTAL OUTCOMES ###
  
  # build the pool of sires: each will be used in at least one cross
  sire_pool = data.frame(
    sire_id = paste0("sire_", 1:n_total_sires),
    sire_age = sample(c(1,3,4,5), n_total_sires, p_sire_age, replace = TRUE),
    sire_effect = rnorm(n_total_sires, 0, sigma_sire),
    times_used = sample(as.numeric(names(p_n_cross)), n_total_sires, p_n_cross, replace = TRUE)
  )
  sire_pool$sire_age = factor(sire_pool$sire_age, levels = c(1, 3, 4, 5))
  
  # ensure there is at least one sire of each age
  # the rest of the code will fail if one of the ages is missing
  # this is VERY rare, but we are running many simulations
  sire_pool$sire_age[1:4] = c(1,3,4,5)
  
  # build the pool of dams: each should be used in at least one cross, but not guaranteed
  dam_pool = data.frame(
    dam_id = paste0("dam_", 1:n_total_dams),
    dam_effect = rnorm(n_total_dams, 0, sigma_dam)
  )
  
  # build crosses: for each male, randomly assign 'times_used' females to cross with
  crosses = lapply(1:nrow(sire_pool), function(i) {
    tmp = data.frame(
      sire_id = sire_pool$sire_id[i],
      dam_id = sample(dam_pool$dam_id, sire_pool$times_used[i], replace = FALSE)
    )
    tmp$cross_id = paste0(tmp$sire_id, "-", tmp$dam_id)
    tmp
  })
  crosses = do.call(rbind, crosses)
  
  # generate number of progeny per cross (number of binomial trials)
  # assume the number of cross-specific progeny is a Poisson random variable
  # enforce at least two progeny
  crosses$n_progeny = rpois(nrow(crosses), mean_progeny_per_cross)
  crosses$n_progeny[crosses$n_progeny < 2] = 2
  
  # build the full data set: generates 1 row per progeny, 
  # with parent attributes (age, ids, and REs) and individual progeny weight included
  dat = lapply(1:nrow(crosses), function(i) {
    
    tmp = data.frame(
      sire_id = crosses$sire_id[i],
      sire_age = as.character(sire_pool[sire_pool$sire_id == crosses$sire_id[i],"sire_age"]),
      dam_id = crosses$dam_id[i],
      cross_id = crosses$cross_id[i],
      sire_effect = sire_pool$sire_effect[sire_pool$sire_id == crosses$sire_id[i]],
      dam_effect = dam_pool$dam_effect[dam_pool$dam_id == crosses$dam_id[i]],
      progeny_wt = round(rlnorm(crosses$n_progeny[i], meanLog_progeny_wt - 0.5 * sdLog_progeny_wt^2, sdLog_progeny_wt))
    )
    
    # build the design matrix for this cross
    # for obtaining fixed-effect minijack rate for each progeny based on sire age and weight
    DM = DM_age[tmp$sire_age,]; DM[,"progeny_wt"] = tmp$progeny_wt
    
    # produce progeny-specific minijack rate: fixed-effect (sire age & weight) + parent random effects
    tmp$minijack_rate = as.numeric(plogis(DM %*% betas + tmp$sire_effect + tmp$dam_effect))
    
    # add a random minijack status variable
    tmp$minijack = rbinom(nrow(tmp), 1, tmp$minijack_rate)
    
    tmp
  })
  dat = do.call(rbind, dat)
  
  # create 3 data sets, differ only in which sire age is the reference group
  # this is so all contrasts can be obtained without bootstrap
  # bootstrap is too time-intensive for stochastic simulation study
  dat1 = dat; dat1$sire_age = factor(dat1$sire_age, levels = c(1, 3, 4, 5))  # age 1 is the reference
  dat3 = dat; dat3$sire_age = factor(dat3$sire_age, levels = c(3, 1, 4, 5))  # age 3 is the reference
  dat4 = dat; dat4$sire_age = factor(dat4$sire_age, levels = c(4, 1, 3, 5))  # age 4 is the reference
  
  # fit the models. all are identical except for which sire age is the reference group
  if (fit) {
    fit1 = glmmTMB::glmmTMB(minijack ~ progeny_wt + sire_age + (1|sire_id) + (1|dam_id),
                            family = binomial, data = dat1)
    fit3 = glmmTMB::glmmTMB(minijack ~ progeny_wt + sire_age + (1|sire_id) + (1|dam_id),
                            family = binomial, data = dat3)
    fit4 = glmmTMB::glmmTMB(minijack ~ progeny_wt + sire_age + (1|sire_id) + (1|dam_id),
                            family = binomial, data = dat4)
    
    # extract the fixed effects coefficients table from each fit
    coefs1 = summary(fit1)$coef$cond
    coefs3 = summary(fit3)$coef$cond
    coefs4 = summary(fit4)$coef$cond
    
    # obtain the p-values of each contrast
    # THIS METHOD DOES NOT ACCOUNT FOR TYPE I ERROR RATE INFLATION
    p_vals = data.frame(p1v3 = coefs1["sire_age3","Pr(>|z|)"],
                        p1v4 = coefs1["sire_age4","Pr(>|z|)"],
                        p1v5 = coefs1["sire_age5","Pr(>|z|)"],
                        p3v4 = coefs3["sire_age4","Pr(>|z|)"],
                        p3v5 = coefs3["sire_age5","Pr(>|z|)"],
                        p4v5 = coefs4["sire_age5","Pr(>|z|)"]
    )
    
    # obtain fixed-effect fitted values -- estimated minijack rates by sire age at the average progeny weight
    newdata1 = data.frame(sire_age = factor(c(1,3,4,5), levels = levels(dat1$sire_age)),
                          progeny_wt = exp(meanLog_progeny_wt), sire_id = NA, dam_id = NA)
    preds = predict(fit1, newdata1, se.fit = TRUE)
    
    # fitted values
    means = data.frame(
      age_1_mjr_mean = plogis(preds$fit[1]),
      age_3_mjr_mean = plogis(preds$fit[2]),
      age_4_mjr_mean = plogis(preds$fit[3]),
      age_5_mjr_mean = plogis(preds$fit[4])
    )
    
    # lower 95% CL
    lwr95cl = data.frame(
      age1_mjr_lwr95 = plogis(preds$fit[1] + qnorm(0.025) * preds$se.fit[1]),
      age3_mjr_lwr95 = plogis(preds$fit[2] + qnorm(0.025) * preds$se.fit[2]),
      age4_mjr_lwr95 = plogis(preds$fit[3] + qnorm(0.025) * preds$se.fit[3]),
      age5_mjr_lwr95 = plogis(preds$fit[4] + qnorm(0.025) * preds$se.fit[4])
    )
    
    # upper 95% CL
    upr95cl = data.frame(
      age1_mjr_upr95 = plogis(preds$fit[1] + qnorm(0.975) * preds$se.fit[1]),
      age3_mjr_upr95 = plogis(preds$fit[2] + qnorm(0.975) * preds$se.fit[2]),
      age4_mjr_upr95 = plogis(preds$fit[3] + qnorm(0.975) * preds$se.fit[3]),
      age5_mjr_upr95 = plogis(preds$fit[4] + qnorm(0.975) * preds$se.fit[4])
    )
    
    # true values
    true = data.frame(
      age_1_mjr_true = minijack_probs["1"],
      age_3_mjr_true = minijack_probs["3"],
      age_4_mjr_true = minijack_probs["4"],
      age_5_mjr_true = minijack_probs["5"]
    )
    
    # random effect SDs
    re_output = glmmTMB::VarCorr(fit1)
    sigma = data.frame(
      sigma_sire_est = sqrt(re_output$cond$sire_id[1]),
      sigma_dam_est = sqrt(re_output$cond$dam_id[1]),
      sigma_sire_true = sigma_sire,
      sigma_dam_true = sigma_dam
    )
    
    # combine output into a data.frame
    output = cbind(n_crosses = nrow(crosses), true, means, lwr95cl, upr95cl, sigma, p_vals)
    rownames(output) = NULL
    
    # return the output
    return(output)
  } else {
    # if not fitting, just return the data set
    return(dat1)
  }
}

##### BUILD SCENARIOS #####

# study design factor 1: number of sires & dams used: increases the number of crosses
parents1 = c(n_total_sires = 34, n_total_dams = 23)  # baseline: approximate 2015/2016 setting
parents2 = parents1 * 2                              # doubled
parents3 = parents1 * 4                              # quadrupled
parent_scenarios = rbind(Base = parents1, Double = parents2, Quadruple = parents3)

# study design factor 2: number of progeny per cross
progeny1 = 22     # baseline: approximate 2015/2016 setting
progeny2 = 44     # doubled
progeny3 = 88     # quadrupled
progeny_scenarios = c(Base = progeny1, Double = progeny2, Quadruple = progeny3)

# biology factor 1: strength of sire age effect on minijack rate
minijack_probs1 = c("1" = 0.50, "3" = 0.50, "4" = 0.50, "5" = 0.50) # no effect
minijack_probs2 = c("1" = 0.50, "3" = 0.45, "4" = 0.40, "5" = 0.35) # "slight" effect
minijack_probs3 = c("1" = 0.50, "3" = 0.40, "4" = 0.30, "5" = 0.20) # "strong" effect
minijack_scenarios = rbind(None = minijack_probs1, Weak = minijack_probs2, Strong = minijack_probs3)

# biology factor 2: magnitude of random effect variability
sigma1 = c("sire" = 0.85, "dam" = 1)    # average SD across 2014-2016
sigma2 = c("sire" = 0.425, "dam" = 0.5) # half SD
sigma_scenarios = rbind(Base = sigma1, Half = sigma2)

# scenario combos
scenarios = expand.grid(parents = rownames(parent_scenarios),
                        progeny = names(progeny_scenarios),
                        effect = rownames(minijack_scenarios),
                        sigma = rownames(sigma_scenarios))
scenarios = cbind(scenario_ID = 1:nrow(scenarios), scenarios)

##### RUN SIMULATIONS #####

# CALL SCRIPT FROM COMMAND LINE (i.e., terminal) TO RUN MULTIPLE PARALLEL ITERATIONS
# Rscript power-simulation.R 1 10

first_replicate = as.numeric(commandArgs(trailingOnly = TRUE)[1])
n_replicates = as.numeric(commandArgs(trailingOnly = TRUE)[2])

# OR, SET THESE MANUALLY HERE TO RUN INTERACTIVELY
# IT TAKES APPROXIMATELY 40 MINS TO RUN ONE SIMULATION THROUGH ALL 54 SCENARIOS
# first_replicate = 1
# n_replicates = 10  

# find the last replicate ID to run in this script
last_replicate = first_replicate + n_replicates - 1

# container to store all output
output = NULL

# start a timer
starttime = Sys.time()

# loop through scenarios
for (i in scenarios$scenario_ID) {
  
  # print a progress update
  cat("\n", "Scenario: ", i, "\n")
  
  # loop through replicates of this scenario
  for (j in first_replicate:last_replicate) {
    
    # print a progress update
    cat("\r", "Replicate:", j)
    
    # set the random seed for reproducibility: iteration # becomes the seed
    set.seed(j)
    
    # build basic identifier information for this replicate
    temp_id = data.frame(scenario_ID = i, replicate = j)
    
    # carry out simulation/estimation for this replicate using correct scenario settings for each factor
    temp_output = simulate_experiment( 
      n_total_sires = parent_scenarios[scenarios[i,"parents"],"n_total_sires"],
      n_total_dams = parent_scenarios[scenarios[i,"parents"],"n_total_dams"],
      mean_progeny_per_cross = unname(progeny_scenarios[scenarios[i,"progeny"]]),
      minijack_probs = minijack_scenarios[scenarios[i,"effect"],],
      sigma_sire = sigma_scenarios[scenarios[i,"sigma"],"sire"],
      sigma_dam = sigma_scenarios[scenarios[i,"sigma"],"dam"],
      fit = TRUE
    )
    
    # combine ID info with estimation results
    temp_output = cbind(temp_id, temp_output)
    
    # combine results from this replicate with results from previous replicates
    output = rbind(output, temp_output)
  }
}

# stop the timer and determine total time elapsed
stoptime = Sys.time()
cat("\n")
format(round(stoptime - starttime))

# combine with scenario IDs
output = merge(scenarios, output, by = "scenario_ID")

# save output
write.csv(output, paste0("power-simulation-output-", first_replicate, "-", last_replicate, ".csv"), row.names = FALSE)

# CODE TO PLOT THE OUTPUT IS FOUND IN supplemental-material.Rmd
