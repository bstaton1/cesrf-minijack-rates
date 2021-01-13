# ::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT THAT HOUSES FUNCTIONS USED IN MAIN-TEXT ANALYSES #
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

##### CREATE_PRED_DATA() #####

# creates a prediction data set

# fit: a model fit object
# progeny_wt: character; "mean" will use the mean progeny weight, any other value will give a sequence from smallest to largest

create_pred_data = function(fit, progeny_wt = "mean") {
  # determine which sire ages were used in this year
  # if they were in 2014, use only ages 3-5. If not 2015, use age 1,3-5
  if (any(substr(fit$frame$sire_id, 1, 3) == "M14")) {
    sire_ages = c(3,4,5)
  } else {
    sire_ages = c(1,3,4,5)
  }
  if (progeny_wt == "mean") {
    pred_progeny_wt = mean(fit$frame$progeny_wt)
  } else {
    pred_progeny_wt = seq(min(fit$frame$progeny_wt), max(fit$frame$progeny_wt), length = 30)
  }
  
  # build a prediction data set
  # may or may not want to include columns for the ids
  pred_data = expand.grid(sire_age = sire_ages, progeny_wt = pred_progeny_wt, sire_id = NA, dam_id = NA)
  
  # return it
  return(pred_data)
}

##### CREATE_CONTRASTS() #####

# based on the sire ages used in fitting the model, create the contrasts of interest
# fit: a fitted model object

create_contrasts = function(fit) {
  # extract the ages used in this model fit
  ages = as.numeric(as.character(create_pred_data(fit)$sire_age))
  
  # build a matrix with all combinations
  m = matrix(NA, length(ages), length(ages))
  colnames(m) = rownames(m) = ages
  
  # build the age names for each cell
  age1 = rep(ages, each = length(ages))
  age2 = rep(ages, length(ages))
  
  # get only unique contrasts
  x = as.logical(lower.tri(m))
  
  # return the data frame with needed contrasts
  data.frame(age1 = age1[x], age2 = age2[x], contrast = paste0(age1[x], "-", age2[x]))
}

##### EXPIT() #####

# performs the inverse logit transformation
# lp: numeric; a value on the logit scale

expit = function(lp) {
  exp(lp)/(1 + exp(lp))
} 

##### ODDS() #####

# calculates odds given a probability
# p: numeric; a value on the probability scale

odds = function(p) {
  p/(1 - p)
} 

##### GET_PROBS() #####

# obtains fitted probabilities from a model
# fit: a fitted model object
# progeny_wt: passed to create_pred_data()
# prob_only: logical; if TRUE, only probabilities will be returned, if FALSE, confidence intervals will be returned

get_probs = function(fit, progeny_wt = "mean", prob_only = T) {
  # produce the data set to predict from
  pred_data = create_pred_data(fit, progeny_wt = progeny_wt)
  
  # produce the predictions
  preds = predict(fit, newdata = pred_data, se.fit = T, re.form = NA)
  
  # create the output
  out = pred_data
  out$prob = expit(preds$fit)
  out$lwr95ci = expit(preds$fit + qnorm(0.025) * preds$se.fit)
  out$upr95ci = expit(preds$fit + qnorm(0.975) * preds$se.fit)
  
  # drop the CIs if requested
  if (prob_only) {
    probs = out$prob
    names(probs) = out$sire_age
    out = probs
  }
  
  # return the output
  return(out)
}

##### GET_ODDS_RATIOS() #####

# based on a model fit, calculate odds ratios between all sire age contrasts
# fit: a fitted model object

get_odds_ratios = function(fit) {
  
  # obtain the mean probability
  probs = get_probs(fit, prob_only = F)
  
  # which contrasts are needed?
  contrasts = create_contrasts(fit)
  contrasts$age1 = as.character(contrasts$age1)
  contrasts$age2 = as.character(contrasts$age2)
  
  # container for odds ratios
  out = numeric(nrow(contrasts))
  
  # calculate odds ratios between each age contrast
  for (i in 1:nrow(contrasts)) {
    out[i] = odds(probs[probs$sire_age == contrasts$age1[i],"prob"])/odds(probs[probs$sire_age == contrasts$age2[i],"prob"])
  }
  
  # assign them names for the contrast type
  names(out) = contrasts$contrast
  
  # return the output
  return(out)
}

##### SIM_FIT_FROM_NULL #####

# simulates data from the null model and fits the non-null model
# fit_null: a fitted model object, the model should not contain a sire_age effect

sim_fit_from_null = function(fit_null) {
  # simulate data out of the null model
  new_response = simulate(fit_null)[[1]][,1]
  
  # determine which year was used to fit model
  yr = str_extract(fit_null$frame$sire_id, "[:digit:][:digit:]-")
  yr = str_remove(yr[1], "-")
  yr = as.numeric(paste0("20", yr))
  
  # add the new response to the data from this year
  new_data = subset(dat, year == yr)
  new_data$minijack = new_response
  
  # fit the non-null model
  new_fit = glmmTMB(minijack ~ sire_age + progeny_wt + (1|sire_id) + (1|dam_id), family = binomial, data = new_data)
  
  # return the summary information
  c(get_probs(new_fit), get_odds_ratios(new_fit))
}

##### GET_P_VAL() #####

# from bootstrapped output, calculate a two-sided p-value
# obs_odds_ratio: the odds ratio from a contrast calculated from original data
# null_odds_ratios: the vector of bootstrapped odds ratios obtained from fitting non-null model to data simulated from the null model

get_p_val = function(obs_odds_ratio, null_odds_ratios) {
  mean(abs(log(null_odds_ratios)) >= abs(log(obs_odds_ratio)))
}

##### SUMMARIZE_BOOT() #####

# given bootstrapped output, calculate relevant summaries
# boot_out: output of running bootMer() on the non-null model
# boot_out_null: output of running bootMer() on the null model

summarize_boot = function(boot_out, boot_out_null) {
  # start the data frame with the names of each quantity
  df = data.frame(id = colnames(boot_out$t))
  
  # add a type: odds ratio or probability
  df$type = ifelse(grepl("-", df$id), "odds_ratio", "probability")
  
  # add the estimate from the original model and data
  df$estimate = unname(boot_out$t0)
  
  # obtain the 95% confidence limits based on bootstrapped data sets
  df$lwr95 = sapply(1:ncol(boot_out$t), function(i) boot.ci(boot_out, index = i, type = "perc")$percent[4])
  df$upr95 = sapply(1:ncol(boot_out$t), function(i) boot.ci(boot_out, index = i, type = "perc")$percent[5])
  
  # calculate p-values of the odds ratios
  ids = df$id[df$type == "odds_ratio"]
  p_val = sapply(ids, function(i) {
    get_p_val(obs_odds_ratio = boot_out$t0[i], null_odds_ratios = boot_out_null$t[,i])
  })
  
  # add it to the data frame
  df$p_val = NA
  df$p_val[df$type == "odds_ratio"] = p_val
  
  # return the output
  return(df)
}
