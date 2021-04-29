# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT THAT CONDUCTS THE MINIJACK ANALYSIS FOUND IN THE MAIN TEXT #
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

# NOTE: THE BOOTSTRAP IN THIS SCRIPT TAKES APPROXIMATELY 2.5 HOURS TO RUN

# load necessary packages
source("00-packages.R")

# read/format data file
source("01-data-prep.R")

# load in the functions
source("02-functions.R")

##### FIT GLMMs #####

# fit GLMMs: fixed effects of sire age and progeny weight, random effects for parent IDs
# separate model for each year
fit_14 = glmmTMB(minijack ~ sire_age + progeny_wt + (1|sire_id) + (1|dam_id),
                 data = subset(dat, year == 2014), family = binomial)
fit_15 = glmmTMB(minijack ~ sire_age + progeny_wt + (1|sire_id) + (1|dam_id),
                 data = subset(dat, year == 2015), family = binomial)
fit_16 = glmmTMB(minijack ~ sire_age + progeny_wt + (1|sire_id) + (1|dam_id),
                 data = subset(dat, year == 2016), family = binomial)

# fit GLMMs: fixed effects of progeny weight, random effects for parent IDs
# separate model for each year
# these are the "null" models, used for testing hypothesis that sire age is important
fit_null_14 = glmmTMB(minijack ~ progeny_wt + (1|sire_id) + (1|dam_id),
                      data = subset(dat, year == 2014), family = binomial)
fit_null_15 = glmmTMB(minijack ~ progeny_wt + (1|sire_id) + (1|dam_id),
                      data = subset(dat, year == 2015), family = binomial)
fit_null_16 = glmmTMB(minijack ~ progeny_wt + (1|sire_id) + (1|dam_id),
                      data = subset(dat, year == 2016), family = binomial)

##### PERFORM PARAMETRIC BOOTSTRAP #####

# start a clock
starttime = Sys.time()

# number of bootstrapped samples per year
nsim = 2000

# number of cores assigned to the cluster
# this assumes your computer has at least 9 cores that can accept jobs
# fewer cores will make the code run more slowly
# run: parallel::detectCores()
# to figure out how many your computer has. leave at least one free for other tasks
ncpus = 10

# initialize a parallel computing cluster
my_cluster = makeSOCKcluster(ncpus)

# send needed packages to the cluster
clusterEvalQ(my_cluster, {library("lme4"); library("glmmTMB"); library("stringr")})

# send needed environmental variables to the cluster
clusterExport(my_cluster, ls())

# set up the random number generator on the cluster
clusterSetupRNG(my_cluster, type = "RNGstream", seed = rep(1, 6))

# perform the bootstrap for 2014: non-null model
cat("\nRunning Bootstrap: 2014 (non-null model)")
boot_out_14 = bootMer(fit_14, FUN = function(rand_fit) c(get_probs(rand_fit), get_odds_ratios(rand_fit)), nsim = nsim,
                      parallel = "snow", cl = my_cluster, ncpus = ncpus, seed = 1)

# perform the bootstrap for 2014: null model
cat("\nRunning Bootstrap: 2014 (null model)")
boot_out_null_14 = bootMer(fit_null_14, FUN = sim_fit_from_null, nsim = nsim,
                           parallel = "snow", cl = my_cluster, ncpus = ncpus, seed = 1)

# perform the bootstrap for 2015: non-null model
cat("\nRunning Bootstrap: 2015 (non-null model)")
boot_out_15 = bootMer(fit_15, FUN = function(rand_fit) c(get_probs(rand_fit), get_odds_ratios(rand_fit)), nsim = nsim,
                      parallel = "snow", cl = my_cluster, ncpus = ncpus, seed = 1)

# perform the bootstrap for 2015: null model
cat("\nRunning Bootstrap: 2015 (null model)")
boot_out_null_15 = bootMer(fit_null_15, FUN = sim_fit_from_null, nsim = nsim,
                           parallel = "snow", cl = my_cluster, ncpus = ncpus, seed = 1)

# perform the bootstrap for 2016: non-null model
cat("\nRunning Bootstrap: 2016 (non-null model)")
boot_out_16 = bootMer(fit_16, FUN = function(rand_fit) c(get_probs(rand_fit), get_odds_ratios(rand_fit)), nsim = nsim,
                      parallel = "snow", cl = my_cluster, ncpus = ncpus, seed = 1)

# perform the bootstrap for 2016: null model
cat("\nRunning Bootstrap: 2016 (null model)")
boot_out_null_16 = bootMer(fit_null_16, FUN = sim_fit_from_null, nsim = nsim,
                           parallel = "snow", cl = my_cluster, ncpus = ncpus, seed = 1)

# stop the cluster
stopCluster(my_cluster)

# end the clock
stoptime = Sys.time()

# calculate the time difference
format(stoptime - starttime, digits = 2)

# save the output objects
if (!dir.exists("model-output")) dir.create("model-output")
saveRDS(boot_out_14, "model-output/boot_out_14.rds")
saveRDS(boot_out_15, "model-output/boot_out_15.rds")
saveRDS(boot_out_16, "model-output/boot_out_16.rds")
saveRDS(boot_out_null_14, "model-output/boot_out_null_14.rds")
saveRDS(boot_out_null_15, "model-output/boot_out_null_15.rds")
saveRDS(boot_out_null_16, "model-output/boot_out_null_16.rds")
