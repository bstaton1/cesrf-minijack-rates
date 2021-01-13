# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT THAT CONDUCTS THE MINIJACK ANALYSIS FOUND IN THE MAIN TEXT #
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

# NOTE: THE BOOTSTRAP IN THIS SCRIPT TAKES APPROXIMATELY 2.2 HOURS TO RUN

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

