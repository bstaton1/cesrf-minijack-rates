# ::::::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT THAT PREPARES THE RAW DATA FILE FOR ANALYSIS #
# ::::::::::::::::::::::::::::::::::::::::::::::::::: #

# read in the raw data file
dat = read.csv("minijack-data.csv")

# discard age 5 female crosses
dat = subset(dat, dam_age != 5)

# construct a unique cross variable
dat$cross = paste(dat$dam_id, dat$sire_id, sep = "--")

# discard progeny that do not have weight data associated
dat = subset(dat, !is.na(progeny_wt))

# convert the "progeny_mat_status" variable to binary
dat$minijack = ifelse(dat$progeny_mat_status == "Minijack", 1, 0)

# coerce the "sire_age" variable to factor: models will treat it as a categorical variable
dat$sire_age = as.factor(dat$sire_age)

# discard any variables that won't be used
discard = c("spawn_date", "dam_age", "dam_wt", "sire_wt", "progeny_mat_status", "cross_type")
dat = dat[,-which(colnames(dat) %in% discard)]

# create a progeny weight category: useful in plotting
dat$wt_cat = cut(dat$progeny_wt, breaks = seq(5, 105, by = 5))

# remove crosses that produced less than 10 progeny
dat = dat[!(dat$cross %in% names(which(table(dat$cross) < 10))),]
