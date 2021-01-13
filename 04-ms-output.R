# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
# SCRIPT TO PRODUCE FIGURES/TABLES SHOWING UP IN MANUSCRIPT #
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

# clear the workspace
rm(list = ls(all = T))

# load necessary packages
source("00-packages.R")

# read/format data file
source("01-data-prep.R")

# load in the functions
source("02-functions.R")

# directory for output
out_dir = "ms-output"
if (!dir.exists(out_dir)) dir.create(out_dir)

# resolution for figures
ppi = 600

# read the bootstrap output objects: created in "03-main-analysis.R"
boot_out_14 = readRDS("model-output/boot_out_14.rds")
boot_out_15 = readRDS("model-output/boot_out_15.rds")
boot_out_16 = readRDS("model-output/boot_out_16.rds")
boot_out_null_14 = readRDS("model-output/boot_out_null_14.rds")
boot_out_null_15 = readRDS("model-output/boot_out_null_15.rds")
boot_out_null_16 = readRDS("model-output/boot_out_null_16.rds")

##### CREATE TABLE SUMMARIZING BOOTSTRAP OUTPUT #####

# summarize the bootstrap output from each year separately
ests_14 = summarize_boot(boot_out_14, boot_out_null_14); ests_14 = cbind(year = 2014, ests_14)
ests_15 = summarize_boot(boot_out_15, boot_out_null_15); ests_15 = cbind(year = 2015, ests_15)
ests_16 = summarize_boot(boot_out_16, boot_out_null_16); ests_16 = cbind(year = 2016, ests_16)

# combine into one 
ests = rbind(ests_14, ests_15, ests_16)

# reorder columns
ests = ests[,c("year", "type", "id", "estimate", "lwr95", "upr95", "p_val")]

# save to a csv file
write.csv(ests, file.path(out_dir, "bootstrap-summaries.csv"), row.names = F)

