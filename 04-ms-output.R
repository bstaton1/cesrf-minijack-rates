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

##### CREATE FIGURE: FIXED-EFFECT MINIJACK RATE BY PROGENY WEIGHT, UNCERTAINTY, AND DATA #####

# refit the null models
fit_null_14 = glmmTMB(minijack ~ progeny_wt + (1|sire_id) + (1|dam_id),
                      data = subset(dat, year == 2014), family = binomial)
fit_null_15 = glmmTMB(minijack ~ progeny_wt + (1|sire_id) + (1|dam_id),
                      data = subset(dat, year == 2015), family = binomial)
fit_null_16 = glmmTMB(minijack ~ progeny_wt + (1|sire_id) + (1|dam_id),
                      data = subset(dat, year == 2016), family = binomial)

# calculate summary statistics for each weight category (for plotting data only)
agg1_14 = aggregate(minijack ~ wt_cat, data = subset(dat, year == 2014), mean); colnames(agg1_14) = c("wt_cat", "p_minijack")
agg2_14 = aggregate(minijack ~ wt_cat,data = subset(dat, year == 2014), length); colnames(agg2_14) = c("wt_cat", "n_smolt")
agg_14 = merge(agg1_14, agg2_14, by = "wt_cat"); agg_14 = cbind(year = 2014, agg_14)

agg1_15 = aggregate(minijack ~ wt_cat, data = subset(dat, year == 2015), mean); colnames(agg1_15) = c("wt_cat", "p_minijack")
agg2_15 = aggregate(minijack ~ wt_cat,data = subset(dat, year == 2015), length); colnames(agg2_15) = c("wt_cat", "n_smolt")
agg_15 = merge(agg1_15, agg2_15, by = "wt_cat"); agg_15 = cbind(year = 2015, agg_15)

agg1_16 = aggregate(minijack ~ wt_cat, data = subset(dat, year == 2016), mean); colnames(agg1_16) = c("wt_cat", "p_minijack")
agg2_16 = aggregate(minijack ~ wt_cat,data = subset(dat, year == 2016), length); colnames(agg2_16) = c("wt_cat", "n_smolt")
agg_16 = merge(agg1_16, agg2_16, by = "wt_cat"); agg_16 = cbind(year = 2016, agg_16)

agg = rbind(agg_14, agg_15, agg_16)

# calculate the midpoint of each weight category
x = stringr::str_split(agg$wt_cat, ",")
x = lapply(x, function(y) stringr::str_remove_all(y, "\\]"))
x = lapply(x, function(y) stringr::str_remove_all(y, "\\("))
agg$mp = unlist(lapply(x, function(y) sum(as.numeric(y))/2))

# set point shapes and sizes based on sample size
agg$pch = ifelse(agg$n_smolt < 30, 22, 21)
agg$cex = ifelse(agg$n_smolt < 30, 1.25, agg$n_smolt * 0.01)

png(file.path(out_dir, "fixed-effect-progeny-wt.png"), h = 6 * ppi, w = 2.75 * ppi, res = ppi)
par(mar = c(1,1,1,1), oma = c(2,3,0,0), mfrow = c(3,1), mgp = c(2,0.35,0), tcl = -0.15, cex.axis = 1.2)

# 2014 plot
plot(p_minijack ~ mp, data = subset(agg, year == 2014), las = 1, type = "n", xlim = range(dat$progeny_wt) + c(-5,5), ylim = c(0,1))
probs = get_probs(fit_null_14, "contrast", F)
probs = probs[probs$sire_age == 3,]
polygon(x = c(probs$progeny_wt, rev(probs$progeny_wt)),
        y = c(probs$lwr95ci, rev(probs$upr95ci)), col = scales::alpha("grey25", 0.25), border = NA)
lines(lwr95ci ~ progeny_wt, data = probs, col = "grey")
lines(upr95ci ~ progeny_wt, data = probs, col = "grey")
lines(prob ~ progeny_wt, data = probs, col = "black", lwd = 2)
box()
points(p_minijack ~ mp, data = subset(agg, year == 2014), 
       pch = pch, col = "black", bg = scales::alpha("grey25", 0.25), cex = cex)
usr = par("usr"); xdiff = diff(usr[1:2]); ydiff = diff(usr[3:4])
text(usr[2] + xdiff * 0.015, usr[4] - ydiff * 0.125, pos = 2, labels = 2014, font = 2, cex = 1.2)
legend("bottomright", legend = c("<30", 60, 100, 200, 300, 400), title = "# Smolt",
       pt.cex = c(1.25, c(60, 100, 200, 300, 400) * 0.01), pch = c(22, 21, 21, 21, 21, 21), 
       col = "black", pt.bg = scales::alpha("grey25", 0.25), bty = "n")

# 2015 plot
plot(p_minijack ~ mp, data = subset(agg, year == 2015), las = 1, type = "n", xlim = range(dat$progeny_wt) + c(-5,5), ylim = c(0,1))
probs = get_probs(fit_null_15, "contrast", F)
probs = probs[probs$sire_age == 3,]
polygon(x = c(probs$progeny_wt, rev(probs$progeny_wt)),
        y = c(probs$lwr95ci, rev(probs$upr95ci)), col = scales::alpha("grey25", 0.25), border = NA)
lines(lwr95ci ~ progeny_wt, data = probs, col = "grey")
lines(upr95ci ~ progeny_wt, data = probs, col = "grey")
lines(prob ~ progeny_wt, data = probs, col = "black", lwd = 2)
box()
points(p_minijack ~ mp, data = subset(agg, year == 2015), 
       pch = pch, col = "black", bg = scales::alpha("grey25", 0.25), cex = cex)
usr = par("usr"); xdiff = diff(usr[1:2]); ydiff = diff(usr[3:4])
text(usr[2] + xdiff * 0.015, usr[4] - ydiff * 0.05, pos = 2, labels = 2015, font = 2, cex = 1.2)

# 2016 plot
plot(p_minijack ~ mp, data = subset(agg, year == 2016), las = 1, type = "n", xlim = range(dat$progeny_wt) + c(-5,5), ylim = c(0,1))
probs = get_probs(fit_null_16, "contrast", F)
probs = probs[probs$sire_age == 3,]
polygon(x = c(probs$progeny_wt, rev(probs$progeny_wt)),
        y = c(probs$lwr95ci, rev(probs$upr95ci)), col = scales::alpha("grey25", 0.25), border = NA)
lines(lwr95ci ~ progeny_wt, data = probs, col = "grey")
lines(upr95ci ~ progeny_wt, data = probs, col = "grey")
lines(prob ~ progeny_wt, data = probs, col = "black", lwd = 2)
box()
points(p_minijack ~ mp, data = subset(agg, year == 2016), 
       pch = pch, col = "black", bg = scales::alpha("grey25", 0.25), cex = cex)
mtext(side = 1, outer = T, line = 0.75, "Progeny Weight (g)")
mtext(side = 2, outer = T, line = 1.5, "Minijack Rate")
usr = par("usr"); xdiff = diff(usr[1:2]); ydiff = diff(usr[3:4])
text(usr[2] + xdiff * 0.015, usr[4] - ydiff * 0.05, pos = 2, labels = 2016, font = 2, cex = 1.2)
dev.off()

