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

##### CREATE FIGURE: FIXED-EFFECT MINIJACK RATE BY SIRE AGE, UNCERAINTY, AND DATA SPREAD #####

# convert data into aggregated format. each row is a cross
agg1 = aggregate(minijack ~ cross, data = dat, function(x) length(x)); colnames(agg1)[2] = "male_smolt"
agg2 = aggregate(minijack ~ cross, data = dat, function(x) sum(x)); colnames(agg2)[2] = "minijack_smolt"
agg3 = aggregate(year ~ cross, data = dat, function(x) unique(x)); colnames(agg3)[2] = "year"
agg4 = aggregate(dam_id ~ cross, data = dat, function(x) unique(x))
agg5 = aggregate(sire_id ~ cross, data = dat, function(x) unique(x))
agg6 = aggregate(sire_age ~ cross, data = dat, function(x) unique(x))

# combine these data sets into one
agg = merge(agg1, agg2)
agg = merge(agg, agg3)
agg = merge(agg, agg4)
agg = merge(agg, agg5)
agg = merge(agg, agg6)

# calculate the observed minijack rate for each cross
agg$p_minijack = agg$minijack_smolt/agg$male_smolt

# plotting function: create a blank plot
blank_plot = function(yr) {
  plot(1,1, type = "n", xlim = c(0.5,4.5), ylim = c(0,1), xaxt = "n", las = 1, xlab = "", ylab = "")
  axis(side = 1, at = 1:4, labels = c(1,3:5))
  usr = par("usr"); xdiff = diff(usr[1:2]); ydiff = diff(usr[3:4])
  text(usr[2] + xdiff * 0.015, usr[4] - ydiff * 0.125, pos = 2, labels = yr, font = 2, cex = 1.2)
}

# plotting function: add observed data
add_points = function(age, yr, at_x, col) {
  # at_x = ifelse(age == 1, 1, age - 1)
  agg_sub = subset(agg, sire_age == age & year == yr)
  points(p_minijack ~ I(rep(at_x, nrow(agg_sub)) + runif(nrow(agg_sub), -0.2, 0.2)),
         pch = 16, col = scales::alpha(col, 0.25), data = agg_sub,
         cex = 0.075 * male_smolt)
}

# plotting function: add fixed-effect fitted values
add_fit = function(age, yr, at_x, col) {
  # at_x = ifelse(age == 1, 1, age - 1)
  ests_sub = subset(ests, id == age & type == "probability" & year == yr)
  
  if (nrow(ests_sub) > 0) {
    points(x = at_x, y = ests_sub$estimate, cex = 2, pch = 18, col = col)
    segments(at_x, ests_sub$lwr95, at_x, ests_sub$upr95, col = col, lwd = 1)
  }
}

# plotting function: add a legend
add_legend = function() {
  mult = 0.075
  sizes = c(10, 20, 30, 40)
  legend("topleft", legend = sizes, pt.cex = mult * sizes, pch = 16, col = scales::alpha("grey20", 0.25), bty = "n",
         title = "Male Smolt")
}

# plotting function: add "compact letter display" - significance letters
add_cld = function(yr, alpha = 0.05) {
  ests_sub = subset(ests, year == yr & type == "odds_ratio")
  
  contrasts = ests_sub$id
  reject_null = sapply(1:nrow(ests_sub), function(i) ests_sub$p_val[i] < alpha)
  
  names(reject_null) = contrasts
  
  cld = multcompView::multcompLetters(reject_null)$Letters
  
  usr = par("usr"); xdiff = diff(usr[1:2]); ydiff = diff(usr[3:4])
  
  ages = as.numeric(names(cld))
  at_x = ifelse(ages == 1, 1, ages - 1)
  text(at_x, usr[4] - ydiff * 0.05, labels = cld)
}

# make the plot
ages = c(1,3:5)
years = 2014:2016
at_x = 1:4
cols = c("salmon", "royalblue", "forestgreen", "orange")

png(file.path(out_dir, "fixed-effect-sire-age.png"), h = 6 * ppi, w = 2.75 * ppi, res = ppi)
par(mar = c(1,1,1,1), oma = c(2,3,0,0), mfrow = c(3,1), mgp = c(2,0.35,0), tcl = -0.15, cex.axis = 1.2)
for (y in 1:3) {
  blank_plot(years[y])
  for (a in 1:4) {
    add_points(ages[a], years[y], at_x[a], cols[a])
    add_fit(ages[a], years[y], at_x[a], cols[a])
    if (y == 1 & a == 1) add_legend()
    if (y == 1 & a == 1) text(x = 1, y = 0.2, labels = "No Data", font = 3)
  }
  add_cld(years[y])
}

mtext(side = 1, outer = T, line = 0.75, "Sire Age")
mtext(side = 2, outer = T, line = 1.5, "Minijack Rate")
dev.off()


##### CREATE FIGURE: RANDOM-EFFECT MINIJACK RATE BY SIRE_ID, UNCERTAINTY, AND DATA #####

# refit models
fit_14 = glmmTMB(minijack ~ sire_age + progeny_wt + (1|sire_id) + (1|dam_id),
                 data = subset(dat, year == 2014), family = binomial)
fit_15 = glmmTMB(minijack ~ sire_age + progeny_wt + (1|sire_id) + (1|dam_id),
                 data = subset(dat, year == 2015), family = binomial)
fit_16 = glmmTMB(minijack ~ sire_age + progeny_wt + (1|sire_id) + (1|dam_id),
                 data = subset(dat, year == 2016), family = binomial)

create_male_re_plot = function(fit, yr, legend1, legend2) {
  # extract data the model was fitted to
  frame = fit$frame
  
  # retain only unique sire_ids
  x = frame[!duplicated(frame$sire_id),]
  
  # drop the minijack column
  x = x[,-1]
  
  # drop the progeny column
  x = x[,-2]
  
  # set dam_id to NA
  x$dam_id = NA
  
  # set progeny length the mean for that male
  mn_wt = aggregate(progeny_wt ~ sire_id, data = subset(dat, year == yr), mean)
  x = merge(x, mn_wt, by = "sire_id")
  
  # obtain fitted values that account for male fixed and random effect, but not female
  preds = predict(fit, newdata = x, type = "link", se.fit = T)
  
  x$estimate = expit(preds$fit)
  x$lwr95ci = expit(preds$fit + qnorm(0.025) * preds$se.fit)
  x$upr95ci = expit(preds$fit + qnorm(0.975) * preds$se.fit)
  
  # sort by fitted value
  x = x[order(x$estimate),]
  
  # create a blank plot
  par(mar = c(0,0,1,0), tcl = -0.15, mgp = c(2, 0.35, 0))
  inds = 1:nrow(x)
  plot(x$estimate ~ inds, ylim = c(0,1), xaxt = "n", yaxt = "n", type = "n", xlab = "", ylab = "", las = 1, main = yr)
  
  # set the color for points representing males of different ages
  cols = c("1" = "salmon", "3" = "royalblue", "4" = "forestgreen", "5" = "orange")
  tcols = scales::alpha(cols, 0.4)
  
  # loop through sire_ids, drawing the fitted value and observed proportions
  for (j in 1:nrow(x)) {
    # extract the raw data for this ID
    agg_sub = subset(agg, sire_id == x$sire_id[j])
    
    # draw the point estimate
    points(x$estimate[j] ~ j, pch = 18, col = cols[as.character(x$sire_age[j])], cex = 1)
    
    # draw the approximated 95%ci
    segments(j, x$lwr95ci[j], j, x$upr95ci[j], col = cols[as.character(x$sire_age[j])])
    
    # draw the observed proportions
    points(agg_sub$p_minijack ~ rep(j, nrow(agg_sub)), col = tcols[as.character(agg_sub$sire_age)], pch = 16, cex = 0.075 * agg_sub$male_smolt)
  }
  
  # draw boundaries to help separate individuals
  # abline(v = c(1:nrow(x), nrow(x) + 1) - 0.5, col = "grey90")
  
  usr = par("usr"); xdiff = diff(usr[1:2]); ydiff = diff(usr[3:4])
  # text(x = usr[2], y = usr[4] - ydiff * 0.025, pos = 2, labels = yr, font = 2)
  # 
  # draw a legend if requested
  if (legend1) {
    legend(x = usr[1], y = usr[4], ncol = 2, title = "Male Age", legend = names(cols), pch = 15, pt.cex = 1.5, col = cols, bty = "n")
  }
  if (legend2) {
    mult = 0.075
    sizes = c(10, 20, 30, 40)
    legend(x = usr[1], y = usr[4], ncol = 2, legend = sizes, pt.cex = mult * sizes, pch = 16, col = scales::alpha("grey20", 0.25), bty = "n",
           title = "Male Smolt")
  }
  
  # draw a box
  box()
}

png(file.path(out_dir, "random-effect-sire-id.png"), h = 2.5 * ppi, w = 5.62 * ppi, res = ppi)
par(mfrow = c(1,3), oma = c(2.5,3,0,0.5), lend = "square")
create_male_re_plot(fit_14, 2014, T, F)
axis(side = 2, at = seq(0, 1, 0.2), labels = T, las = 2)
create_male_re_plot(fit_15, 2015, F, T)
create_male_re_plot(fit_16, 2016, F, F)
mtext(side = 1, outer = T, line = 0.5, "Individual Sire")
mtext(side = 1, outer = T, line = 1.35, "(Ranked by Random Effect Size)", cex = 0.75, font = 1)
mtext(side = 2, outer = T, line = 1.65, "Minijack Rate")
dev.off()

##### CREATE FIGURE: RANDOM-EFFECT MINIJACK RATE BY DAM_ID, UNCERTAINTY, AND DATA #####

prep_female_re_for_plot = function(fit, yr, sire_age) {
  
  # extract data the model was fitted to
  frame = fit$frame
  
  # retain only unique dam_ids
  x = frame[!duplicated(frame$dam_id),]
  
  # drop the minijack column
  x = x[,-1]
  
  # drop the progeny_length column
  x = x[,-2]
  
  # set sire_id to NA
  x$sire_id = NA
  
  # calculate mean length of progeny
  age = sire_age
  
  if (yr == 2014 & sire_age == 1) {
    mn_wt = data.frame(dam_id = unique(x$dam_id), progeny_wt = 30)
  } else {
    mn_wt = aggregate(progeny_wt ~ dam_id, data = subset(dat, year == yr & sire_age == age), mean)
  }
  x = merge(x, mn_wt, by = "dam_id")
  
  # obtain fitted values that account for male fixed and random effect, but not female
  preds = predict(fit, newdata = x, type = "link", se.fit = T)
  
  x$estimate = expit(preds$fit)
  x$lwr95ci = expit(preds$fit + qnorm(0.025) * preds$se.fit)
  x$upr95ci = expit(preds$fit + qnorm(0.975) * preds$se.fit)
  
  # sort by fitted value
  x = x[order(x$estimate),]
  x
}

female_re_plot = function(fit, sire_age, yr, legend, yaxis = TRUE) {
  x = prep_female_re_for_plot(fit, yr, sire_age)
  
  if (sire_age == 1 & yr == 2014) {
    x = x[-(1:nrow(x)),]
  }
  
  # create a blank plot
  par(mar = c(0,0,0,0), tcl = -0.15, mgp = c(2, 0.35, 0))
  inds = 1:nrow(x)
  
  if (nrow(x) == 0) {
    plot(1,1, type = "n", xaxt = "n", ylim = c(0,1), xlim = c(0,1), yaxt = "n")
    text(x = 0.5, y = 0.4, "No Data", font = 3)
  } else {
    plot(x$estimate ~ inds, ylim = c(0,1), xaxt = "n", yaxt = "n", type = "n", xlab = "", ylab = "", las = 1)
    
    # set the color for points representing males of different ages
    # cols = c("1" = "grey20", "3" = "grey20", "4" = "grey20", "5" = "grey20")
    # tcols = scales::alpha(cols, 0.4)
    cols = c("1" = "salmon", "3" = "royalblue", "4" = "forestgreen", "5" = "orange")
    tcols = scales::alpha(cols, 0.4)
    
    # loop through sire_ids, drawing the fitted value and observed proportions
    age = sire_age
    for (j in 1:nrow(x)) {
      # extract the raw data for this ID
      agg_sub = subset(agg, dam_id == x$dam_id[j] & sire_age == age)
      
      # draw the point estimate
      if (nrow(agg_sub) > 0) {
        points(x$estimate[j] ~ j, pch = 18, cex = 1, col = cols[as.character(sire_age)])
        
        # draw the approximated 95%ci
        segments(j, x$lwr95ci[j], j, x$upr95ci[j], col = cols[as.character(sire_age)])
      }
      
      # draw the observed proportions
      points(agg_sub$p_minijack ~ rep(j, nrow(agg_sub)), col = tcols[as.character(sire_age)], pch = 16, cex = 0.075 * agg_sub$male_smolt)
    }
    
    # draw boundaries to help separate individuals
    # abline(v = c(1:nrow(x), nrow(x) + 1) - 0.5, col = "grey90")
    
    # usr = par("usr"); xdiff = diff(usr[1:2]); ydiff = diff(usr[3:4])
    # text(x = usr[2], y = usr[4] - ydiff * 0.025, pos = 2, labels = yr, font = 2)
  }
  
  if (yaxis) {
    axis(side = 2, at = seq(0.0, 1, 0.2), labels = T, las = 2, cex.axis = 0.9)
  }
  
  # draw a legend if requested
  if (legend) {
    mult = 0.075
    sizes = c(10, 20, 30, 40)
    legend("top", ncol = 2, legend = sizes, pt.cex = mult * sizes, pch = 16, col = scales::alpha("grey20", 0.25), bty = "n",
           title = "Male Smolt in Cross")
    # legend(x = usr[1] + xdiff * 0.175, y = usr[4], title = "Male Age", legend = names(cols), pch = 15, pt.cex = 1.5, col = cols, bty = "n")
  }
  
  # draw a box
  box()
}

cols = c("1" = "salmon", "3" = "royalblue", "4" = "forestgreen", "5" = "orange")

png(file.path(out_dir, "random-effect-dam-id.png"), h = 4 * ppi, w = 5.62 * ppi, res = ppi)
par(mfrow = c(3,4), oma = c(2.5,3.5,3,3))
for (y in 2014:2016) {
  if (y == 2014) {
    fit = fit_14
  } else {
    if (y == 2015) {
      fit = fit_15
    } else {
      if (y == 2016) {
        fit = fit_16
      }
    }
  }
  for (a in c(1,3,4,5)) {
    female_re_plot(fit, a, y, ifelse(a == 1 & y == 2014, T, F), ifelse(a == 1, T, F))
    if (y == 2014) {
      mtext(side = 3, outer = F, line = 0.25, a, font = 1, col = cols[as.character(a)])
    }
    if (a == 5) {
      mtext(side = 4, outer = F, line = 0.25, y, font = 1)
    }
  }
}
mtext(side = 1, outer = T, line = 0.5, "Individual Dam", font = 1)
mtext(side = 1, outer = T, line = 1.35, "(Ranked by Random Effect Size)", cex = 0.75, font = 1)
mtext(side = 2, outer = T, line = 2, "Minijack Rate")
mtext(side = 3, outer = T, line = 1.5, "Sire Age", font = 2)
mtext(side = 4, outer = T, line = 2, "Brood Year", font = 2)
dev.off()
