# ________________________________________________________________________________
# load packages

if('pacman' %in% rownames(installed.packages()) == FALSE) {
  
  install.packages('pacman', repos = "http://cran.case.edu")
  
}

pacman::p_load(tidyverse, forecast, lubridate)

# ________________________________________________________________________________
# load soil moisture probe data and snotel data

load(file = "../Library/soilProbes.rda")

load(file = "../Library/dailyPrecipitation.rda")

# ________________________________________________________________________________
# calculate vwc and percent field capacity for control upland sensors

work1 <- soilProbes %>%
  
  filter(treatment=="dry") %>%
  
  group_by(date) %>%
  
  summarise(pfc = mean(pfc))
  
work1$pfc <- ma(work1$pfc, order = 48, centre = TRUE)

work16 <- work1 %>%
  
  filter(strftime(date, format="%H:%M:%S") == "12:00:00") %>%
  
  full_join(dailyPrecipitation, by = "date") %>%

  arrange(date)

# ________________________________________________________________________________
# plot daily rainfall and percent field capacity

pdf(file = "figureS2-S1.pdf", width = 6, height = 8)

cexAll <- .8

par(mfrow = c(3, 1))

par(mar = c(5.1, 4.1, 2.1, 4.1))

work17 <- subset(work16, year(work16$date)==2009 & month(work16$date) %in% c(5, 6, 7, 8, 9, 10))

plot(work17$date, work17$precip, type = "l", lwd = 2, cex.axis = cexAll, xlab = "Day", ylab = "Precipitation (cm/day)", ylim = c(0, 6))

text(x = as.POSIXct("2009-05-04 12:00", format = "%Y-%m-%d %H:%M"), y = 5.8, "(a) 2009", cex = cexAll + .2, font = 2)

work17 <- subset(work16, year(work16$date)==2010 & month(work16$date) %in% c(5, 6, 7, 8, 9, 10))

plot(work17$date, work17$precip, type = "l", lwd = 2, cex.axis = cexAll, xlab = "Day", ylab = "Precipitation (cm/day)", ylim = c(0, 6))

text(x = as.POSIXct("2010-05-04 12:00", format = "%Y-%m-%d %H:%M"), y = 5.8, "(b) 2010", cex = cexAll + .2, font = 2)

par(new = T)

plot(work17$date, work17$pfc, type = "l", lwd = 3, col = "red", axes = F, xlab = NA, ylab = NA, ylim = c(0, 150))

axis(side = 4, at = c(0, 25, 50, 75, 100, 125, 150), cex.axis = cexAll)

mtext(side = 4, line = 3, '% field capacity', cex = cexAll - .1)

abline(h = 100, lwd = 2, lty = 2)

work17 <- subset(work16, year(work16$date)==2011 & month(work16$date) %in% c(5, 6, 7, 8, 9, 10))

plot(work17$date, work17$precip, type = "l", lwd = 2, cex.axis = cexAll, xlab = "Day", ylab = "Precipitation (cm/day)", ylim = c(0, 6))

text(x = as.POSIXct("2011-05-04 12:00", format = "%Y-%m-%d %H:%M"), y = 5.8, "(c) 2011", cex = cexAll + .2, font = 2)

par(new = T)

plot(work17$date, work17$pfc, type = "l", lwd = 3, col = "red", axes = F, xlab = NA, ylab = NA, ylim = c(0, 150))

axis(side = 4, at = c(0, 25, 50, 75, 100, 125, 150), cex.axis = cexAll)

mtext(side = 4, line = 3, '% field capacity', cex = cexAll - .1)

abline(h = 100, lwd = 2, lty = 2)

dev.off()
