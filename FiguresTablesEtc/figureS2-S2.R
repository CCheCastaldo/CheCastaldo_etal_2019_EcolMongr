# ________________________________________________________________________________
# load packages

if('pacman' %in% rownames(installed.packages()) == FALSE) {
  
  install.packages('pacman', repos = "http://cran.case.edu")
  
}

pacman::p_load(tidyverse)

# ________________________________________________________________________________
# load data

load(file = "../Library/hoboPumicePlain.rda")

work6 <- hoboPumicePlain

load(file = "../Library/monthlySnotel.rda") 

work7 <- monthlySnotel %>%
  
  group_by(month) %>%
  
  summarise(rainfall = mean(precip))
  
# ________________________________________________________________________________

temp <- c(work6$temp.m[5], work6$temp.m[6], 0, work6$temp.m[3], work6$temp.m[4], 0, work6$temp.m[1], work6$temp.m[2], 0, work6$temp.m[9], work6$temp.m[10], 0, work6$temp.m[7], work6$temp.m[8])

temp.sd <- c(work6$temp.sd[5], work6$temp.sd[6], 0, work6$temp.sd[3], work6$temp.sd[4], 0, work6$temp.sd[1], work6$temp.sd[2], 0, work6$temp.sd[9], work6$temp.sd[10], 0, work6$temp.sd[7], work6$temp.sd[8])

rh <- c(work6$RH.m[5], work6$RH.m[6], 0, work6$RH.m[3], work6$RH.m[4], 0, work6$RH.m[1], work6$RH.m[2], 0, work6$RH.m[9], work6$RH.m[10], 0, work6$RH.m[7], work6$RH.m[8])

rh.sd <- c(work6$RH.sd[5], work6$RH.sd[6], 0, work6$RH.sd[3], work6$RH.sd[4], 0, work6$RH.sd[1], work6$RH.sd[2], 0, work6$RH.sd[9], work6$RH.sd[10], 0, work6$RH.sd[7], work6$RH.sd[8])

rainfall <- c(work7$rainfall[3], 0, work7$rainfall[2], 0,work7$rainfall[1], 0, work7$rainfall[6], 0, work7$rainfall[5])

# ________________________________________________________________________________
# plot it

pdf(file = "figureS2-S2.pdf", width = 6, height = 8)

cexAll <- 1.1

cexAll2 <- .75

par(mfrow = c(3, 1))

par(pin = c(5, 2), mai = c(.5, 1, .25, 1))

barplot(rainfall, axis.lty = 0, xpd = FALSE, xaxt = "n", ylim = c(0, 25), ylab = expression("Precipitation (cm/mo)"), cex.lab = cexAll, 
  
  cex.axis = cexAll, space = 0, col = "darkgray")

abline(h = 0, lty = 1, lwd = 1)

mtext("Jun", side = 1, line = 0.5, cex = cexAll2, at = 0.5)

mtext("Jul", side = 1, line = 0.5, cex = cexAll2, at = 2.5)

mtext("Aug", side = 1, line = 0.5, cex = cexAll2, at = 4.5)

mtext("Sep", side = 1, line = 0.5, cex = cexAll2, at = 6.5)

mtext("Oct", side = 1, line = 0.5, cex = cexAll2, at = 8.5)

text(x = 0, y = 24, "(a)", cex = cexAll2 + .5, font = 2)


par(pin = c(5, 2), mai = c(.5, 1, .25, 1))

barplot(temp, axis.lty = 0, xpd = FALSE, xaxt = "n", ylim = c(0, 40), ylab = expression("Maximum daily temperature"~(degree~C)), cex.lab = cexAll, cex.axis = cexAll, space = 0, col = c("white", "gray", "white"))

abline(h = 0, lty = 1, lwd = 1)

arrows(c(1:14) * 1 - 0.5, temp - temp.sd, c(1:14) * 1 - 0.5, temp + temp.sd, code = 3, angle = 90, length = 0.05)

mtext("Jun", side = 1, line = 0.5, cex = cexAll2, at = 1)

mtext("Jul", side = 1, line = 0.5, cex = cexAll2, at = 4)

mtext("Aug", side = 1, line = 0.5, cex = cexAll2, at = 7)

mtext("Sep", side = 1, line = 0.5, cex = cexAll2, at = 10)

mtext("Oct", side = 1, line = 0.5, cex = cexAll2, at = 13)

text(x = 0, y = 38, "(b)", cex = cexAll2 + .5, font = 2)

legend(.25, 45, inset = .05, cex = cexAll, title = "", c("upland", "riparian"), fill = c("gray", "white"), bty = "n")

par(pin = c(5, 2), mai = c(.5, 1, .25, 1))

barplot(rh, axis.lty = 0, xpd = FALSE, xaxt = "n", ylim = c(0, 100), ylab = expression("Minimum daily RH (%)"), cex.lab = cexAll, 
  
  cex.axis = cexAll, space = 0, col = c("white", "gray", "white"))

abline(h = 0, lty = 1, lwd = 1)

mtext("Jun", side = 1, line = 0.5, cex = cexAll2, at = 1)

mtext("Jul", side = 1, line = 0.5, cex = cexAll2, at = 4)

mtext("Aug", side = 1, line = 0.5, cex = cexAll2, at = 7)

mtext("Sep", side = 1, line = 0.5, cex = cexAll2, at = 10)

mtext("Oct", side = 1, line = 0.5, cex = cexAll2, at = 13)

arrows(c(1:14) * 1 - 0.5, rh - rh.sd, c(1:14) * 1 - 0.5, rh + rh.sd, code = 3, angle = 90, length = 0.05)

text(x = 0, y = 96, "(c)", cex = cexAll2 + .5, font = 2)

legend(.25, 113, inset = .05, cex = cexAll, title = "", c("upland", "riparian"), fill = c("gray", "white"), bty = "n")

dev.off()