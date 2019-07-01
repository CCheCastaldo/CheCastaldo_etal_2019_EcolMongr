library(coda)
library(tikzDevice)

version <- "i-009f9c75a9206a0ac"

#######################################################################################################################################
# load data

load(file = paste0("../ModelBuild/GlobalModel/", version, "/MCMCsamplesCore.rda"))

chains <- rbind(as.matrix(MCMCsamplesCore[[1]]), as.matrix(MCMCsamplesCore[[2]]), as.matrix(MCMCsamplesCore[[3]]))

summary(MCMCsamplesCore)

###############################################################################################################################################################################

CoeffPlot <- function(samples, rows, yspan = c(-5, 7), ylab = NULL, xlab = "values", ymar = 15, cex = .8, col = "black"){

  plotData <- summary(samples)[[2]][rows,]

  par(mar=c(5, ymar, 2, 1))

  plot(c(yspan[1], yspan[2]), c(0.5, length(rows) + 0.5), type = "n", xlab = xlab, ylab = " ", axes = FALSE, cex.lab = cex, cex.axis = cex)

  abline(v = 0, col = "gray")

  axis(1, cex.lab = cex, cex.axis = cex) 

  if (is.null(ylab)) axis(2, at = length(rows):1, labels = row.names(plotData), padj = 0.5, las = 1, cex.lab = cex, cex.axis = cex)

  else axis(2, at = length(rows):1, labels = ylab, padj = 0.5, las = 1, cex.lab = cex, cex.axis = cex)

  segments(x0 = plotData[,1], x1 = plotData[,5], y0 = length(rows):1, y1 = length(rows):1, col = col)

  segments(x0 = plotData[,2], x1 = plotData[,4], y0 = length(rows):1, y1 = length(rows):1, lwd = 3, col = col)

  points(plotData[,3], length(rows):1, col = col)

invisible()

}

###############################################################################################################################################################################

tikz(file = "fig4.tex", width = 6, height = 8)

par(mfrow = c(2,1))

ylab =
  
  c(
    
    "intercept $\\mu_{\\gamma, \\textnormal{riparian}, \\textnormal{male}}$",
    "intercept $\\mu_{\\gamma, \\textnormal{riparian}, \\textnormal{female}}$",
    "intercept $\\mu_{\\gamma, \\textnormal{upland}, \\textnormal{male}}$",
    "intercept $\\mu_{\\gamma, \\textnormal{upland}, \\textnormal{male}}$",
    "reproduction $\\beta_{2, \\textnormal{riparian}, \\textnormal{male}}$",
    "reproduction $\\beta_{2, \\textnormal{riparian}, \\textnormal{female}}$",
    "reproduction $\\beta_{2, \\textnormal{upland}, \\textnormal{male}}$",
    "reproduction $\\beta_{2, \\textnormal{upland}, \\textnormal{female}}$",
    "RGR $\\beta_{3, \\textnormal{riparian}}$",
    "RGR $\\beta_{3, \\textnormal{upland}}$",
    "diameter $\\beta_{4, \\textnormal{riparian}}$",
    "diameter $\\beta_{4, \\textnormal{upland}}$")

rows <- c(25, 26, 27, 28, 4, 5, 6, 7, 12, 13, 16, 17)

co = "black"

CoeffPlot(samples = MCMCsamplesCore, ylab = ylab, rows = rows, yspan = c(-4, 8), xlab = "Posterior credible intervals", ymar = 12, cex = .8, col = co)

mtext("$\\textbf{(a)}$", side = 3, line = 1, at = -10.5, cex = 1)

ylab =
  
  c(
    
    "intercept $\\mu_{\\phi, \\textnormal{riparian}, \\textnormal{male}}$",
    "intercept $\\mu_{\\phi, \\textnormal{riparian}, \\textnormal{female}}$",
    "intercept $\\mu_{\\phi, \\textnormal{upland}, \\textnormal{male}}$",
    "intercept $\\mu_{\\phi, \\textnormal{upland}, \\textnormal{male}}$",
    "reproduction $\\beta_{5, \\textnormal{riparian}, \\textnormal{male}}$",
    "reproduction $\\beta_{5, \\textnormal{riparian}, \\textnormal{female}}$",
    "reproduction $\\beta_{5, \\textnormal{upland}, \\textnormal{male}}$",
    "reproduction $\\beta_{5, \\textnormal{upland}, \\textnormal{female}}$",
    "RGR $\\beta_{6, \\textnormal{riparian}}$",
    "RGR $\\beta_{6, \\textnormal{upland}}$",
    "diameter $\\beta_{7, \\textnormal{riparian}}$",
    "diameter $\\beta_{7, \\textnormal{upland}}$")

rows <- c(29, 30, 31, 32, 8, 9, 10, 11, 14, 15, 18, 19)

co = "black"

CoeffPlot(samples = MCMCsamplesCore, ylab = ylab, rows = rows, yspan = c(-4, 8), xlab = "Posterior credible intervals", ymar = 12, cex = .8, col = co)

mtext("$\\textbf{(b)}$", side = 3, line = 1, at = -10.5, cex = 1)

dev.off()

###############################################################################################################################################################################

pdf(file = "fig4.pdf", width = 6, height = 8)

ylab =
  
  c(
    
    "Intercept (riparian, male)",
    "Intercept (riparian, female)",
    "Intercept (upland, male)",
    "Intercept (upland, female)",
    "Reproduction (riparian, male)",
    "Reproduction (riparian, female)",
    "Reproduction (upland, male)",
    "Reproduction (upland, female)",
    "RGR (riparian)",
    "RGR (upland)",
    "Size (riparian)",
    "Size (upland)")

par(mfrow = c(2,1))

rows <- c(25, 26, 27, 28, 4, 5, 6, 7, 12, 13, 16, 17)

co = "black"

CoeffPlot(samples = MCMCsamplesCore, ylab = ylab, rows = rows, yspan = c(-4, 8), xlab = "Posterior credible intervals", ymar = 12, cex = .8, col = co)

mtext("a", side = 3, line = 0, at = -14.5, cex = 1)
      
rows <- c(29, 30, 31, 32, 8, 9, 10, 11, 14, 15, 18, 19)

co = "black"

CoeffPlot(samples = MCMCsamplesCore, ylab = ylab, rows = rows, yspan = c(-4, 8), xlab = "Posterior credible intervals", ymar = 12, cex = .8, col = co)

mtext("b", side = 3, line = 0, at = -14.5, cex = 1)

dev.off()

