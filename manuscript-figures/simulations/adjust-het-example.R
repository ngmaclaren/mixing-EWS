## This file should
## - open ba-doublewell-u-TRUE-TRUE-030 in the main environment
## - open ba-doublewell-u-FALSE-FALSE-005 in a separate environment
## - extract the xrange and yrange from the separate environment
## - make the figure in the main environment, using the extracted xrange and yrange
palette("Tableau 10")
save_plots <- TRUE

net <- "ba"
dyn <- "doublewell"
bparam <- "u"
varyv <- "TRUE"
varys <- "TRUE"
whichrun <- 30

datadir <- "/projects/academic/naokimas/neil/early-warning-theory/"
filestub <- paste(c(net, dyn, bparam, varyv, varys), collapse = "-")
datafile <- paste0(filestub, ".RData")
load(paste0(datadir, datafile))

library(parallel)
library(igraph)
library(latex2exp)
source("calc-functions.R")

filename <- paste0(net, "/", filestub, "-", "030.pdf")

## make the separate environment here
otherfig.env <- new.env() # correct syntax?
with(otherfig.env, {
    net <- "ba"
    dyn <- "doublewell"
    bparam <- "u"
    varyv <- "FALSE"
    varys <- "FALSE"
    whichrun <- 5
    filestub <- paste(c(net, dyn, bparam, varyv, varys), collapse = "-")
    datafile <- paste0(filestub, ".RData")
    load(paste0(datadir, datafile))

    result <- results[[whichrun]]

    S_d <- result$d
    S_ews <- result$tau
    S_stdev <- result$stdev
    largesd <- result$largesd
})
    
## will need to write out the plot_result() code here, adjusting it as necessary
result <- results[[whichrun]]
simrun <- simruns[[whichrun]]

S_d <- result$d
S_ews <- result$tau
S_stdev <- result$stdev
largesd <- result$largesd

neach <- lengths(S_d)
ptcol <- 1:length(S_d)
ptcols <- unlist(mapply(function(n1, n2) rep(n1, n2), ptcol, neach))

ht <- 8
wd <- 12
labelsize <- 2.5
ticksize <- 2
if(save_plots) {
    pdf(filename, height = ht, width = wd, family = "sans")
} else dev.new(height = ht, width = wd)
par(mfrow = c(2, 3), mar = c(4.5, 4.5, 0.5, 0.5) + 0.5)
                                        # changes are here 
xlim <- range(c(with(otherfig.env, range(unlist(S_ews))), range(unlist(S_ews))))
ylim <- range(c(with(otherfig.env, range(unlist(S_d))), range(unlist(S_d))))
for(i in seq_along(S_d)) {
    x <- S_ews[[i]]
    y <- S_d[[i]]

    maxd <- which.max(y)
    max_y <- y[maxd]
    max_x <- x[maxd]
    
    plot(
        x, y,
        xlab = TeX(r"($\tau$)", italic = TRUE), ylab = TeX(r"($d$)", italic = TRUE),
        xlim = xlim, ylim = ylim,
        pch = 1, col = 1, cex = 1, lwd = 2,
        cex.axis = ticksize, cex.lab = labelsize
    )
    points(
        x = largesd$tau[[i]],
        y = largesd$d[[i]],
        pch = 4, cex = 2, lwd = 2, col = 2
    )

    abline(v = max_x, col = "gray", lwd = .9, lty = 2)
    abline(h = max_y, col = "gray", lwd = .9, lty = 2)

    if(i == 1) {
        legend("topright", bty = "n", legend = c("Sampled node sets", "Largest SD node set"),
               col = 1:2, pch = c(1, 4), pt.lwd = 2, cex = 1.75)
    }
    mtext(paste0("(", letters[i], ")"), line = -2.5, adj = 0.02, font = 2, cex = 2)
}
if(save_plots) dev.off()
