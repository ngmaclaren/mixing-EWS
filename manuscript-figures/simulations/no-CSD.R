library(parallel)
ncores <- detectCores() - 1
RNGkind("L'Ecuyer-CMRG")
set.seed(2) # given s = 0.094, discarded one random seed for not enough sample (wanted at least 30)
library(igraph)
library(latex2exp)
library(sfsmisc)

source("sim-functions.R")
source("calc-functions.R")
load("./data/networks.rda")

save_plots <- TRUE

net <- "ba"
dyn <- "doublewell"
bparam <- "u"

g <- networks[[net]]
A <- as_adj(g, "both", sparse = FALSE)
N <- vcount(g)
r <- c(1, 2, 5) # for noise-induced transition
x.init <- rep(r[1], N)
nsamples <- 100
spacing <- 1
equil_time <- 100
D <- 0.05
u.init <- rep(0, N)
s <- 0.094
ntimesteps <- 80*200/DeltaT # approx num. of sims * length of typical sim / Δt

X <- doublewell(x.init, A, r, D, u.init, s, ntimesteps, DeltaT)
first_transition <- which(apply(X, 1, max) > 2)[1] - 1

## check for transition here
## pdf("./img/noCSD-noise_transition.pdf")
## matplot(##seq_len(first_transition + 100), X[seq_len(first_transition + 100), ],
##     1:40000, X[1:40000, ],
##     type = "l", lty = 1, lwd = .25, col = "black",
##     xlab = "t", ylab = "x"
## )
## dev.off()

tipping_node <- which.max(X[first_transition + 1, ])
samples <- seq(from = first_transition, to = 1, by = -spacing/DeltaT)
.samples <- seq(1, by = 100, to = length(samples))
samplelist <- lapply(.samples, function(x) rev(samples[seq(from = x, by = 1, length.out = 100)]))
                                        # drop the last one so we have complete samples of 100
samplelist <- rev(samplelist[-length(samplelist)])

                                        # add some extra for plotting
plotsamples <- seq(from = first_transition + spacing/DeltaT, by = spacing/DeltaT, length.out = 500)
                                        # concatenating
plotsamples <- unlist(c(samplelist, list(plotsamples)))
                                        # downsampling
plotsamples <- plotsamples[seq(from = 100, to = length(plotsamples), by = 100)]
plotX <- X[plotsamples, ]

Xs <- lapply(samplelist, function(samp) X[samp, ])
Cs <- lapply(Xs, cov)

k1 <- round(.1*length(Xs))
k2 <- round(.9*length(Xs))

Xs_var <- do.call(rbind, lapply(Cs, diag))

                                        # find the maximizer of d
n <- 5
L <- formals(simulate_doublewell)$nsamples
C1 <- Cs[[k1]]
C2 <- Cs[[k2]]
result <- optimize_nodeset(n, C1, C2, L)

                                        # Panel A
ns_opt <- result$nodeset
d_opt <- result$d
avgvar_opt <- rowMeans(Xs_var[, ns_opt])
us <- seq_len(length(Xs) + 5)
labelsize <- ticksize <- 1.75
## pdf("./img/stochastic-forcing-demonstration.pdf", width = 8)
## palette("Tableau 10")
## par(mar = c(4, 5, 1, 4)+0.8)
## matplot(seq_len(nrow(plotX)), plotX, type = "l", lty = 1, lwd = .5, col = "black",
##         xlab = "", ylab = "", axes = FALSE, ylim = c(0.5, 5.5), xaxs = "i", yaxs = "i")
## axis(4, pretty(c(0.5, 5.5)),
##      cex.axis = ticksize)
## mtext(TeX(r"($x_i^*$)"), 4, 3, cex = labelsize)
## par(new = TRUE)
## plot(NULL, xlim = range(us), ylim = c(0, 0.008),
##      xlab = "Sample index, k", ylab = TeX(r"(Average variance, $\hat{V}_S$)"),
##      cex.axis = ticksize, cex.lab = labelsize, xaxs = "i", yaxs = "i")
## lines(us, c(Xs_var[, tipping_node], rep(NA, 5)), lty = 1, lwd = 3, col = 2)
## lines(us, c(avgvar_opt, rep(NA, 5)), lty = 1, lwd = 3, col = 1)
## legend("topleft", bty = "n", lwd = 3, col = c(1, 2, "black"), cex = .75*labelsize,
##        legend = c(TeX(r"($\hat{V}_S$, optimized node set)"),
##                   TeX(r"($\hat{V}_S$, tipping node)"),
##                   TeX(r"($x_i^*$)")))
## dev.off()

                                        # Panel B
sr_result <- analyze_stable_range(seq_along(Xs), Xs, Cs, seq_along(Xs))

save(
    plotX, ticksize, labelsize, us, tipping_node, Xs_var, avgvar_opt, sr_result, Xs, Cs,
    file = "./figureS32/stochastic-forcing.rda"
)

## with(
##     list(
##         result = sr_result,
##         simrun = list(Xs, Cs),
##         filename = "./img/stochastic-forcing-sixpanel.pdf"
##     ), {
##         S_d <- result$d
##         S_ews <- result$tau
##         S_stdev <- result$stdev
##         largesd <- result$largesd
##         neach <- lengths(S_d)
##         ptcol <- 1:length(S_d)
##         ptcols <- unlist(mapply(function(n1, n2) rep(n1, n2), ptcol, neach))
##         plotNs <- c(1:5, N)
##         ht <- 8
##         wd <- 12
##         labelsize <- 2.5
##         ticksize <- 2
##         if(save_plots) {
##             pdf(filename, height = ht, width = wd, family = "sans")
##         } else dev.new(height = ht, width = wd)
##         par(mfrow = c(2, 3), mar = c(4.5, 4.5, 0.5, 0.5) + 0.5)
##         xlim <- range(unlist(S_ews))
##         ylim <- range(unlist(S_d))
##         for(i in seq_along(S_d)) {
##             x <- S_ews[[i]]
##             y <- S_d[[i]]
##             maxd <- which.max(y)
##             max_y <- y[maxd]
##             max_x <- x[maxd]
##             plot(
##                 x, y,
##                 xlab = TeX(r"($\tau$)", italic = TRUE), ylab = TeX(r"($d$)", italic = TRUE),
##                 xlim = xlim, ylim = ylim,
##                 pch = 1, col = 1, cex = 1, lwd = 2,
##                 cex.axis = ticksize, cex.lab = labelsize#, font.lab = 3
##             )
##             abline(v = max_x, col = "gray", lwd = .9, lty = 2)
##             abline(h = max_y, col = "gray", lwd = .9, lty = 2)
##             mtext(paste("n =", plotNs[i]), line = -2.5, adj = 0.02, font = 1, cex = .75*labelsize)
##         }
##         if(save_plots) dev.off()
##     }
## )

### From here, do similar simulations but with a step function in u. I get to decide how long the simulation is.
set.seed(2)
r <- c(1, 3, 5)
s <- 0.05
u.pulse <- 5

u <- 0
W <- preallocate_noise(s, N, ntimesteps)
x <- x.init
X <- matrix(0, nrow = ntimesteps, ncol = N)
for(timestep in seq_len(ntimesteps)) {
    X[timestep, ] <- x
    if(timestep == 10000/DeltaT) u <- u.pulse
    x <- x +
        (-(x - r[1])*(x - r[2])*(x - r[3]) + D*colSums(A*x) + u)*DeltaT +
        W[timestep, ]*sqrt(DeltaT)
}
first_transition2 <- which(apply(X, 1, max) > 2)[1] - 1

samples <- seq(from = first_transition2, to = 1, by = -spacing/DeltaT)
.samples <- seq(1, by = 100, to = length(samples))
samplelist <- lapply(.samples, function(x) rev(samples[seq(from = x, by = 1, length.out = 100)]))
                                        # drop the last one so we have complete samples of 100
samplelist <- rev(samplelist[-length(samplelist)])

                                        # add some extra for plotting
plotsamples <- seq(from = first_transition2 + spacing/DeltaT, by = spacing/DeltaT, length.out = 1000)
                                        # concatenating
plotsamples <- unlist(c(samplelist, list(plotsamples)))
                                        # downsampling
plotsamples <- plotsamples[seq(from = 100, to = length(plotsamples), by = 100)]
plotX <- X[plotsamples, ]

Xs <- lapply(samplelist, function(samp) X[samp, ])
Cs <- lapply(Xs, cov)

k1 <- round(.1*length(Xs))
k2 <- round(.9*length(Xs))

Xs_var <- do.call(rbind, lapply(Cs, diag))

                                        # find the maximizer of d
n <- 5
L <- formals(simulate_doublewell)$nsamples
C1 <- Cs[[k1]]
C2 <- Cs[[k2]]
result <- optimize_nodeset(n, C1, C2, L)

                                        # Panel C
ns_opt <- result$nodeset
d_opt <- result$d
avgvar_opt <- rowMeans(Xs_var[, ns_opt])
us <- seq_len(length(Xs) + 10)
labelsize <- ticksize <- 1.75
## pdf("./img/sudden-stress-demonstration.pdf", width = 8)
## palette("Tableau 10")
## par(mar = c(4, 5, 1, 4)+0.8)
## matplot(us, plotX, type = "l", lty = 1, lwd = .5, col = "black",
##         xlab = "", ylab = "", axes = FALSE, ylim = c(0.5, 6), xaxs = "i", yaxs = "i")
## axis(4, pretty(c(0.5, 5.5)),
##      cex.axis = ticksize)
## mtext(TeX(r"($x_i^*$)"), 4, 3, cex = labelsize)
## par(new = TRUE)
## plot(NULL, xlim = range(us),
##      ylim = c(0, 0.006),
##      xlab = "Sample index, k", ylab = TeX(r"(Average variance, $\hat{V}_S$)"),
##      cex.axis = ticksize, cex.lab = labelsize, xaxs = "i", yaxs = "i")
## lines(us, c(Xs_var[, tipping_node], rep(NA, 10)), lty = 1, lwd = 3, col = 2)
## lines(us, c(avgvar_opt, rep(NA, 10)), lty = 1, lwd = 3, col = 1)
## legend("topleft", bty = "n", lwd = 3, col = c(1, 2, "black"), cex = .75*labelsize,
##        legend = c(TeX(r"($\hat{V}_S$, optimized node set)"),
##                   TeX(r"($\hat{V}_S$, tipping node)"),
##                   TeX(r"($x_i^*$)")))
## dev.off()

                                        # Panel D
sr_result <- analyze_stable_range(seq_along(Xs), Xs, Cs, seq_along(Xs))
save_plots <- TRUE

save(
    plotX, ticksize, labelsize, us, tipping_node, Xs_var, avgvar_opt, sr_result, Xs, Cs,
    file = "./figureS32/sudden-stress.rda"
)

## with(
##     list(
##         result = sr_result,
##         simrun = list(Xs, Cs),
##         filename = "./img/sudden-stress-sixpanel.pdf"
##     ), {
##         S_d <- result$d
##         S_ews <- result$tau
##         S_stdev <- result$stdev
##         largesd <- result$largesd
##         neach <- lengths(S_d)
##         ptcol <- 1:length(S_d)
##         ptcols <- unlist(mapply(function(n1, n2) rep(n1, n2), ptcol, neach))
##         plotNs <- c(1:5, N)
##         ht <- 8
##         wd <- 12
##         labelsize <- 2.5
##         ticksize <- 2
##         if(save_plots) {
##             pdf(filename, height = ht, width = wd, family = "sans")
##         } else dev.new(height = ht, width = wd)
##         par(mfrow = c(2, 3), mar = c(4.5, 4.5, 0.5, 0.5) + 0.5)
##         xlim <- range(unlist(S_ews))
##         ylim <- c(0, 3.4)#range(unlist(S_d))
##         for(i in seq_along(S_d)) {
##             x <- S_ews[[i]]
##             y <- S_d[[i]]
##             maxd <- which.max(y)
##             max_y <- y[maxd]
##             max_x <- x[maxd]
##             plot(
##                 x, y,
##                 xlab = TeX(r"($\tau$)", italic = TRUE), ylab = TeX(r"($d$)", italic = TRUE),
##                 xlim = xlim, ylim = ylim,
##                 pch = 1, col = 1, cex = 1, lwd = 2,
##                 yaxt = "n",
##                 cex.axis = ticksize, cex.lab = labelsize#, font.lab = 3
##             )
##             axis(2, axTicks(2)[c(TRUE, FALSE)], cex.axis = ticksize)
##             abline(v = max_x, col = "gray", lwd = .9, lty = 2)
##             abline(h = max_y, col = "gray", lwd = .9, lty = 2)
##             mtext(paste("n =", plotNs[i]), line = -2.5, adj = 0.02, font = 1, cex = .75*labelsize)
##         }
##         if(save_plots) dev.off()
##     }
## )
