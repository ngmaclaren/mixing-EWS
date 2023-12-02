## Rewrite this to demonstrate each of the functions in calc-functions.R
## calc_mu
## calc_nu
## calc_d
## chooser
## optimize_nodeset
## node_selector
## optimize_nodeset_greedy
## optimize_nodeset_size

## library(parallel)
## ncores <- detectCores() - 1
## RNGkind("L'Ecuyer-CMRG")

library(igraph, lib.loc = "/user/neilmacl/rlocal")

source("sim-functions.R")
source("calc-functions.R")

g <- largest_component(sample_gnm(10, 20))
A <- as_adj(g, "both", sparse = FALSE)
N <- vcount(g)

system.time({
    dl <- simulate_doublewell(A)# , u.step = 0.01
})

ks <- data.frame(
    k1 = floor(c(.1, .1, .3, .1, .3, .5, .1, .3, .5, .7)*length(dl$Cs)),
    k2 = ceiling(c(.9, .7, .9, .5, .7, .9, .3, .5, .7, .9)*length(dl$Cs))
)
percentiles <- data.frame(
    k1 = c(.1, .1, .3, .1, .3, .5, .1, .3, .5, .7),
    k2 = c(.9, .7, .9, .5, .7, .9, .3, .5, .7, .9)
)

n <- 5
k1 <- floor(.1*length(dl$Cs))
k2 <- length(dl$Cs) - (k1 - 1)
L <- formals(simulate_doublewell)$nsamples

C1 <- dl$Cs[[k1]]
C2 <- dl$Cs[[k2]]
## result <- optimize_nodeset(n, C1, C2, L)
## result <- optimize_nodeset_greedy(n, C1, C2, L)

print("Value of d using N nodes")
optimize_nodeset(N, C1, C2, L)$d
print("")
print("Optimizing n, greedy alg")
temp <- optimize_nodeset_size(C1, C2, L, "greedy")
print(temp$nodeset)
print(length(temp$nodeset))
print(temp$d)
print("")
print("Optimizing n, stochastic alg")
temp <- optimize_nodeset_size(C1, C2, L, "stochastic", control = list(maxn = 5000))
print(temp$nodeset)
print(length(temp$nodeset))
print(temp$d)

Xs_var <- do.call(rbind, lapply(dl$Xs, function(X) apply(X, 2, var)))

lapply(
    seq_len(nrow(ks)),
    function(i) {
        result <- optimize_nodeset(n, dl$Cs[[ks[i, 1]]], dl$Cs[[ks[i, 2]]], L)
        nodeset <- result$nodeset
        tau <- cor(seq_along(dl$Cs), rowMeans(Xs_var[, nodeset]), method = "kendall")
        c(k1 = percentiles[i, 1], k2 = percentiles[i, 2], d = result$d, tau = tau)
    }
)

## Make a plot, the x-axis of which is bparam, up to the val b/f bifurcation
## There are several y-axes. One is state.

## plotX <- rowMeans(do.call(rbind, lapply(dl$Xs, colMeans))) # y-axis 1
plotX <- do.call(rbind, lapply(dl$Xs, colMeans)) # y-axis 1

## Then, avg node set std dev
ns_opt <- result$nodeset
d_opt <- result$d
ns_rand <- replicate(1, sample(1:N, n, FALSE), FALSE)
d_rand <- lapply(ns_rand, function(nodes) {
    calc_d(calc_mu(nodes, C1), calc_mu(nodes, C2), calc_nu(nodes, C1, L), calc_nu(nodes, C2, L))
})

avgsd_opt <- sapply(dl$Xs, function(X) mean(apply(X[, ns_opt], 2, sd)))

avgsd_rand <- lapply(
    ns_rand,
    function(nodes) {
        sapply(dl$Xs, function(X) mean(apply(X[, nodes], 2, sd)))
    }
)

xrange <- range(as.numeric(plotX))
sdrange <- range(unlist(c(avgsd_opt, avgsd_rand)))
plot1_ylim <- c(
    xrange[1],
    (sdrange[2]/sdrange[1])*xrange[1]
)

## now, plot
pdf("./img/demonstration.pdf", width = 8)
par(mar = c(4, 4, 1, 4))
## plot(seq_along(dl$Xs), plotX, type = "l", lty = 1, lwd = 3, col = 1,
##      xlab = "Currently the bifurcation parameter index, will fix later", ylab = "x*")
matplot(seq_along(dl$Xs), plotX, type = "l", lty = 1, lwd = .5, col = 1,
        ylim = plot1_ylim,
        xlab = "Currently the bifurcation parameter index, will fix later", ylab = "x*")
par(new = TRUE)
plot(NULL, xlim = range(seq_along(dl$Xs)), ylim = range(unlist(c(avgsd_opt, avgsd_rand))),
     axes = FALSE, xlab = "", ylab = "")
for(i in seq_along(avgsd_rand)) {
    lines(seq_along(dl$Xs), avgsd_rand[[i]], lty = 1, lwd = 3, col = adjustcolor(3, 1))#lwd = 1, .25
}
rm(i)
lines(seq_along(dl$Xs), avgsd_opt, lty = 1, lwd = 3, col = 2)
axis(4, pretty(range(unlist(c(avgsd_opt, avgsd_rand)))))
mtext("Average standard deviation", 4, 3)
legend("topleft", bty = "n", lwd = 2, col = 1:3,
       legend = c("Node state", "Avg. SD, optimized node set", "Avg. SD, random node set"))
dev.off()
