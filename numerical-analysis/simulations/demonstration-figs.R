library(parallel)
ncores <- detectCores() - 1
RNGkind("L'Ecuyer-CMRG")
set.seed(12345)
library(igraph, lib.loc = "/user/neilmacl/rlocal/")
library(latex2exp, lib.loc = "/user/neilmacl/rlocal/")
library(sfsmisc)

source("sim-functions2.R")
source("calc-functions2.R")
load("./data/networks.rda")

net <- "ba"
dyn <- "doublewell"
bparam <- "u"

g <- networks[[net]]
A <- as_adj(g, "both", sparse = FALSE)
N <- vcount(g)
r <- c(1, 3, 5)
x.init <- rep(r[1], N)
nsamples <- 100
spacing <- 1
equil_time <- 100
D <- 0.05
u.init <- rep(0, N)
u.step <- 0.025
s <- 0.05

ntimesteps <- determine_ntimesteps(nsamples, spacing, equil_time, DeltaT)
samples <- determine_samples(nsamples, spacing, ntimesteps, DeltaT)

i <- 1
Xs <- list()
Cs <- list()

nlowerstate <- sum(x.init < r[2])
nlowerstates <- numeric()
u <- u.init
condition <- function(nlowerstate, threshold = 0.80) nlowerstate > floor(threshold*N)

while(condition(nlowerstate)) {
    X <- doublewell(x.init, A, r, D, u, s, ntimesteps, DeltaT)[samples, ]
    nlowerstate <- sum(X[nrow(X), ] < r[2])

    Xs[[i]] <- X
    Cs[[i]] <- cov(X)
    nlowerstates[i] <- nlowerstate

    if(i %% 10 == 0) print(max(u))
    i <- i + 1
    u <- u + u.step
}

                                        # Node set
all_ls <- which(nlowerstates == N)
.Xs <- Xs[all_ls]
.Cs <- Cs[all_ls]
n <- 5
k1 <- floor(.1*length(.Cs))
k2 <- length(.Cs) - (k1 - 1)
L <- formals(simulate_doublewell)$nsamples
C1 <- .Cs[[k1]]
C2 <- .Cs[[k2]]
result <- optimize_nodeset(n, C1, C2, L)
Xs_var <- do.call(rbind, lapply(.Cs, diag))

ns_opt <- result$nodeset
d_opt <- result$d
ns_rand <- sample(1:N, n, FALSE)
d_rand <- calc_d(
    calc_mu(ns_rand, C1), calc_mu(ns_rand, C2), calc_nu(ns_rand, C1, L), calc_nu(ns_rand, C2, L)
)

avgvar_opt <- rowMeans(Xs_var[, ns_opt])
avgvar_rand <- rowMeans(Xs_var[, ns_rand])

                                        # Plotting
plotX <- do.call(rbind, lapply(Xs, colMeans))
us <- seq(0, by = u.step, length.out = length(Xs))
cor(us[1:length(avgvar_opt)], avgvar_opt, method = "kendall") # 0.81
cor(us[1:length(avgvar_rand)], avgvar_rand, method = "kendall") # 0.81
labelsize <- ticksize <- 1.75

                                        # Figure 5
save(
    us, avgvar_opt, avgvar_rand, ticksize, labelsize, all_ls, plotX, file = "./figure5/data.rda"
)


## Using the same simulation, demonstrate stopping criteria
opts <- mclapply(1:N, function(n) {
    optimize_nodeset(n, C1, C2, L)
}, mc.cores = ncores)
dvals <- sapply(opts, `[[`, "d")
stoppingcriterion <- dvals > c(-1, 1.01*dvals[-length(dvals)])
ddiff <- (dvals[2:length(dvals)] - dvals[1:(length(dvals)-1)])/(2:N - 1:(N-1))
marker <- which(stoppingcriterion == FALSE)[1]
varvals <- lapply(
    opts,
    function(opt) sapply(.Xs, function(X) mean(apply(as.matrix(X[, opt$nodeset]), 2, var)))
)
taus <- sapply(varvals, cor, us[all_ls], method = "kendall")

markns <- c(5, marker)

                                        # Figure S30
save(N, dvals, markns, taus, ticksize, labelsize, file = "./figureS30/data.rda")

