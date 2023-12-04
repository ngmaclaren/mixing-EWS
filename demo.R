## library(igraph)

source("sim-functions.R")
source("calc-functions.R") # will call library(arrangements)

set.seed(123)

                                        # to generate a new graph
## g <- sample_pa(30, m = 2, directed = FALSE, start.graph = make_full_graph(3))
## A <- as_adj(g, "both", sparse = FALSE)
                                        # to use the BA network from the manuscript
A <- as.matrix(read.table("./data/ba.tsv", sep = "\t"))
N <- nrow(A)

system.time(dl <- simulate_doublewell(A)) # about 60 sec on a 64-bit Intel Core i3-5010U CPU @ 2.10GHz

n <- 5
k1 <- round(.1*length(dl$Cs))
k2 <- round(.9*length(dl$Cs))

                                        # default values are stored in the simulate_* formals
L <- formals(simulate_doublewell)$nsamples # 100
bparam <- formals(simulate_doublewell)$bparam # u
u.step <- formals(simulate_doublewell)$u.step # 0.025
                                        # a vector of the bparam values
us <- seq(0, by = u.step, length.out = length(dl$Cs))

                                        # Select the covariance matrices at k1 and k2
C1 <- dl$Cs[[k1]]
C2 <- dl$Cs[[k2]]

                                        # value of d using N nodes
calc_d(calc_mu(1:N, C1), calc_mu(1:N, C2), calc_nu(1:N, C1, L), calc_nu(1:N, C2, L))

                                        # Find a good node set with n=5
ns <- optimize_nodeset(n, C1, C2, L)
                                        # For comparison, generate a random node set
                                        # and compute it's d value
                                        # For convenience, make it the same structure as ns
rns <- list(nodeset = sample(1:N, length(ns$nodeset)))
rns$d <- calc_d(
    calc_mu(rns$nodeset, C1),
    calc_mu(rns$nodeset, C2),
    calc_nu(rns$nodeset, C1, L),
    calc_nu(rns$nodeset, C2, L)
)

                                        # Calculate the EWS (avg. variance)
V.hat_all <- sapply(dl$Cs, function(C) mean(diag(C)))
                                        # Subset the Cs to just our nodes
                                        # use as.matrix to cover the case of one node
V.hat_ns <- sapply(dl$Cs, function(C) mean(diag(as.matrix(C[ns$nodeset, ns$nodeset]))))
V.hat_rns <- sapply(dl$Cs, function(C) mean(diag(as.matrix(C[rns$nodeset, rns$nodeset]))))

                                        # Calculate the Kendall's τ values
cor(us, V.hat_all, method = "kendall")
cor(us, V.hat_ns, method = "kendall")
cor(us, V.hat_rns, method = "kendall") # τ of random node set is high when bparam is homogeneous u

                                        # Optimize node set size
system.time(ns_optsize <- optimize_nodeset_size(C1, C2, L, "stochastic")) # about 5 sec
length(ns_optsize$nodeset) # 16

                                        # Greedy
ns_greedy <- optimize_nodeset_greedy(n, C1, C2, L)
setdiff(ns_greedy$nodeset, ns$nodeset) # greedy ns has these nodes (1, 15, 38), stochastic ns does not
setdiff(ns$nodeset, ns_greedy$nodeset) # vice versa (7, 10, 29)

system.time(ns_optsize_greedy <- optimize_nodeset_size(C1, C2, L, "greedy")) # 0.377 sec
length(ns_optsize_greedy$nodeset) # 26
