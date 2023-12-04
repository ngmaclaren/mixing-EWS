# mixing-EWS

Code implementing algorithms in "Anticipating regime shifts by mixing early warning signals from different nodes," by N. Masuda et al. 

This repository includes algorithms described in "Anticipating regime shifts by mixing early warning signals from different nodes," by N. Masuda et al. These algorithms calculate the metric *d*, which quantifies the quality of a node set, and to use this metric to identify good node sets. We also include functions for simulating $x_i$ across a range of a bifurcation parameter for the dynamics discussed in the manuscript.

# Usage

Functions to compute *d* and identify node sets from a network are in `calc-functions.R`. Simulation functions are in `sim-functions.R`. A short demonstration of these procedures is in `demo.R`. We have also included the networks used in the manuscript as adjacency matrices (`.tsv` files) and as an R data file (`networks.rda`) in the `/data` subdirectory.

## Definitions

Want to write the various definitions here. So, define *d*, describe the stochastic and greedy algorithms, and describe the stopping criteria algorithm.

## Example

We have included a file `demo.R` to demonstrate these procedures. We will use the same network here which we use for demonstrations in the manuscript. To run `demo.R`, place the files `demo.R`, `sim-functions.R`, and `calc-functions.R` in a working directory; place the `/data` directory in the same working directory. All of these steps can be accomplished by cloning this repository, for example by downloading and extracting the .zip archive (click the "Code" button, then select "Download ZIP"), or from the command line:

```sh
git clone https://github.com/ngmaclaren/mixing-EWS.git
```

We first simulate the doublewell dynamics on the Barabási-Albert network. Using the functions in `sim-functions.R`, the only required argument is the adjacency matrix.
```R
source("sim-functions.R")
source("calc-functions.R") # will call library(arrangements)

set.seed(123)

A <- as.matrix(read.table("./data/ba.tsv", sep = "\t"))
N <- nrow(A)

system.time(dl <- simulate_doublewell(A)) # about 60 sec on a 64-bit Intel Core i3-5010U CPU @ 2.10GHz

```
The object `dl` is a list of length 2 and includes the samples from $x_i(t)$ ($L = 100$ by default) and the covariance matrices of the $x_i(t)$ at each value of the bifurcation parameter. The function `simulate_doublewell()` simulates the $x_i(t)$ starting from a default state and gradually increasing a bifurcation parameter until at least one node exits its initial basin of attraction. For the double-well dynamics, we use uniform node stress $u$ as the bifurcation parameter. We initially set $u=0$ and $x_i = 1 \forall i$ and increase $u$ by 0.025 in each round of simulations. 

We need to define several parameters in order to find a good node set: the size of the node set (i.e., the number of nodes), the covariance matrices used to compute *d*, and how many samples from $x_i(t)$ were used to compute the covariance matrices. We will use the covariance matrices at the 10th and 90th percentiles of the range of bifurcation parameter values.

```R
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

```
We first compute the value of *d* using all nodes. The inputs are the average variance and the variance of the variances of all nodes at the two bifurcation parameter values we chose above:
```R
calc_d(calc_mu(1:N, C1), calc_mu(1:N, C2), calc_nu(1:N, C1, L), calc_nu(1:N, C2, L))
## 14.02725
```
To find a good node set of size five, we can use the `optimize_nodeset()` function:
```R
ns <- optimize_nodeset(n, C1, C2, L)
```
By default, we use a stochastic algorithm: if the number of possible combinations of nodes in the node set is larger than 5000 (which it is in this case), we randomly select 5000 combinations of nodes (uniformly and without replacement). We compute *d* for each selected node set, and choose the node set which maximizes *d*. The function `optimize_nodeset()` accepts an argument `maxn` that allows a different maximum number of combinations above which sampling occurs. 

We can compare our optimized node set against the node set including all nodes and a random node set using Kendall's τ:

```R
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
## 0.9476932
cor(us, V.hat_ns, method = "kendall")
## 0.8301499
cor(us, V.hat_rns, method = "kendall") 
## 0.7919483
```

We arbitrarily chose $n=5$ above, but we can also select the node set size by iteratively increasing $n$ from 1 until we no longer substantially increase $d$. We use a default tolerance of 0.01 (i.e., we stop increasing $n$ the first time $d$ does not increase by more than 1%):

```R
system.time(ns_optsize <- optimize_nodeset_size(C1, C2, L, "stochastic")) # about 5 sec
length(ns_optsize$nodeset) # 16
```

Until now we have used a stochastic algorithm to select a good node set from among a random sample of node sets. We have also implemented a greedy algorithm. The greedy algorithm starts with $n=1$ and finds the node which maximizes *d*. If $n>1$, for each additional node in the node set, the algorithm finds the node which maximizes *d* among those not yet selected and adds it to the node set. Results will usually differ from the stochastic algorithm.

```R
ns_greedy <- optimize_nodeset_greedy(n, C1, C2, L)
setdiff(ns_greedy$nodeset, ns$nodeset) # greedy ns has these nodes (1, 15, 38), stochastic ns does not
setdiff(ns$nodeset, ns_greedy$nodeset) # vice versa (7, 10, 29)

system.time(ns_optsize_greedy <- optimize_nodeset_size(C1, C2, L, "greedy")) # 0.377 sec
length(ns_optsize_greedy$nodeset) # 26
```

# Dependencies

[arrangements](https://cran.r-project.org/package=arrangements)
