library(optparse)
optionlist <- list(
    make_option(# network
        c("-g", "--network"), type = "character", default = "ba",
        help = "Choose the network on which to run the simulations. Default is %default, options are 'erdosrenyi', 'fitness', 'ba', 'chesapeake', 'coweeta', 'catlins', 'akatore', and 'dolphin'. For the main simulation runs, we use three model networks ('ba', 'erdosrenyi', and 'fitness') and three empirical networks ('chesapeake', 'catlins', and 'dolphin')."
    ),
    make_option(# dynamics
        c("-d", "--dynamics"), type = "character", default = "doublewell",
        help = "Choose the dynamics to simulate on the network. Default is %default, options are 'doublewell', 'SIS', 'mutualistic', and 'genereg'."
    ),
    make_option( # v switch
        c("-v", "--vary-nodestress"), action = "store_true", default = FALSE,
        help = "Use this flag to allow the node specific stress to vary. Default (i.e., when the flag is not included) is %default. Without the flag, stress for all nodes is the same and begins at 0. With the flag, each node has additional stress drawn from U(0, 1)."
    ),
    make_option( # s switch
        c("-s", "--vary-noisestrength"), action = "store_true", default = FALSE,
        help = "Use this flag to allow each node to experience different dynamical noise strength. Without the flag (i.e., the default is %default), all nodes will experience a noise processes σW, with σ the same for all nodes and specific to the chosen dynamics. With the flag, σ will be adjusted by a a constant drawn from U(-σ, σ)."
    ),
    make_option(# bifurcation parameter
        c("-b", "--bparam"), type = "character", default = "u",
        help = "Select the bifurcation parameter (default %default). Options are 'u' or 'D'."
    ),
    make_option( # ntrials
        c("-n", "--ntrials"), type = "integer", default = 3,
        help = "The number of independent simulations to run (default %default)."
    )   
)

args <- parse_args(
    OptionParser(option_list = optionlist),
                                        # Set by hand here for debugging if desired
    ## args = c(
    ##     "--network=ba", "--dynamics=mutualistic", "--ntrials=19", "--bparam=D",
    ##   , "--vary-nodestress", "--vary-noisestrength"
    ## ),
    convert_hyphens_to_underscores = TRUE
)

net <- args$network
dyn <- args$dynamics
varyv <- args$vary_nodestress
varys <- args$vary_noisestrength
bparam <- args$bparam
ntrials <- args$ntrials

outdir <- "/projects/academic/naokimas/neil/early-warning-theory/"
outfile <- paste0(outdir, paste(c(net, dyn, bparam, varyv, varys), collapse = "-"), ".RData")

library(parallel)
ncores <- detectCores() - 1
                                        # Setting the RNG to L'Ecuyer and declaring a seed should be
                                        # sufficient for reproducible results. `mclapply` uses
                                        # parallel's RNG stream setting behavior by default.
RNGkind("L'Ecuyer-CMRG")
                                        # Sets the RNG for all simulations.
                                        # e.g. s.adj will be proportionally the same for all noise strengths
set.seed(1248)
library(igraph)

source("sim-functions2.R")
source("calc-functions2.R")

load("./data/networks.rda")

                                        # network
g <- networks[[net]]
A <- as_adj(g, "both", sparse = FALSE)
N <- vcount(g)
                                        # dynamics
dynamics <- switch(
    dyn,
    doublewell = simulate_doublewell,
    SIS = simulate_SIS,
    genereg = simulate_genereg,
    mutualistic = simulate_mutualistic
)
stepsize <- switch(
    dyn,
    doublewell = switch(bparam, u = 0.025, D = 0.0025),
    SIS = switch(bparam, u = NULL, D = 0.0025),
    genereg = switch(bparam, u = -0.01, D = -0.01),
    mutualistic = switch(bparam, u = -0.1, D = -0.01)
)
                                        # node-specific stress
u.adj <- switch(
    varyv + 1,
    false = rep(0, N),
    true = runif(N, -0.25, 0.25)
)
                                        # node-specific noise strength
s.adj <- switch(
    varys + 1,
    0, # false
    switch(# true
        dyn,
        doublewell = runif(N, -0.045, 0.045), # .1*.05
        SIS = runif(N, -4.5e-4, 4.5e-4), # .1*5e-4
        genereg = runif(N, -4.5e-6, 4.5e-6), # .1*5e-6
        mutualistic = runif(N, -0.225, 0.225) # ±0.25-.1*0.25; eliminates a problematic near-zero variance effect for this dynamics that occurs on the BA network; the rest are reduced to match
    )
)

                                        # Simulation
simruns <- mclapply(seq(ntrials), function(x) {
    dynamics(A, bparam = bparam, u.adj = u.adj, s.adj = s.adj)
}, mc.cores = ncores)

                                        # Analyze
                                        # The below should be faster, but isn't. At least, not in the interactive session. The sbatch version may be less memory-limited, so leave it in place for that use.
results <- mclapply(
    simruns, simanalysis, backoff = c(0.1, 0.9),
    mc.cores = ncores
)
                                        # Use lapply instead
## results <- lapply(simruns, simanalysis)


                                        # then save
save.image(outfile)
