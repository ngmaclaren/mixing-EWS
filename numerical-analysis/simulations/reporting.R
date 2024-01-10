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
    ##     ## "--network=ba", "--dynamics=genereg", "--bparam=D",
    ##     "--network=ba", "--dynamics=doublewell", "--bparam=u",
    ##     ## "--vary-nodestress", "--vary-noisestrength",
    ##     "--ntrials=50"
    ## ),
    convert_hyphens_to_underscores = TRUE
)

net <- args$network
dyn <- args$dynamics
varyv <- args$vary_nodestress
varys <- args$vary_noisestrength
bparam <- args$bparam
ntrials <- args$ntrials
save_plots <- TRUE

datadir <- "/projects/academic/naokimas/neil/early-warning-theory/"
datafile <- paste0(paste(c(net, dyn, bparam, varyv, varys), collapse = "-"), ".RData")
filestub <- gsub(".RData", "", datafile)
load(paste0(datadir, datafile))

library(parallel)
library(igraph, lib.loc = "/user/neilmacl/rlocal/")
source("calc-functions2.R")

                                        # Create network directory if not present
                                        # Output will be stored in that directory
if(!dir.exists(net)) dir.create(net)
filenames <- paste0(
    net, "/", filestub, "-",
    sprintf("%03d", seq_len(ntrials)),
    ".pdf"
)

                                        # Six-panel plots
for(i in seq_along(results)) {
    if(i %% 5 == 0) plot_result(results[[i]], simruns[[i]], filenames[i])
}

                                        # (1) the overall summaries
props <- 2*as.data.frame(get_all_props(results))
colnames(props) <- paste0("n", 1:5)

tds <- as.data.frame(get_all_taudists(results))
colnames(tds) <- paste0("n", 1:5)

bpvals <- as.data.frame(t(sapply(results, `[[`, "bpvals")))
colnames(bpvals) <- c("low", "high")
bpvals$trial <- sprintf("%03d", seq_len(ntrials))
bpvals <- bpvals[, c("trial", "low", "high")]

                                        # (2) the required data for each saved figure
correlations <- as.data.frame(
    do.call(rbind, lapply(results, function(res) mapply(cor, res$d, res$tau)))
)
colnames(correlations) <- paste0("corr_n", 1:ncol(correlations))
correlations$trial <- sprintf("%03d", seq_len(ntrials))
correlations <- correlations[, c("trial", paste0("corr_n", 1:5))]

rankings <- as.data.frame(
    do.call(
        rbind,
        lapply(
            results, function(res) {
                mapply(
                    function(x, y) {
                        idx <- which.max(x)
                        which(sort(y, decreasing = TRUE) == y[idx])[1]
                    },
                    res$d, res$tau
                )
            }
        )
    )
)
colnames(rankings) <- paste0("rank_n", 1:ncol(rankings))
rankings$trial <- sprintf("%03d", seq_len(ntrials))
rankings <- rankings[, c("trial", paste0("rank_n", 1:5))]
        
sink(paste0(net, "/", filestub, ".txt"))

cat("p1 (this is already multiplied by 2)", sep = "\n")
print(
    data.frame(
        mean = colMeans(props),
        lower = sapply(props, quantile, probs = c(0.025)),
        upper = sapply(props, quantile, probs = c(0.975))
    )
)
cat("\n", "p2", sep = "\n")
print(
    data.frame(
        mean = colMeans(tds),
        lower = sapply(tds, quantile, probs = c(0.025)),
        upper = sapply(tds, quantile, probs = c(0.975))
    )
)
cat("\n\n", "Pearson correlations for each frame of each image file", sep = "\n")
print( correlations)
cat("\n", "Rank of max d node set in terms of tau", sep = "\n")
print(rankings)
cat("\n\n", "Bifurcation parameter values", sep = "\n")
print(bpvals)

sink()
