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
    ##     "--network=ba", "--dynamics=doublewell", "--ntrials=50"
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

library(parallel)
NCORES <- detectCores() - 1
library(latex2exp, lib.loc = "/user/neilmacl/rlocal/")

datadir <- "/projects/academic/naokimas/neil/early-warning-theory/"
possfiles <- list.files(datadir, pattern = ".RData")
datafiles <- possfiles[grep(paste(net, dyn, bparam, sep = "-"), possfiles)]
imgfile <- paste0("./img/summaryfig-", net, "-", dyn, ".pdf")

procedure <- function(datafile) {
    load(paste0(datadir, datafile))
    source("calc-functions2.R", local = TRUE)

    props <- get_all_props(results)
    dists <- get_all_taudists(results)

    list(props = props, dists = dists)
}

dl <- lapply(datafiles, procedure)

pdf(imgfile, height = 4, width = 12)
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 0.5, 0.5) + 0.5)
labelsize <- 2.5
ticksize <- 2
palette("Tableau 10")
for(i in seq_along(dl)) {

    props <- dl[[i]]$props
    dists <- dl[[i]]$dists
    
                                        # this is a panel
    plot(
        NULL, xlim = c(.5, 5.5), ylim = c(0, 2),
        xlab = TeX(r"($n$)", italic = TRUE), ylab = TeX(r"($p_1$, $p_2$)", italic = TRUE),
        cex.axis = ticksize, cex.lab = labelsize
    )
    abline(h = 1, lwd = 0.5, lty = 1, col = "black")

    matpoints(
        x = matrix(jitter(rep(1:5 - 0.15, each = nrow(props)), amount = .05), nrow = nrow(props)),
        y = 2*props,
        pch = 0, col = 1, cex = 1
    )

    points(
        x = 1:5 - 0.15, y = 2*colMeans(props),
        pch = 0, col = 1, cex = 5, lwd = 2.5
    )

    matpoints(
        x = matrix(jitter(rep(1:5 + 0.15, each = nrow(dists)), amount = .05), nrow = nrow(dists)),
        y = dists,
        pch = 1, col = 2, cex = 1
    )

    points(
        x = 1:5 + 0.15, y = colMeans(dists),
        pch = 1, col = 2, cex = 5, lwd = 2.5
    )
    
    mtext(paste0("(", letters[i], ")"), line = -2.5, adj = 0.02, font = 2, cex = 2)

    if(i == 1) {
        legend("topright", bty = "n", col = 1:2, pch = c(0, 1), pt.cex = 2, pt.lwd = 2,
               legend = c(TeX(r"($p_1$)", italic = TRUE), TeX(r"($p_2$)", italic = TRUE)), cex = 1.75)
    }
}
dev.off()
