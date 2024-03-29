library(optparse)
optionlist <- list(
    make_option(# bifurcation parameter
        c("-b", "--bparam"), type = "character", default = "u",
        help = "Select the bifurcation parameter (default %default). Options are 'u' or 'D'."
    ),
    make_option(# new sims?
        c("--makedata"), action = "store_true", default = FALSE,
        help = "Are new CSVs needed? If flag is present, script runs, produces the CSVs, and quits. If not, it will make and save the figures."
    )
)

args <- parse_args(
    OptionParser(option_list = optionlist),
                                        # Set by hand here for debugging if desired
    ## args = c(
    ##     "--bparam=u", "--makedata", "0.1", "0.7"
    ## ),
                                        #
    positional_arguments = c(0, 2),
    convert_hyphens_to_underscores = TRUE
)

bparam <- args$options$bparam
makedata <- args$options$makedata
backoff <- as.numeric(args$args)
filetag <- paste(backoff, collapse = "_")

library(parallel)
NCORES <- detectCores() - 1
RNGkind("L'Ecuyer-CMRG")
set.seed(42)
library(latex2exp)

source("calc-functions.R")

load("./data/networks.rda")

procedure <- function(datafile) {
    load(paste0(datadir, datafile))
    source("calc-functions2.R", local = TRUE)

    results <- lapply(simruns, simanalysis, backoff = backoff)
    
    props <- get_all_props(results)
    dists <- get_all_taudists(results)

    list(props = props, dists = dists)
}

if(makedata) {
    dynamics <- switch(
        bparam,
        u = c("doublewell", "mutualistic", "genereg", ""),
        D = c("doublewell", "mutualistic", "genereg", "SIS")
    )

    datadir <- "/projects/academic/naokimas/neil/early-warning-theory/"
    datafiles <- list.files(datadir, pattern = paste0("-", bparam, "-"))

                                        # Parallelize here
                                        # mclapply still very slow. Correction: seems to be faster than lapply would be, but still slower than expected. This comment may only apply in an interactive session.
    pdata <- mclapply(datafiles, procedure, mc.cores = NCORES)
    
    names(pdata) <- gsub(".RData", "", datafiles)
    save(pdata, file = paste0("./data/pdata-", bparam, "-", filetag, ".rda"))
    q(save = "no")
} else {
    ## load(paste0("./data/pdata-", bparam, ".rda"))
    load(paste0("./data/pdata-u-", filetag, ".rda"))
    pdata_u <- pdata
    load(paste0("./data/pdata-D-", filetag, ".rda"))
    pdata_D <- pdata
}

palette("Tableau 10")

plotit <- function(plotdata, show_ylabel = TRUE) {
    ns <- 1:5
    p1 <- do.call(rbind, lapply(plotdata, function(x) colMeans(x$props)))
    p2 <- do.call(rbind, lapply(plotdata, function(x) colMeans(x$dists)))
    if(length(plotdata) == 3) ltys <- 1:3 else ltys = c(1, 3)
    matplot(
        ns, 2*t(p1), type = "o", pch = 0, lty = ltys, lwd = 2, cex = 2, col = 1, ylim = c(0, 1.25),
        axes = FALSE, xlab = "", ylab = ""
        ## xlab = "", ylab = TeX(r"($p_1$, $p_2$)", italic = TRUE),
        ## cex.lab = 1.75, cex.axis = 1.75
    )
    box()
    axis(1, cex.axis = 1.75)
    title(xlab = "n", font.lab = 3, cex.lab = 1.75)
    if(show_ylabel) {
        axis(2, cex.axis = 1.75)
        title(ylab = TeX(r"($p_1$, $p_2$)", italic = TRUE), cex.lab = 1.75)
    } else {
        axis(2, tick = TRUE, labels = FALSE)
    }
    matpoints(ns, t(p2), pch = 1, lwd = 2, cex = 2, col = 2)
    matlines(ns, t(p2), lty = ltys, lwd = 2, col = 2)
}

makeplot <- function(nets, labels, filename, singleletter = FALSE) {
    pdf(filename, height = 12, width = 12)
    par(mfrow = c(4, 4), mar = c(4, 4, 0, 0)+0.5)
    idx <- 1
    for(i in 1:2) {# net
        net <- nets[i]
        simgroup <- list(
            u = paste(net, c("doublewell", "mutualistic", "genereg"), sep = "-"),
            D = paste(net, c("doublewell", "mutualistic", "genereg", "SIS"), sep = "-")
        )
        for(j in 1:2) { # bparam
            for(k in 1:4) {
                if(k %in% 2:4) show_ylabel <- FALSE else show_ylabel <- TRUE
                
                if(k == 4 && j == 1 && i == 1) {
                    plot(NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
                    legend(
                        "left", bty = "n", cex = 1.5, pt.cex = 1.75, pt.lwd = 2, lwd = 2,
                        legend = c(
                            "p1", #expression(italic(p)[1]), #TeX(r"($p_1$)", italic = TRUE),
                            "p2", #expression(italic(p)[2]), #TeX(r"($p_2$)", italic = TRUE),
                            "Homogeneous",
                            "Heterogeneous stress",
                            paste(c("Heterogeneous stress", "and noise"), collapse = "\n")
                        ),
                        pch = c(0, 1, NA, NA, NA), lty = c(NA, NA, 1:3), col = c(1, 2, 10, 10, 10)
                    )
                } else if(k == 4 && j == 1 && i == 2) {
                    plot(NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
                } else {
                    if(j == 1) {
                        plotdata <- pdata_u[grep(simgroup$u[k], names(pdata_u))]
                        plotit(plotdata, show_ylabel = show_ylabel)
                        ##mtext(simgroup$u[k])
                        labeler(idx, labels, singleletter)
                        idx <- idx + 1
                    } else {
                        plotdata <- pdata_D[grep(simgroup$D[k], names(pdata_D))]
                        plotit(plotdata, show_ylabel = show_ylabel)
                        ##mtext(simgroup$D[k])
                        labeler(idx, labels, singleletter)
                        idx <- idx + 1
                    }
                }
            }
        }
    }
    dev.off()
}

labeler <- function(idx, labels, singleletter = FALSE) {
    if(singleletter) {
        mtext(labels[idx], line = -2.5, adj = 0.02, font = 2, cex = 2)
    } else {
        mtext(labels[idx], line = -2.5, adj = 0.02, font = 1, cex = 1.5)
    }
}

formain <- c("ba", "chesapeake")
forSI1 <- c("erdosrenyi", "fitness")
forSI2 <- c("catlins", "dolphin")

makeplot(
    formain,
    paste0("(", letters, ")"),
    paste0("./img/overall-summary-maintext-", filetag, ".pdf"),
    singleletter = TRUE
)
makeplot(
    forSI1,
    paste0("(a", letters, ")"),
    paste0("./img/overall-summary-SI1-", filetag, ".pdf"),
    singleletter = TRUE
)
makeplot(
    forSI2,
    paste0("(b", letters, ")"),
    paste0("./img/overall-summary-SI2-", filetag, ".pdf"),
    singleletter = TRUE
)

