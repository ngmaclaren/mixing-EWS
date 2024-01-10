library(parallel)
ncores <- detectCores() - 1
source("calc-functions2.R")

library(igraph)
load("./data/networks.rda")

datadir <- "/projects/academic/naokimas/neil/early-warning-theory/"
datafiles <- list.files(datadir, pattern = ".RData")

extract_highlow <- function(simrun, model, net) {
                                        # Can assume all nodes are available because of the way the
                                        # simulations are written.
    
    Xs <- simrun$Xs
    ## k <- length(Xs) - (determine_backoff(Xs)[2] - 1) # make sure this is correctly written in the manuscript, or change it to be easier to describe.
    k <- round(0.9*length(Xs)) # hard-coding this to be the same as our standard k2
                                        # I make the node set based on the point closer to the bifurcation
    X <- Xs[[k]]

    g <- networks[[net]]
    N <- vcount(g)
    ##A <- as_adj(g, "both", sparse = FALSE)

    ranks <- sapply(1:N, function(v) {
        sum(colMeans(as.matrix(X[, neighbors(g, v)])))
    })

    df <- data.frame(avail = 1:N, rank = ranks)
    
                                        # from lower to upper: high input---doublewell and SIS
    if(model %in% c("doublewell", "SIS")) {
        df <- df[order(df$rank, decreasing = TRUE), ]
                                        # from upper to lower: low input---genereg and mutualistic
    } else if(model %in% c("genereg", "mutualistic")) {
        df <- df[order(df$rank, decreasing = FALSE), ]
    }

    return(lapply(c(1:5, N), function(n) as.numeric(V(g)[df$avail[1:n]])))
}

extract <- function(result, metric = c("d", "largesd", "input"),
                    simrun = NULL, model = NULL, net = NULL, stepsize = NULL) {
    stopifnot(length(metric) == 1)

    if(metric == "d") {
        idxs <- sapply(result$d, which.max)
        data.frame(
            metric = metric,
            n = 1:6,
            tau = sapply(1:6, function(n) result$tau[[n]][idxs[n]])
        )
    } else if(metric == "input") {
        Xs_var <- do.call(rbind, lapply(simrun$Xs, function(X) apply(X, 2, var)))
        S_input <- extract_highlow(simrun, model, net)
        var_input <- lapply(S_input, function(s) rowMeans(as.matrix(Xs_var[, s])))
        tau_input <- lapply(var_input, function(v) {
            cor(seq(0, by = stepsize, length.out = length(simrun$Xs)), v, method = "kendall")
        })
        data.frame(
            metric = metric,
            n = 1:6,
            tau = unlist(tau_input)
        )
    } else {
        data.frame(
            metric = metric,
            n = 1:6,
            tau = unlist(result$largesd$tau)
        )
        ## idxs <- unlist(result[[metric]])
    }

}

collect_result <- function(datafile, plot = FALSE) {
    load(paste0(datadir, datafile))
    source("calc-functions2.R", local = TRUE)

    strs <- strsplit(datafile, "-")[[1]]
    net <- strs[1]
    model <- strs[2]
    bparam <- strs[3]
    varyv <- strs[4]
    varys <- gsub(".RData", "", strs[5])

    tau_d <- do.call(rbind, lapply(results, extract, "d"))
    ## tau_largesd <- do.call(rbind, lapply(results, extract, "largesd"))
    tau_input <- do.call(
        rbind,
        mapply(
            function(result, simrun) {
                extract(result, "input", simrun = simrun, model = model, net = net, stepsize = stepsize)
            }, results, simruns, SIMPLIFY = FALSE
        )
    )

    taus <- rbind(tau_d, tau_input) #tau_largesd)

    if(plot) {
        dev.new(height = 10, width = 16)
        par(mfrow = c(2, 3), mar = c(4, 4, 1, 1) + 0.5)
        for(n in 1:5) {
            x <- taus$tau[taus$metric == "largesd" & taus$n == n]
            y <- taus$tau[taus$metric == "d" & taus$n == n]
            xlim <- ylim <- range(c(x, y))
            plot(
                x, xlab = "Large SD", y, ylab = "d",
                col = n + 1, pch = 1, cex.axis = 1.5, cex.lab = 1.5,
                xlim = xlim, ylim = ylim
            )
            abline(a = 0, b = 1, lwd = .5, col = 1)
        }
    }

    df <- aggregate(tau ~ metric + n, data = taus, FUN = mean)
    df$net <- net
    df$model <- model
    df$bparam <- bparam
    df$varyv <- varyv
    df$varys <- varys
    
    return(df)
}

## collect_result(datafiles[1])
tautau_results <- do.call(rbind, lapply(datafiles, collect_result))

write.csv(tautau_results, "./data/tau-tau-data-highlowinput.csv", row.names = FALSE)
