library(parallel)
ncores <- detectCores() - 1
source("calc-functions2.R")

datadir <- "/projects/academic/naokimas/neil/early-warning-theory/"
datafiles <- list.files(datadir, pattern = ".RData")

extract <- function(result, metric = c("d", "largesd")) {
    stopifnot(length(metric) == 1)

    if(metric == "d") {
        idxs <- sapply(result$d, which.max)
        data.frame(
            metric = metric,
            n = 1:6,
            tau = sapply(1:6, function(n) result$tau[[n]][idxs[n]])
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
    tau_largesd <- do.call(rbind, lapply(results, extract, "largesd"))

    taus <- rbind(tau_d, tau_largesd)

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

write.csv(tautau_results, "./data/tau-tau-data.csv", row.names = FALSE)
