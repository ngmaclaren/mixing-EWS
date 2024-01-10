library(parallel)
ncores <- detectCores() - 1
source("calc-functions.R")

datadir <- "/projects/academic/naokimas/neil/early-warning-theory/"
datafiles <- list.files(datadir, pattern = ".RData")

datafiles <- datafiles[grep("ba-", datafiles)]

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

    ## results <- mclapply(
    ##     simruns,
    ##     function(x) simanalysis(x),# backoff = determine_backoff(x$Xs)),
    ##     mc.cores = ncores
    ## )

    tau_d <- do.call(rbind, lapply(results, extract, "d"))
    tau_largesd <- do.call(rbind, lapply(results, extract, "largesd"))

    return(list(tau_d = tau_d, tau_largesd = tau_largesd))
}

tautau_results <- lapply(datafiles, collect_result)

summary(
    apply(
        data.frame(a = tautau_results$tau_d[, "tau"], b = tautau_results$tau_largesd[, "tau"]),
        1,
        function(row) row[1] - row[2]
    )
)

