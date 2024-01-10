### Need to decide if remove plotting functions to a different file
### Can happen later.

## μ is the average variance of nodes within a node set.
## C is a covariance matrix, so this function takes the mean of the main diagonal of the submatrix of
## C which includes only the rows and columns corresponding to the selected nodes.
## See Eq. (8) in the main text and the Methods subsection "Calculation and performance assessment of
## early warning signals".
calc_mu <- function(nodes, C) mean(diag(as.matrix(C[nodes, nodes])))

## ν is the variance of variances of nodes within a node set.
## Required inputs are: the nodes in the set (selected by index), the covariance matrix, and L, which
## is the number of samples from x_i(t).
## See Eq. (11) in the main text and the Methods subsection "Calculation and performance assessment of
## early warning signals".
calc_nu <- function(nodes, C, L) (2*sum(as.matrix(C[nodes, nodes])^2))/((length(nodes)^2)*(L - 1))

## Find the distance between two distributions of variances (with mean variance μ and variance of
## variance ν):
## the absolute difference between the two means, scaled by the square root of the sum of the two
## variances.
## Use calc_mu and calc_nu to compute the variance means and variances.
calc_d <- function(mu1, mu2, nu1, nu2) {
    abs(mu2 - mu1)/sqrt(nu1 + nu2)
}

## This .name version of calc_d is a wrapper for calc_d with arguments that make sense for the
## analyses in the main text.
.calc_d <- function(nodes, k1, k2, Cs, L) {
    calc_d(
        calc_mu(nodes, Cs[[k1]]),
        calc_mu(nodes, Cs[[k2]]),
        calc_nu(nodes, Cs[[k1]], L),
        calc_nu(nodes, Cs[[k2]], L)
    )
}

## Chooses node sets. If the number of possible combinations is below `maxn`, will find the solution
## exactly. Otherwise, takes a random sample of `maxn` different possible node sets. 
## n is the size of the node set.
## N is the number of nodes in the network (and the number of rows/columns in the covariance matrix)
## The function will output a n x ? matrix, the columns of which are unique combinations of node
## indices. First, the function finds out how many possible unique sets of size n nodes there are
## given the total number of nodes. If that number is less than `maxn`, it returns all combinations.
## If the number is larger than `maxn`, it returns `maxn` random subsets of the possible combinations.
## Nodes are in rows, node sets are in columns. 
chooser <- function(n, N, maxn = 5000) {
                                        # will need to adjust syntax on cluster, here or somewhere else
    require(arrangements, lib.loc = "/user/neilmacl/rlocal/") 

    census_n <- choose(N, n)
    
    if(census_n < maxn) {
        ## combn(1:N, n)
        t(combinations(1:N, n)) # to make the same dims as combn()
    } else {
        ## combn(1:N, n)[, sample(1:census_n, maxn)]
        t(combinations(1:N, n, nsample = maxn))
    }
}

## Given a desired node set size `n`, two covariance matrices `C1` and `C1`, and the number of samples
## from x_i(t) used to compute the covariance matrices `L`, return the best node set of size `n` along
## with its `d` value.
optimize_nodeset <- function(n, C1, C2, L, maxn = 5000) {
    stopifnot(nrow(C1) == nrow(C2))
    stopifnot(ncol(C1) == ncol(C2))
    stopifnot(isSymmetric(C1))
    stopifnot(isSymmetric(C2))
    
    N <- nrow(C1)
    nodesets <- chooser(n, N, maxn)

    ds <- apply(nodesets, 2, function(nodes) {
        mu1 <- calc_mu(nodes, C1)
        mu2 <- calc_mu(nodes, C2)
        nu1 <- calc_nu(nodes, C1, L)
        nu2 <- calc_nu(nodes, C2, L)
        calc_d(mu1, mu2, nu1, nu2)
    })

    best <- which.max(ds)
    nodeset <- as.numeric(nodesets[, best])
    d <- ds[best]

    return(list(nodeset = sort(nodeset), d = d))
}

## Given a possible set of nodes, return the single node that has the largest d
node_selector <- function(possible, C1, C2, L) {
    ds <- sapply(possible, function(node) {
        calc_d(C1[node, node], C2[node, node], calc_nu(node, C1, L), calc_nu(node, C2, L))
    })
    
    possible[which.max(ds)]
}

## A greedy version of optimize_nodeset(). Beginning with the set of all nodes, select the node
## associated with the largest d, add it to the node set, and remove it from the list of possible
## nodes. Continue until n nodes are in the nodeset. 
optimize_nodeset_greedy <- function(n, C1, C2, L) {
    stopifnot(nrow(C1) == nrow(C2))
    stopifnot(ncol(C1) == ncol(C2))
    stopifnot(isSymmetric(C1))
    stopifnot(isSymmetric(C2))

    N <- nrow(C1)

    nodeset <- numeric(n)
    possible <- seq_len(N)
    
    i <- 1
    while(i <= n) {
        bestnode <- node_selector(possible, C1, C2, L)
        nodeset[i] <- bestnode
        possible <- possible[-which(possible == bestnode)]
        i <- i + 1
    }

    d <- calc_d(
        calc_mu(nodeset, C1),
        calc_mu(nodeset, C2),
        calc_nu(nodeset, C1, L),
        calc_nu(nodeset, C2, L)
    )

    return(list(nodeset = sort(nodeset), d = d))
}

optimize_nodeset_size <- function(C1, C2, L, alg = c("stochastic", "greedy"), control = list()) {
    alg <- match.arg(alg, c("stochastic", "greedy"))
    func <- switch(
        alg,
        stochastic = optimize_nodeset,
        greedy = optimize_nodeset_greedy
    )

    if(alg == "stochastic") require(arrangements)

    current_d <- 0
    last_d <- -1
    n <- 1
    nodeset <- numeric(n)

    if("tolerance" %in% names(control)) tolerance <- control$tolerance else tolerance <- 0.01

    while(current_d > (1 + tolerance)*last_d) {
        ## print(c(n, current_d, last_d))
        if("maxn" %in% names(control)) {
            result <- func(n, C1, C2, L, control$maxn)
        } else {
            result <- func(n, C1, C2, L)
        }

        last_d <- current_d
        current_d <- result$d
        n <- n + 1
    }

    return(result)
}
    
### Functions specific to our analysis
analyze_stable_range <- function(rng, Xs, Cs, Ks, backoff = c(0.1, 0.9)) {
    L <- nrow(Xs[[1]])
    .backoff <- round(max(rng)*backoff)
    k1 <- rng[.backoff[1]]
    k2 <- rng[.backoff[2]]

                                        # collect node sets
    S <- lapply(c(1:5, N), function(x) chooser(x, N))

                                        # d score for all node sets
    S_d <- lapply(S, function(x) apply(x, 2, function(col) .calc_d(col, k1, k2, Cs, L)))
                                        # sample variance of x_i(t) over the L samples at each k∈{1...K}
    Xs_var <- do.call(rbind, lapply(Xs, function(X) apply(X, 2, var)))[rng, ]
                                        # average variance for all node sets
    S_var <- lapply(S, function(x) apply(x, 2, function(col) rowMeans(as.matrix(Xs_var[, col]))))
                                        # Select the n nodes with the largest sd at k2
    S_largesd <- lapply(c(1:5, N), function(n) {
                                        # collect the std devs at k2
        sds <- sqrt(Xs_var[k2, ])
                                        # need node indices of the nodes with the largest sd
        which(sds %in% sort(sds, decreasing = TRUE)[seq(n)])
    })
    d_largesd <- lapply(S_largesd, function(s) .calc_d(s, k1, k2, Cs, L))
    var_largesd <- lapply(S_largesd, function(s) rowMeans(as.matrix(Xs_var[, s])))
    tau_largesd <- lapply(var_largesd, function(v) cor(Ks[rng], v, method = "kendall"))
                                        # Poorly named variable: this is tau of EWS vs. bif param
    S_ews <- lapply(
        S_var,
        function(x) apply(x, 2, function(col) cor(Ks[rng], col, method = "kendall"))
    )

    return(list(S = S, d = S_d, tau = S_ews,
                largesd = list(S = S_largesd, d = d_largesd, var = var_largesd, tau = tau_largesd),
                bpvals = c(Ks[k1], Ks[k2])))
}

analyze_stable_range_alt <- function(rng, Xs, Cs, Ks, backoff = c(0.1, 0.9)) {
    L <- nrow(Xs[[1]])
    .backoff <- round(max(rng)*backoff)
    k1 <- rng[.backoff[1]]
    k2 <- rng[.backoff[2]]

                                        # collect node sets
    S <- lapply(c(1:5, N), function(x) chooser(x, N))

                                        # d score for all node sets
    S_d <- lapply(S, function(x) apply(x, 2, function(col) .calc_d(col, k1, k2, Cs, L)))
                                        # sample sd of x_i(t) over the L samples at each k∈{1...K}
    Xs_sd <- do.call(rbind, lapply(Xs, function(X) apply(X, 2, sd)))[rng, ]
                                        # average sd for all node sets
    S_sd <- lapply(S, function(x) apply(x, 2, function(col) rowMeans(as.matrix(Xs_sd[, col]))))
                                        # Poorly named variable: this is tau of EWS vs. bif param
    S_ews <- lapply(
        S_sd,
        function(x) apply(x, 2, function(col) cor(Ks[rng], col, method = "kendall"))
    )

    return(list(S = S, d = S_d, tau = S_ews, bpvals = c(Ks[k1], Ks[k2])))
}

simanalysis <- function(dl, backoff = NULL, ...) {
    Xs <- dl$Xs
    Cs <- dl$Cs

    if(dyn %in% c("doublewell", "SIS")) {
        Ks <- seq(from = 0, by = stepsize, length.out = length(Xs))
    } else {
        Ks <- switch(
            bparam,
            u = seq(from = 0, by = stepsize, length.out = length(Xs)),
            D = seq(from = 1, by = stepsize, length.out = length(Xs))
        )
    }

    if(is.null(backoff)) backoff <- c(0.1, 0.9)

    analyze_stable_range(rng = seq_along(Xs), Xs = Xs, Cs = Cs, Ks = Ks,
                         backoff = backoff)
}

                                        # This is τ against d
                                        # v2
plot_result <- function(result, simrun, filename = NULL) {
    library(latex2exp, lib.loc = "/user/neilmacl/rlocal/")
    palette("Tableau 10")
    S_d <- result$d
    S_ews <- result$tau
    S_stdev <- result$stdev
    largesd <- result$largesd
    
    neach <- lengths(S_d)
    ptcol <- 1:length(S_d)
    ptcols <- unlist(mapply(function(n1, n2) rep(n1, n2), ptcol, neach))

    ht <- 8
    wd <- 12
    labelsize <- 2.5
    ticksize <- 2
    if(save_plots) {
        pdf(filename, height = ht, width = wd, family = "sans")
    } else dev.new(height = ht, width = wd)
    par(mfrow = c(2, 3), mar = c(4.5, 4.5, 0.5, 0.5) + 0.5)
    xlim <- range(unlist(S_ews))
    ylim <- range(unlist(S_d))
    for(i in seq_along(S_d)) {
        x <- S_ews[[i]]
        y <- S_d[[i]]

        maxd <- which.max(y)
        max_y <- y[maxd]
        max_x <- x[maxd]
        
        plot(
            x, y,
            xlab = TeX(r"($\tau$)", italic = TRUE), ylab = TeX(r"($d$)", italic = TRUE),
            xlim = xlim, ylim = ylim,
            pch = 1, col = 1, cex = 1, lwd = 2,
            cex.axis = ticksize, cex.lab = labelsize#, font.lab = 3
        )
        points(
            x = largesd$tau[[i]],
            y = largesd$d[[i]],
            pch = 4, cex = 2, lwd = 2, col = 2
        )
                                        # Uncomment to compare node selected by greedy algorithm
        ## if(i < length(S_d)) {
        ##     add_greedy_point(i, simrun)
        ## } else {
        ##     add_greedy_point(N, simrun)
        ## }

        abline(v = max_x, col = "gray", lwd = .9, lty = 2)
        abline(h = max_y, col = "gray", lwd = .9, lty = 2)

        if(i == 1) {
            legend("topright", bty = "n", legend = c("Sampled node sets", "Largest SD node set"),
                   col = 1:2, pch = c(1, 4), pt.lwd = 2, cex = 1.75)
        }
        mtext(paste0("(", letters[i], ")"), line = -2.5, adj = 0.02, font = 2, cex = 2)
    }
    if(save_plots) dev.off()
}

## add_greedy_point <- function(n, dl) { # (i, simrun), will need to pass it in from the main function
##     Xs <- dl$Xs
##     Cs <- dl$Cs
##     us <- seq(from = 0, by = stepsize, length.out = length(Xs))

##     rng = seq_along(Xs)
##     backoff <- determine_backoff(Xs)
##     L <- nrow(Xs[[1]])
##     k1 <- rng[backoff[1]] # ith from the start
##     k2 <- rng[length(rng) - (backoff[2]-1)] # ith from the end, including the last

##     Xs_var <- do.call(rbind, lapply(Xs, function(X) apply(X, 2, var)))
    
##     result <- optimize_nodeset_greedy(n, Cs[[k1]], Cs[[k2]], L)
##                                         # tau of node set, so cor(us, rowMeans(Xs_var[, nodeset]...
##     x <- cor(us, rowMeans(as.matrix(Xs_var[, result$nodeset])), method = "kendall")
##     y <- result$d# d of node set

##     points(x, y, pch = 3, cex = 2, lwd = 2, col = 3)
## }
    

### Convenience and reporting functions

get_zscore <- function(result) {
    dvals <- result$d[[1]]
    d <- dvals[which.max(dvals)]
    (d - mean(dvals))/sd(dvals)
}

get_prop <- function(i, result) {
    idx <- which.max(result$d[[i]])

    tau <- result$tau[[i]][idx]
    taus <- result$tau[[i]]

    if(dyn %in% c("doublewell", "SIS")) {
        sum(taus > tau)/length(taus)
    } else {
        sum(taus < tau)/length(taus)
    }
}

get_all_props <- function(results) {
    do.call(
        rbind,
        mclapply(
            results, function(result) sapply(1:5, get_prop, result),
            mc.cores = detectCores() - 1
        )
    )
}

get_taudist <- function(i, result) {
    idx <- which.max(result$d[[i]])

    tau <- result$tau[[i]][idx]
    taus <- result$tau[[i]]
    meantau <- mean(taus)

    if(dyn %in% c("doublewell", "SIS")) {
        maxtau <- max(taus)
        (maxtau - tau)/(maxtau - meantau)
    } else {
        mintau <- min(taus)
        (tau - mintau)/(meantau - mintau)
    }
}

get_all_taudists <- function(results) {
    do.call(
        rbind,
        mclapply(
            results, function(result) sapply(1:5, get_taudist, result),
            mc.cores = detectCores() - 1
        )
    )
}

get_greedy <- function(n, simrun) {
                                        # modify here to expand to other backoff values
    k1 <- floor(.1*length(simrun$Cs))
    k2 <- length(simrun$Cs) - (k1 - 1)
    L <- formals(simulate_doublewell)$nsamples

    C1 <- simrun$Cs[[k1]]
    C2 <- simrun$Cs[[k2]]

    res <- optimize_nodeset_greedy(n, C1, C2, L)

    Xs_var <- do.call(rbind, lapply(simrun$Xs, function(X) apply(X, 2, var)))

    list(
        S = res$nodeset,
        d = res$d,
        tau = cor(seq_along(simrun$Xs), rowMeans(as.matrix(Xs_var[, res$nodeset])), method = "kendall")
    )
}

    
get_prop_greedy <- function(n, simrun, result) {
    greedy <- get_greedy(n, simrun)
    if(greedy > 0) {
        sum(result$tau[[n]] > greedy)/length(result$tau[[n]])
    } else {
        sum(result$tau[[n]] < greedy)/length(result$tau[[n]])
    }
}

get_taudist_greedy <- function(n, simrun, result) {
    greedy <- optimize_nodeset_greedy(n, simrun)
    taus <- result$tau[[n]]
    b <- which.max(abs(taus))

    (taus[b] - greedy)/(taus[b] - mean(taus))
}
    
