library(arrangements)

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
## from x_i(t) `L` used to compute the covariance matrices, return the best node set of size `n` along
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

## Using either the stochastic or greedy algorithms above, increase the size of the node set `n` from 1 until `d` no longer increases by more than a tolerance. Default tolerance is 1%. Use the `control` argument to pass tolerance and `maxn` if desired, like this: control = list(tolerance = 0.005, maxn = 1000). 
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
