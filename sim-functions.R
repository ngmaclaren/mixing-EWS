## Define a default Δt value to use in all integration
DeltaT <- 0.01

## Generate a matrix of Gaussian white noise with sd = `strength` for use in simulating an SDE.
## Generating the random numbers ahead of time saves computation time.
## The matrix has ntimesteps rows and nnodes columns; total simulation time is determined elsewhere.
## Accepts either homogeneous or node-specific noise strength, depending on the length of `strength`.
## Noise will be added to the computed derivative in each time step by looking up values from this
## matrix.
preallocate_noise <- function(strength, nnodes, ntimesteps) {
    if(length(strength) == 1) {
        W <- matrix(rnorm(nnodes*ntimesteps, sd = strength), ncol = nnodes)
    } else {
        stopifnot(length(strength) == nnodes)
        W <- do.call(
            cbind,
            lapply(strength, function(x) rnorm(ntimesteps, sd = x))
        )
    }

    return(W)
}

## Determine the correct number of time steps for a simulation
## Input the number of samples to be taken from the end of the time series, the spacing (in user time
## units) at which those samples should be taken, and a desired equilibration time. DeltaT is set
## above.
## Output is a scalar. Used to preallocate noise and in the integration loop. 
determine_ntimesteps <- function(nsamples, spacing, equil_time, DeltaT) {
    (equil_time/DeltaT) + ((nsamples*spacing)/DeltaT)
}

## Given the number of desired samples, sample spacing and the output of determine_ntimesteps(),
## returns the indices of the samples in the SDE simulation output matrix.
determine_samples <- function(nsamples, spacing, ntimesteps, DeltaT) {
    seq(from = ntimesteps, by = -spacing/DeltaT, length.out = nsamples)
}

### Models and drivers

## Each model has a function that simulates the model given some parameters and a driver function that
## provides a user interface, some defaults, and calls the simulation function over a succession of
## control parameter values.
##
## Although there is a `direction` argument for the driver functions, the simulations are only done in a standard direction for each model: from the lower basin for doublewell and SIS, upper basin for mutualistic and genereg.


doublewell <- function(x.init, A, r, D, u, s, ntimesteps, DeltaT) {
    N <- length(x.init) #N is available in the top environment, but not necessarily inside the function
    stopifnot(N == ncol(A)) #x.init needs to have the same length (for each node) as A has columns
    stopifnot(isSymmetric(unname(A))) # This function uses a calculation that requires a symmetric matrix

    W <- preallocate_noise(s, N, ntimesteps) # generate the Gaussian noise

    x <- x.init # initial conditions
    X <- matrix(0, nrow = ntimesteps, ncol = N) # preallocate the matrix
    for(timestep in seq_len(ntimesteps)) {
                                        # fills initial conditions and the state of the system as
                                        # computed by the previous iteration.
        X[timestep, ] <- x 
        x <- x + # current state
                                        # Add the deterministic part.
                                        # Note that D*colSums(<matrix>) is correct, but only for
                                        # symmetric matrices and because the x in A*x is the x_js only
            (-(x - r[1])*(x - r[2])*(x - r[3]) + D*colSums(A*x) + u)*DeltaT + 
                                        # Add the noise part, looking up a value from the preallocated
                                        # matrix and scaling by the square root of Δt. 
            W[timestep, ]*sqrt(DeltaT) 
    }

    return(X) # This is the whole simulation history, with time in rows and node states in columns.
}

## The driver function for doublewell(). The adjacency matrix, `A`, comes from the network. This is
## the only required argument---the rest has defaults chosen for the paper.
simulate_doublewell <- function(A, bparam = "u", N = ncol(A), # arguments are read sequentially
                                r = c(1, 3, 5), 
                                x.init = rep(r[1], N),
                                u.adj = rep(0, N), s.adj = 0,
                                u.init = rep(0, N) + u.adj, u.step = 0.025,
                                D = switch(bparam, u = 0.05, D = 0),
                                D.step = 0.0025, s = 0.05 + s.adj,
                                nsamples = 100, spacing = 1, equil_time = 100,
                                ...) {
    stopifnot(length(x.init) == length(u.init)) # make sure vectors are of equal length. Useful for other functions.
    stopifnot(length(x.init) == ncol(A))
    
    ntimesteps <- determine_ntimesteps(nsamples, spacing, equil_time, DeltaT)
    samples <- determine_samples(nsamples, spacing, ntimesteps, DeltaT)

    i <- 1
    Xs <- list() # will be a list of sample matrices (samples in rows, nodes in columns)
    Cs <- list() # will be a list of covariance matrices, computed from the Xs. 
    nlowerstate <- sum(x.init < r[2]) # for stopping condition

    u <- u.init
    while(nlowerstate == N) { # simulate until the first node leaves the lower basin of attraction
                                        # simulate the model, retaining only the samples from x_i(t)
        X <- doublewell(x.init, A, r, D, u, s, ntimesteps, DeltaT)[samples, ]
                                        # Recompute nlowerstate, for the stopping condition
        nlowerstate <- sum(X[nrow(X), ] < r[2])

        Xs[[i]] <- X
        Cs[[i]] <- cov(X)

                                        # Iterate the value of the control parameter and the list index
        if(bparam == "u") {
            if(i %% 10 == 0) print(max(u))
            u <- u + u.step
        } else if(bparam == "D") {
            if(i %% 10 == 0) print(D)
            D <- D + D.step
        }
        i <- i + 1 
    }

                                        # Return the x_i(t) samples and the covariance matrices
                                        # Discard the final simulation, in which one or more nodes
                                        # left the lower basin of attraction.
    return(list(Xs = Xs[-length(Xs)], Cs = Cs[-length(Cs)])) 
}

## Each of the remaining function pairs follows the same pattern as the doublewell()/simulate_doublewell() pair. 

SIS <- function(x.init, A, mu, D, s, ntimesteps, DeltaT) {
    N <- length(x.init)
    stopifnot(N == ncol(A))

    W <- preallocate_noise(s, N, ntimesteps)

    x <- x.init
    X <- matrix(0, nrow = ntimesteps, ncol = N)
    for(timestep in seq_len(ntimesteps)) {
        X[timestep, ] <- x
        x <- x +
                                        # The trick used for doublewell() doesn't work, because we
                                        # need both x_i and x_j. Instead, use outer() to make a matrix
                                        # with the same shape as A. The first argument is x_i, the
                                        # second is x_j. Then, use rowSums to sum over the x_js for
                                        # each x_i.
            (-mu*x + D*rowSums(A*outer(1 - x, x)))*DeltaT +
            W[timestep, ]*sqrt(DeltaT)
                                        # Set any negative values, which are unphysical, to zero.
        x <- ifelse(x < 0, 0, x)
    }

    return(X)
}

## This function (for the stopping condition below) allows for approaching the bifurcation point from either direction. Returns TRUE/FALSE for the while condition.
condition <- function(nlowerstate, N, direction) {
    switch(
        direction,
        fromlower = nlowerstate == N,
        fromupper = nlowerstate == 0
    )
}

## It only makes sense to simulate the SIS model with D (often called λ) as the control parameter, so the arguments and defaults reflect this. 
simulate_SIS <- function(A, bparam = "D", N = ncol(A),
                         x.init = rep(0.01, N),
                         s.adj = 0,
                         D = 0, D.step = 0.0025,
                         mu = 1, s = 0.0005 + s.adj,
                         nsamples = 100, spacing = 1, equil_time = 100,
                         direction = "fromlower",
                         ...) {
    stopifnot(bparam == "D")
    stopifnot(length(x.init) == ncol(A))
    
    N <- length(x.init)
    ntimesteps <- determine_ntimesteps(nsamples, spacing, equil_time, DeltaT)
    samples <- determine_samples(nsamples, spacing, ntimesteps, DeltaT)

    i <- 1
    Xs <- list()
    Cs <- list()
    nlowerstate <- sum(x.init < .1)

    while(condition(nlowerstate, N, direction)) {
        if(i %% 10 == 0) print(max(D))
        
        X <- SIS(x.init, A, mu, D, s, ntimesteps, DeltaT)[samples, ]
        nlowerstate <- sum(X[nrow(X), ] < .1)

        Xs[[i]] <- X
        Cs[[i]] <- cov(X)

        D <- D + D.step
        i <- i + 1
    }

    return(list(Xs = Xs[-length(Xs)], Cs = Cs[-length(Cs)]))
}

                                        # Gene regulatory/Michaelis-Menten
genereg <- function(x.init, A, params, D, u, s, ntimesteps, DeltaT) {
    N <- length(x.init)
    stopifnot(N == ncol(A))
    stopifnot(isSymmetric(unname(A))) # This function uses a calculation that requires a symmetric matrix

    B <- params["B"]
    f <- params["f"]
    h <- params["h"]

    W <- preallocate_noise(s, N, ntimesteps)

    x <- x.init
    X <- matrix(0, nrow = ntimesteps, ncol = N)
    for(timestep in seq_len(ntimesteps)) {
        X[timestep, ] <- x
        x <- x +
            (-B*(x^f) + D*colSums(A*((x^h)/(1 + (x^h)))) + u)*DeltaT + # only x_j, same as doublewell()
            W[timestep, ]*sqrt(DeltaT)
        x <- ifelse(x < 0, 0, x)
    }

    return(X)
}

simulate_genereg <- function(A, bparam = "u", N = ncol(A),
                             x.init = rep(5, N),
                             u.adj = rep(0, N), s.adj = 0,
                             ## Initial D is set inside the function, based on value of B (1 or 1/5)
                             u.init = rep(0, N) + u.adj, u.step = 0.01, D.step = 0.01,
                             params = c(B = 1),
                             s = 5e-6 + s.adj, # preferred
                             nsamples = 100, spacing = 1, equil_time = 100,
                             direction = "fromupper",
                             ...) {
    stopifnot(length(x.init) == length(u.init))
    stopifnot(length(x.init) == ncol(A))
    
    N <- length(x.init)
    params["f"] <- 1
    params["h"] <- 2
    ntimesteps <- determine_ntimesteps(nsamples, spacing, equil_time, DeltaT)
    samples <- determine_samples(nsamples, spacing, ntimesteps, DeltaT)

    i <- 1
    Xs <- list()
    Cs <- list()
    nlowerstate <- sum(x.init < .1)

    if(params["B"] == 1) D <- 1 else if(params["B"] == 1/5) D <- 0.35
    u <- u.init
    while(condition(nlowerstate, N, direction)) {
        X <- genereg(x.init, A, params, D, u, s, ntimesteps, DeltaT)[samples, ]
        nlowerstate <- sum(X[nrow(X), ] < .1)

        Xs[[i]] <- X
        Cs[[i]] <- cov(X)

        if(bparam == "u") {
            if(i %% 10 == 0) print(max(u))
            u <- u - u.step #note subtraction here. Decrementing bparam because starting from upper.
        } else if(bparam == "D") {
            if(i %% 10 == 0) print(D)
            D <- D - D.step
        }
        i <- i + 1
    }

    return(list(Xs = Xs[-length(Xs)], Cs = Cs[-length(Cs)]))
}

                                        # Mutualistic species
mutualistic <- function(x.init, A, params, D, u, s, ntimesteps, DeltaT) {
    N <- length(x.init)
    stopifnot(N == ncol(A))
    B <- params["B"]
    K <- params["K"]
    C <- params["C"]
    D. <- params["D."] # \tilde{D} in the paper
    E <- params["E"]
    H <- params["H"]

    W <- preallocate_noise(s, N, ntimesteps)

    x <- x.init
    X <- matrix(0, nrow = ntimesteps, ncol = N)
    for(timestep in seq_len(ntimesteps)) {
        X[timestep, ] <- x
        coupling <- rowSums(A*(outer(x, x)/(D. + outer(E*x, H*x, `+`)))) # rowSums b/c f(x_i, x_j)
        x <- x +
            (B + x*(1 - (x/K))*((x/C) - 1) + D*coupling + u)*DeltaT +
            W[timestep, ]*sqrt(DeltaT)
        x <- ifelse(x < 0, 0, x)
    }

    return(X)
}

simulate_mutualistic <- function(A, bparam = "u", N = ncol(A),
                                 x.init = rep(5, N),
                                 u.adj = rep(0, N), s.adj = 0,
                                 u.init = rep(0, N) + u.adj, u.step = 0.1,
                                 params = c(B = 0.1, K = 5, C = 1, D. = 5, E = 0.9, H = 0.1),
                                 D = 1, D.step = 0.01,
                                 s = 0.25 + s.adj,
                                        # note different time scale
                                 nsamples = 100, spacing = .1, equil_time = 10,
                                 ...) {
    stopifnot(length(x.init) == length(u.init))
    stopifnot(length(x.init) == ncol(A))
    
    N <- length(x.init)
    ntimesteps <- determine_ntimesteps(nsamples, spacing, equil_time, DeltaT)
    samples <- determine_samples(nsamples, spacing, ntimesteps, DeltaT)

    i <- 1
    Xs <- list()
    Cs <- list()
    nlowerstate <- sum(x.init < .1)

    u <- u.init
    while(nlowerstate == 0) {
        X <- mutualistic(x.init, A, params, D, u, s, ntimesteps, DeltaT)[samples, ]
        nlowerstate <- sum(X[nrow(X), ] < .1)

        Xs[[i]] <- X
        Cs[[i]] <- cov(X)

        if(bparam == "u") {
            if(i %% 10 == 0) print(min(u))
            u <- u - u.step
        } else if(bparam == "D") {
            if(i %% 10 == 0) print(D)
            D <- D - D.step
            if(D < 0) break
        }
        i <- i + 1
    }

    return(list(Xs = Xs[-length(Xs)], Cs = Cs[-length(Cs)]))
}
