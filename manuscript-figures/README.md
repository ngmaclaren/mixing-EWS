This folder contains the data and code to reproduce the figures in the manuscript.

## Manuscript figures

The `./figure*` directories contain the data and code to reproduce the named figures in the manuscript. Note that these are "raw" figures (i.e., without formatting for consistent style, etc.). Each directory is self contained: the `.py` or `.R` files in each directory use the data files located in the same directory to generate the corresponding figure.

Python code in the `./figure*` directories depends on `numpy`, `matplotlib`, and `scipy`. R code depends on `sfsmisc`, `latex2exp`, and `optparse`, as well as `./manuscript-figures/plot-functions.R`. Shell script files (e.g., `./figure8_S2--S28_S31/make-plots.sh`) are written in Bash.

## Simulating from scratch

The `./simulations` directory contains the code required to generate the simulations we used in Section 4, "Numerical results for larger networks and various dynamical systems", of the manuscript. Details for running those simulations can be found in that directory's README.
