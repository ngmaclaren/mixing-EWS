## File flow

The source files are calc-functions2.R and sim-functions2.R. I would like to change the names of these files to calc-functions.R and sim-functions.R, and write a plot-functions.R file that collects the plotting functions in one place (and takes them out of calc-functions.R). A demonstration of the functions discussed in the paper can be found in mwe.R (at the top level of this repo). 

With the above files in the working directory, the analyses in the manuscript proceed as follows:

- sbatch run-all-simulations.sh (calls run-simulations.R); must be complete before running anything else
- sbatch analyze-all-results.sh (calls reporting.R)
- sbatch summary-fig.sh (calls summary-fig.R)
- sh run-all-osumfigs.sh, which will call sbatch overall-summary-fig.sh for the appropriate k1/k2 combinations. 
- sbatch tau-tau-scatter-data.sh
- sh run-remaining-figs.sh: tau-tau-scatter.R, tau-tau-scatter-highlowinput.R, overall-summary-fig-alt.R

The analysis files as written assume they will run on a slurm-based high-performance compute cluster. It is likely the .sh files will not run without modification for your particular environment. At the very least, directories and file names may need to be changed. However, the .R files can be used interactively either by setting variables directly (e.g., args$network <- "ba") or by setting the shell args using optparse's tools (e.g., in the file, set args <- parse_args(..., args = c("--network=ba"), ...) or call the .R file from the commandline with e.g. Rscript reporting.R --network=ba). Be sure to set the number of trials and number of cores to reasonable values.
