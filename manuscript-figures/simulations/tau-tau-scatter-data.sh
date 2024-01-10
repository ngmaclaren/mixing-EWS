#!/bin/bash

#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=26
#SBATCH --mem=50000

#SBATCH --job-name=tautau
#SBATCH --output=tautau.out
#SBATCH --mail-user=neilmacl@buffalo.edu
#SBATCH --mail-type=ALL
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc

module load gcc
module load openmpi
module load r

Rscript tau-tau-scatter-data.R
Rscript tau-tau-scatter-data-highlowinput.R
