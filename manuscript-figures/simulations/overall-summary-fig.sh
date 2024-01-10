#!/bin/bash

#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=26
#SBATCH --mem=50000

#SBATCH --job-name=ewt-osumfig
#SBATCH --output=ewt-osumfig.out
#SBATCH --mail-user=neilmacl@buffalo.edu
#SBATCH --mail-type=ALL
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc

module load gcc
module load openmpi
module load r

k1=${1}
k2=${2}

Rscript overall-summary-fig.R --bparam=u --makedata ${k1} ${k2}
Rscript overall-summary-fig.R --bparam=D --makedata ${k1} ${k2}

Rscript overall-summary-fig.R ${k1} ${k2}
