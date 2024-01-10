#!/bin/bash

#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=26
#SBATCH --mem=50000

#SBATCH --job-name=summaryfig
#SBATCH --output=summaryfig.out
#SBATCH --mail-user=neilmacl@buffalo.edu
#SBATCH --mail-type=ALL
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc

module load gcc
module load openmpi
module load r

declare -a networks=("ba" "erdosrenyi" "fitness" "chesapeake" "catlins" "dolphin")
declare -a dynamics=("doublewell" "SIS" "mutualistic" "genereg")

for i in "${networks[@]}"; do
    for j in "${dynamics[@]}"; do
	Rscript summary-fig.R --network=$i --dynamics=$j
    done
done
