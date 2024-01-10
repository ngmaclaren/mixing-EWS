#!/bin/bash

#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=26
#SBATCH --mem=50000

#SBATCH --job-name=ewt-sims
#SBATCH --output=ewt-sims.out
#SBATCH --mail-user=neilmacl@buffalo.edu
#SBATCH --mail-type=ALL
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc

module load gcc
module load openmpi
module load r

declare -a networks=("ba" "erdosrenyi" "fitness" "chesapeake" "catlins" "dolphin")

declare -a dynamics=("doublewell" "mutualistic" "genereg") # no SIS for u

n=50

# outfile is <network>-<dynamics>-<b param>-<vary u?>-<vary σ?>

for i in "${networks[@]}"; do
    for j in "${dynamics[@]}"; do
	Rscript run-simulations.R --network=$i --dynamics=$j --ntrials=$n --bparam=u
	Rscript run-simulations.R --network=$i --dynamics=$j --ntrials=$n --bparam=u -v
	Rscript run-simulations.R --network=$i --dynamics=$j --ntrials=$n --bparam=u -v -s

	Rscript run-simulations.R --network=$i --dynamics=$j --ntrials=$n --bparam=D
	Rscript run-simulations.R --network=$i --dynamics=$j --ntrials=$n --bparam=D -v
	Rscript run-simulations.R --network=$i --dynamics=$j --ntrials=$n --bparam=D -v -s
    done
done

# deal with SIS separately
for i in "${networks[@]}"; do
    Rscript run-simulations.R --network=$i --dynamics=SIS --ntrials=$n --bparam=D
    Rscript run-simulations.R --network=$i --dynamics=SIS --ntrials=$n --bparam=D -s
done
