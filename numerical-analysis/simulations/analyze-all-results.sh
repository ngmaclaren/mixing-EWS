#!/bin/bash

#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=26
#SBATCH --mem=50000

#SBATCH --job-name=ewt-analysis
#SBATCH --output=ewt-analysis.out
#SBATCH --mail-user=neilmacl@buffalo.edu
#SBATCH --mail-type=ALL
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc

module load gcc
module load openmpi
module load r

declare -a networks=("ba" "erdosrenyi" "fitness" "chesapeake" "catlins" "dolphin")

declare -a dynamics=("doublewell" "mutualistic" "genereg")

n=50

for i in "${networks[@]}"; do
    for j in "${dynamics[@]}"; do
	Rscript reporting.R --network=$i --dynamics=$j --ntrials=$n --bparam=u
	Rscript reporting.R --network=$i --dynamics=$j --ntrials=$n --bparam=u -v
	Rscript reporting.R --network=$i --dynamics=$j --ntrials=$n --bparam=u -v -s

	Rscript reporting.R --network=$i --dynamics=$j --ntrials=$n --bparam=D
	Rscript reporting.R --network=$i --dynamics=$j --ntrials=$n --bparam=D -v
	Rscript reporting.R --network=$i --dynamics=$j --ntrials=$n --bparam=D -v -s
    done
done

for i in "${networks[@]}"; do
    Rscript reporting.R --network=$i --dynamics=SIS --ntrials=$n --bparam=D
    Rscript reporting.R --network=$i --dynamics=SIS --ntrials=$n --bparam=D -s
done
