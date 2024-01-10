#!/bin/bash

# The primary results are for k1=0.1 k2=0.9 and are already available
# the rest of the combinations are robustness checks
k1s=(0.1 0.1 0.3 0.1 0.3 0.5 0.1 0.3 0.5 0.7)
k2s=(0.9 0.7 0.9 0.5 0.7 0.9 0.3 0.5 0.7 0.9)

for i in ${!k1s[@]}; do
    k1=${k1s[i]}
    k2=${k2s[i]}
    sbatch overall-summary-fig.sh ${k1} ${k2}
    # for debugging:
    # echo ${k1} ${k2}
    # sh overall-summary-fig.sh ${k1} ${k2}
done
