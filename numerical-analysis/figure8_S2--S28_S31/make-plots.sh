#!/bin/bash

# The primary results (i.e., Figures 8 and S31) are for k1=0.1 and k2=0.9
# the rest of the combinations are robustness checks and make Figures S2--S28
k1s=(0.1 0.1 0.3 0.1 0.3 0.5 0.1 0.3 0.5 0.7)
k2s=(0.9 0.7 0.9 0.5 0.7 0.9 0.3 0.5 0.7 0.9)

for i in ${!k1s[@]}; do
    k1=${k1s[i]}
    k2=${k2s[i]}
    Rscript plot.R ${k1} ${k2}
done
