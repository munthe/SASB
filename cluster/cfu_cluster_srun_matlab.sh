#!/bin/bash
# cfu_cluster_srun_matlab.sh
# Run a specified matlab function with a specified set of arguments
# Automatically check for number/non-number type of parameters
# Argument 0 is the name of this script
# Argument 1 is the name of the matlab script
# Arguments 2 and onwards are parameters for MATLAB. Even numbered
# arguments hold the name of a variable that will be present in the
# MATLAB workspace, while the following odd numbered argument holds
# the value of that variable.

hostname
date
args="'$1'"
for ((i=2;i<=$#;i++))
do
    j=$(expr $i + 1)
    args="$args,'${!i}'"
    if ! [[ "${!j}" =~ ^-?[0-9]+([.][0-9]+)?([eE]-?[0-9]+)?$ ]] # Match a possibly negative, possibly decimal number, possible scientific notation
    then
	args="$args,'${!j}'"
    else
	args="$args,${!j}"
    fi
    i=j
done
srun nice -n19 matlab -nodisplay -nosplash -r "cfu_cluster_run_matlab($args)"
