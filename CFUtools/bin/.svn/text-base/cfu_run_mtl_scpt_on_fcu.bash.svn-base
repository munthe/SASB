# run_on_fcfu.sh --- 
# 
# Filename: run_on_fcfu.sh
# Description: 
# Author: Morten Fischer Rasmussen
# Maintainer: Morten Fischer Rasmussen (mofi)
# Created: Mon Apr  9 21:28:24 2012 (+0200)
# Version:  1.0
# Last-Updated: Tue Apr 10 09:10:56 2012 (+0200)
#           By: Morten Fischer
#     Update #: 1
# URL: 
# Keywords: 
# Compatibility: 
# 
# 

# Commentary: This script runs the input argument (file name) in MATLAB on each FCFU-node.
# 
# 
# 
# 

# Change Log:
# 
# 
# 

# Code:


#!/bin/bash

CURRENT_DIR=$(pwd);


printf "\n"

# test that first input is a file
if [ ! -f $1 ]
then
    echo "Error: first input must be a file."
    exit 1
fi

if [ ! -d 'logs' ]
then
    printf "Creating directory: %s/logs\n" $CURRENT_DIR
    mkdir logs;
fi


# Test for 2nd and 3rd argument
if [[ -n "$2" ]]
then
	START="$2"
else
	START=1
fi

if [[ -n "$3" ]]
then
	END="$3"
else
	END=12
fi



echo "Running file \"$1\" in MATLAB on FCFU-nodes: "$START" to "$END"."
echo "The logs will be put in the logs sub-directory.";

echo 'Starting MATLAB on: '
for ((i="$START";i<="$END";i+=1));
do  ID=$RANDOM   # get a random ID
    printf "fcfu%i  ID=%i\n"  $i $ID
    ssh  fcfu$i -f "cd $CURRENT_DIR ; nice -n 19 matlab_cluster -r $1 \
> logs/"$1"_fcfu"$i"_"$ID".log \
2> logs/"$1"_fcfu"$i"_"$ID".err ";
    sleep 0.5s
done

printf  "\nDone.\n"
# 
# run_on_fcfu.sh ends here
