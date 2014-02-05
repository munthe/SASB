# cfu_run_cmd_on_fcfu.bash --- 
# 
# Filename: cfu_run_cmd_on_fcfu.bash
# Description: 
# Author: Morten Fischer Rasmussen
# Maintainer: Morten Fischer Rasmussen (mofi)
# Created: Thur Oct 31 12:26:20 2013 (+0200)
# Version:  1.0
# Last-Updated: Thur Oct 31 12:26:20 2013 (+0200)
#           By: Morten Fischer
#     Update #: 1
# URL: 
# Keywords: 
# Compatibility: 
# 
# 

# Commentary: This script runs the input argument (command name) in a shell on each FCFU-node.
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

echo "Running the command \"$1\" in a shell on FCFU-nodes: "$START" to "$END"."
echo 'Running on: '
for ((i="$START";i<="$END";i+=1));
do  printf "fcfu%i\n"  $i $ID
    ssh  fcfu$i -f "$1";
    sleep 0.5s
done
printf  "\nDone.\n"
# 
# cfu_run_cmd_fcfu.bash ends here
