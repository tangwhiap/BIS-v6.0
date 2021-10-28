#!/bin/bash

# Author: Wenhan TANG - 06/2021


Me=$( readlink -m $( type -p $0 ))
MyDir=`dirname $Me`
cd $MyDir

#-- Make sure that init.sh will be executed before running the program. --#
./init.sh

#-- Running BIS program in the background using command "nohup". --#
#-- Write the screen output information to run.log. --#

date=`date +%Y-%m-%d_%H:%M:%S`
LOGFILE=./logs/run_${date}.log

nohup ./BIS_run.py &> $LOGFILE &

#-- Screen output information. --#
echo ""
echo "=============================================="
echo "---       Bayesian Inversion System        ---"
echo "   The BIS is running in the background ...   "
echo "=============================================="
echo ""


