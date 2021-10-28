#!/bin/bash

# Authors:
#   Wenhan TANG - 06/2021
#   ...

# This script is used to create the root directory of the program and
# write it to the configuration file "configure".

# It must be executed before running the program
# if the directory of it has been changed.

# This script must be located in the root directory of the program.

#-- Get directory of this script --#
Me=$( readlink -m $( type -p $0 ))
MyDir=`dirname $Me`
cd $MyDir

#-- Get the root directory where the script in here --#
ROOTDIR=$MyDir

#-- Write the directory in to configure.py --#
sed -i "/ROOTDIR = /cROOTDIR = \"${ROOTDIR}\" # by init.sh" $ROOTDIR/Bayesian/main/configure.py
