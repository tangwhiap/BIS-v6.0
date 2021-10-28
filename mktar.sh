#!/bin/bash

# Authors:
#   Wenhan TANG - 06/2021
#   ...

#-- Make sure the current directory is the ROOTDIR of the program. --#
Me=$( readlink -m $( type -p $0 ))
MyDir=`dirname $Me`
RootDir=$MyDir
cd $MyDir

#-- Where to save the tar file. --#
TarDir=/home/tangwh/packages
#-- Prefix of the tar file name. --#
Prefix=BIS_v6.0
#-- Author name. --#
Author=TWH
#-- The date (time) when the tar file created. --#
Tnow=`date +%Y-%m-%d_%H-%M`
#-- Suffix of the tar file name. --#
Suffix=TAR.gz

#-- The name of this program's main folder. --#
RootName=BIS_v6.0

#-- Format of tar file name <Prefix>.<Author>.<Date&Time>.<Suffix>
TarName=${Prefix}.${Author}.${Tnow}.${Suffix}

#-- Go to the parrent directory. --#
cd ..

#-- Create the tar file. --#
tar -zcvf $TarDir/$TarName $RootName

