#!/bin/sh

#cd /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.18.04/x86_64-centos7-gcc48-opt/bin/
#source thisroot.sh
#cd -

#source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.24.00/x86_64-centos7-gcc48-opt/bin/thisroot.sh

source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.36.00/x86_64-almalinux9.5-gcc115-opt/bin/thisroot.sh

export LD_LIBRARY_PATH=./lib:DynamicTTree/lib/:CfgManager/lib/:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=./lib:DynamicTTree/lib/:CfgManager/lib/:$DYLD_LIBRARY_PATH
export ROOT_INCLUDE_PATH=./interface:DynamicTTree/interface/:CfgManager/interface/:$ROOT_INCLUDE_PATH
