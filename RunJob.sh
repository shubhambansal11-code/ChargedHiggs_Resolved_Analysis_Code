#!/bin/bash

HOMEDIR=/cephfs/user/s6subans/ChargedHiggsAna/Code

cd ${HOMEDIR}

source /etc/profile
#$BUDDY/.bashrc_CentOS7
#setupATLAS 
#lsetup "root 6.14.04-x86_64-slc6-gcc73-opt"
#export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
#. ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh
setupATLAS
#module load root/6.14.04
#lsetup "root 6.20.06-x86_64-centos7-gcc8-opt"
lsetup  "root 6.18.04-x86_64-centos7-gcc8-opt"
#module load root/6.18.04

IN_PATH=${1}
FILE=${2}
WP=${3}
OUTDIR=${4}
#USEBATCHMODE=${5}
Cluster=${5}
Process=${6}

## Run in $TMPDIR for fast I/O - it will be something like /tmp/7993986.1.short, based on the task ID and the queue name
#testdir=$TMPDIR/run_dir

#if [ -d $testdir ]
#then
    #rm -r $testdir
#fi

#mkdir $testdir
#cd $testdir



echo "${OUTDIR}"
echo "${IN_PATH}"
echo "${FILE}"
echo "${WP}"
#echo $5

./execute.exe ${IN_PATH} ${FILE} ${WP} ${OUTDIR}
#./execute /eos/user/s/shbansal/ChargedHiggsAna/Charged_Higgs_Ntups/MC16a/ sig_Hplus_Wh_m400-0.root 70p /afs/cern.ch/work/s/shbansal/chargedHiggs_Ana/MockOutput 0



#cp *.root /afs/cern.ch/work/s/shbansal/chargedHiggs_Ana/MockOutput/PlotFiles/
#cd /tmp/ 
#cp slurm.%j.out /eos/user/s/shbansal/chargedHiggs_Ana/MockOutput/LogFiles/
#OUTDIR = /cephfs/user/s6subans/ChargedHiggs_Mock/
echo "Write Output File to: ${OUTDIR} ....."

#root -h
#mkdir -p ${OUTDIR}

#cd ${OUTDIR}
#ls
#cp -rv *.root* ${OUTDIR}