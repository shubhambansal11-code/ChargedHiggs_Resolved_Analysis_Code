#!/bin/bash

#Submit script for gbb Tuple Ana
echo "Submit jobs for Charged Higgs Analysis......"

#OUTPATH="/afs/cern.ch/work/s/shbansal/chargedHiggs_Ana/MockOutput"
LOG_FOLDER="/cephfs/user/s6subans/ChargedHiggsLog/"

echo "Logs in: ${LOG_FOLDER}"
echo "Output in: ${OUTPATH}"

OUTDIR="/cephfs/user/s6subans/ChargedHiggs_Mock/"
#USE_BATCH_MODE=1

WP=(
#"77p"
"70p"
#"85p"
#"60p"
)

export EOS_MGM_URL=root://eosuser.cern.ch
echo "reading files !!! "
Paths=(
#"ChargedHiggsNtups/1Lep/Data/"
"ChargedHiggsNtups/1Lep/MC16a/"
"ChargedHiggsNtups/1Lep/MC16d/"
"ChargedHiggsNtups/1Lep/MC16e/"
#"/afs/cern.ch/work/s/shbansal/chargedHiggs_Ana/chargedHiggs_MC16a/"
)

File=(
#data15.root
#data16.root
#data17.root 
#data18.root   
diboson.root
sig_Hplus_Wh_m400-0.root
sig_Hplus_Wh_m800-0.root
sig_Hplus_Wh_m1600-0.root
ttbar.root
Wjets.root
singleTop.root
Zjets.root    
)

rm -rf cluster_pack.tar.gz
tar -czf cluster_pack.tar.gz Makefile_batch main_batch.C  main_RunMVATraining.C dataset/ main/ TH1Fs/ utilis/ LatexOutput/ python/ style/

echo "sucessfuly opened tar files"
for wp in "${WP[@]}"
do
   for path in "${Paths[@]}"
   do
      for file in "${File[@]}"
      do
      ##./execute $path $file $wp $OUTDIR $USE_BATCH_MODE
      condor_submit IN_PATH="${path}" FILE="${file}" WP="${wp}" OUTDIR="${OUTDIR}" /cephfs/user/s6subans/ChargedHiggsAna/Code/run_Hplus.sub
      done
   done
done
echo "all done !!! "