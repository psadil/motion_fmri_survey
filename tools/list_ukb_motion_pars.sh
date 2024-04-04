#!/bin/bash

# list 

OUT=/fastscratch/myscratch/pssadil/ukb_mcf.par
ROOT=/dcs07/smart/data/ukb/rawdata
mc="ses-2/non-bids/20249/fMRI/tfMRI.feat/mc"

for s in ${ROOT}/*; do 
  if [[ -d ${s}/${mc} ]]; then 
    echo ${s}/${mc}/prefiltered_func_data_mcf.par >> $OUT
  fi
done

# add column with index
# helpful for keeping track of trs
for f in data/ukb/*/ses-2/non-bids/20249/fMRI/tfMRI.feat/mc/prefiltered_func_data_mcf.par; do 
  awk '{printf "%s\t%s\n",$0,NR}' ${f} > tmp && mv -f tmp ${f}
done

for f in ukb/*/ses-2/non-bids/20249/fMRI/tfMRI.feat/mc/prefiltered_func_data_mcf.par; do 
  lines=$(cat $f | wc -l)
  if (( $lines > 332)); then
    echo ${f}
  fi
done

