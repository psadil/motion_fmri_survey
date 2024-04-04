#!/bin/bash

# list 

OUT=/fastscratch/myscratch/pssadil/hcp_movement_regressors
ROOT=/dcs07/smart/data/human-connectome-project-openaccess/HCP1200

for s in ${ROOT}/*; do
  for ped in LR RL; do
    for task in ${s}/MNINonLinear/Results/tfMRI*${ped}; do
      regs=${task}/Movement_Regressors.txt
      if [[ -f $regs  ]]; then
        echo $regs >> $OUT
      fi
    done
  done
done

EV=/fastscratch/myscratch/pssadil/hcp_evs

for s in ${ROOT}/*; do
  for ped in LR RL; do
    for task in ${s}/MNINonLinear/Results/tfMRI*${ped}; do
      evs=${task}/EVs
      if [[ -d $evs  ]]; then
        echo $evs >> $EV
      fi
    done
  done
done


# add column with index
# helpful for keeping track of trs
for f in hcp/*/MNINonLinear/Results/*/Movement_Regressors.txt; do 
  cut -w -f 2,3,4,5,6,7,8,9,10,11,12,13 ${f} | cat -b > tmp && mv -f tmp ${f}
done

# same for ev files
for f in hcp_evs/*/MNINonLinear/Results/*/EVs/*; do 
  cat -b ${f} > ev_tmp && mv -f ev_tmp ${f}
done

