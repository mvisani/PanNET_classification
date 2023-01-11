#!/usr/bin/env bash

#creates a random forest model based on the training data and a number of probes
# specified by the user. All results of the pipeline go in the "out" directory 
# 
if [[ -z "${1}" ]]; then
  echo "Please enter the top n probes you want to select" 
  exit 1
fi

Rscript 03_training.R out_${1}_probes ${1};
Rscript 04_cross_validation.R out_${1}_probes ${1};
Rscript 05_calibration.R out_${1}_probes ${1}