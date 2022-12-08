#!/usr/bin/env bash

if [[ -z "${1}" ]]; then
  echo "Please enter the top n probes you want to select" 
  exit 1
fi

Rscript 03_training.R out_${1}_probes ${1};
Rscript 04_cross_validation.R out_${1}_probes ${1};
Rscript 05_calibration.R out_${1}_probes ${1}