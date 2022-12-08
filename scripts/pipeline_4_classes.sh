#!/usr/bin/env bash

if [[ -z "${1}" ]]; then
  echo "Please enter the top n probes you want to select" 
  exit 1
fi

Rscript 06_training_four_classes.R out_${1}_probes_4_classes ${1};
Rscript 07_CV_four_classes.R out_${1}_probes_4_classes ${1};
Rscript 08_calibration_four_classes.R out_${1}_probes_4_classes ${1}