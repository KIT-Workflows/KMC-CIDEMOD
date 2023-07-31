#!/bin/bash

export NANOVER=V4
source $NANOMATCH/$NANOVER/configs/simstack_kit.config
# source the conda environment
source /home/ws/ab5528/miniconda3/etc/profile.d/conda.sh
conda activate cidemod

python3 CIDEMOD.py
