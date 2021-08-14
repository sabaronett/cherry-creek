#!/usr/bin/bash

cd $PBS_O_WORKDIR
module load conda
conda activate rebx-3.4.1
python args.py $1
