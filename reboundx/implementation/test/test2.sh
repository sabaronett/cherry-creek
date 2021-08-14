#!/usr/bin/bash

cd $PBS_O_WORKDIR
/home/barons2/.conda/envs/rebx-3.4.1/bin/python args.py $1
