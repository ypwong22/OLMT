#!/bin/bash
#SBATCH --time=24:0:00
#SBATCH -J UQ_extract
#SBATCH --nodes=1
#SBATCH -A CLI185
#SBATCH -p batch_ccsi
#SBATCH --ntasks-per-node 128

cd ${HOME}/models/OLMT/runscripts
srun -n 128 ${HOME}/.conda/envs/olmt/bin/python -u run_SPRUCE_MIP_UQ_extraction.py
