#!/bin/bash
#SBATCH --time=24:0:00
#SBATCH -J MCMC
#SBATCH --nodes=1
#SBATCH -A CLI185
#SBATCH -p batch_ccsi
#SBATCH --ntasks-per-node 1

cd ${HOME}/models/OLMT/runscripts
python run_ERW_HBR_LAI_ensemble_MCMC.py