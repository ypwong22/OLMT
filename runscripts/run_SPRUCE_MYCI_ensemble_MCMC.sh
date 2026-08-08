#!/bin/bash
#SBATCH -p serial
#SBATCH -q normal
#SBATCH --mem=64gb
#SBATCH --time=1-00:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 64
#SBATCH -J MCMC
#SBATCH -o out/%x-%J.out
#SBATCH -e out/%x-%J.err

cd ${HOME}/models/OLMT/runscripts
~/.conda/envs/OLMT_pf/bin/python run_SPRUCE_MYCI_ensemble_MCMC.py