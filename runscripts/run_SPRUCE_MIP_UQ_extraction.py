#!/usr/bin/env python3

# Default script to extract data is too slow
# Use this one

import sys, os, time
import numpy as np
from scipy.stats import linregress, t
from mpi4py import MPI
from netCDF4 import Dataset
import pandas as pd
import time

# delete when debug
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

workdir = os.getcwd()

# number of simulations
N = 1920

PREFIX = '20240723'

path_out = os.path.join(os.environ['PROJDIR'], 'ELM_Phenology', 'output', 'extract', PREFIX)
os.makedirs(path_out, exist_ok=True)

# number of ensembles to save in each bin file
# this avoids having difficulty in dumping file
BLOCK = 160

if np.mod(N, BLOCK) != 0:
    raise Exception("N must be a multiply of BLOCK")

RUNROOT = "/gpfs/wolf2/cades/cli185/proj-shared/zdr/SPRUCE/e3sm_run/UQ"
niter = int(N / BLOCK)


PLOT_LIST = ['TAMB','T0.00','T0.00eCO2','T2.25','T2.25eCO2','T4.50','T4.50eCO2',
             'T6.75','T6.75eCO2','T9.00','T9.00eCO2']

# variables for Xiaoying Shi
VAR_COL = ['NEE','NPP','CH4PROD','FCH4','HR']
VAR_PFT = ['NPP', 'TLAI', 'TOTVEGC']

pft_list = [2, 3, 11, 12]

nvars = len(VAR_COL) + len(pft_list) * len(VAR_PFT)

# Function to perform post-processing for one ensemble member
def postproc(thisjob, collection):

    for pind, plot in enumerate(PLOT_LIST):

        casename = f"{PREFIX}_US-SPR_ICB20TRCNPRDCTCBC_{plot}"
        baserundir = os.path.join(RUNROOT, casename, f"g{thisjob:05g}")

        ##DEBUG
        ##print(casename, thisjob, plot)

        # extract column variables
        flist_col = [
             os.path.join(baserundir, f'{casename}.elm.h1.{year}-01-01-00000.nc')
             for year in range(2015, 2022)
        ]
        temp = np.full([len(VAR_COL), len(flist_col)], np.nan)
        for f,file in enumerate(flist_col):
            nc = Dataset(file, 'r')
            for v, var in enumerate(VAR_COL):
                temp[v,f] = np.mean(nc[var][:, 0]) * 0.64 + np.mean(nc[var][:, 1]) * 0.36
            nc.close()
        collection[:len(VAR_COL), pind] = np.nanmean(temp, axis = 1)

        # extract pft variables
        flist_pft = [
             os.path.join(baserundir, f'{casename}.elm.h2.{year}-01-01-00000.nc')
             for year in range(2015, 2022)
        ]
        temp = np.full([len(VAR_PFT)*len(pft_list), len(flist_pft)], np.nan)
        for f,file in enumerate(flist_pft):
            nc = Dataset(file, 'r')
            for v, var in enumerate(VAR_PFT):
                for t, pft in enumerate(pft_list):
                    temp[v*len(pft_list) + t,f] = np.mean(nc[var][:, pft]) * 0.64 + \
                                                  np.mean(nc[var][:, pft+32]) * 0.36
            nc.close()
        collection[len(VAR_COL):, pind] = np.nanmean(temp, axis = 1)


## Debug
#thisjob = 3852
#collection = np.empty([nvars, len(PLOT_LIST)])
#postproc(thisjob, collection)


for b in range(niter): # range(niter):
    print("rank = ", rank, "b = ", b, flush = True)

    if rank == 0:
        # process a whole block
        start = b * BLOCK
        collection_all = np.full([BLOCK, nvars, len(PLOT_LIST)], np.nan)

        # --------------------------Collect and save data---------------------
        n_done = 0

        # send first np-1 jobs where np is number of processes
        for n_job in range(1, size):
            comm.send(n_job, dest=n_job, tag=1)
            comm.send(0, dest=n_job, tag=2)
            comm.send(start, dest=n_job, tag=3)

        # assign rest of jobs on demand
        for n_job in range(size, BLOCK + 1):
            process = comm.recv(source=MPI.ANY_SOURCE, tag=4)
            thisjob = comm.recv(source=process, tag=5)
            collection = comm.recv(source=process, tag=6)
            collection_all[thisjob - 1, :, :] = collection
            n_done = n_done + 1
            comm.send(n_job, dest=process, tag=1)
            comm.send(0, dest=process, tag=2)
            comm.send(start, dest=process, tag=3)

        # receive remaining messages and finalize
        while n_done < BLOCK:
            process = comm.recv(source=MPI.ANY_SOURCE, tag=4)
            thisjob = comm.recv(source=process, tag=5)
            collection = comm.recv(source=process, tag=6)
            collection_all[thisjob - 1, :, :] = collection
            n_done = n_done + 1
            comm.send(-1, dest=process, tag=1)
            comm.send(-1, dest=process, tag=2)
            comm.send(-1, dest=process, tag=3)

        #collection_all.dump(
        #    os.path.join(path_out, f"xys_uncertainty_extraction_part{b:03g}.bin")
        #)
        collection_all.dump(
            os.path.join(path_out, f"uncertainty_extraction_part{b:03g}.bin")
        )

    # --------------------- Slave process (get data and calculate) --------------
    else:
        # process individual members inside the block
        collection = np.full([nvars, len(PLOT_LIST)], np.nan)

        status = 0
        while status == 0:
            myjob = comm.recv(source=0, tag=1)
            status = comm.recv(source=0, tag=2)
            start = comm.recv(source=0, tag=3)

            if status == 0:
                ierr = postproc(start + myjob, collection)
                comm.send(rank, dest=0, tag=4)
                comm.send(myjob, dest=0, tag=5)
                comm.send(collection, dest=0, tag=6)

MPI.Finalize()
