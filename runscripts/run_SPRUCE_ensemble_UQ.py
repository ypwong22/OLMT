import sys,os, time
sys.path.append('..')
import numpy as np
import subprocess
import pickle
import model_ELM
from optparse import OptionParser
import pandas as pd
import glob
from netCDF4 import Dataset
import matplotlib.pyplot as plt


#Python code used to manage the ensemble simulations 
#  and perform post-processing of model output.

parser = OptionParser()

caseid='UQ_20231118'
#caseid='FACE_r241231_CalibrationMCMC_RD'
#caseid2='FACE_r240107_CalibrationMCMCe_RD'
compset='ICB20TRCNPRDCTCBC'
suffix=''
site='US-SPR'
sitegroup='AmeriFlux'
runroot=os.path.join(os.environ['E3SM_ROOT'], 'output')
caseroot=os.path.join(os.environ['E3SM_ROOT'],'case_dirs')
if len(suffix) > 0:
  casename=caseid+'_'+site+'_'+compset +'_'+suffix
else:
  casename=caseid+'_'+site+'_'+compset

#Make sure to back up the old pkl files before this step!
#Create case object
mycase = model_ELM.ELMcase(caseid=caseid,compset=compset,site=site, \
        sitegroup=sitegroup, machine='cases-baseline', suffix=suffix, \
        runroot=runroot, caseroot=caseroot)
mycase.casename=casename

print(mycase.casename)

mycase.startyear=2015   #Starting year of run; doesn't have to be same as actual run
mycase.run_n=7         #Number of years for the run (to post process)
mycase.postproc_pfts=[]  #PFTs to postprocess
mycase.postproc_vars=[]
mycase.postproc_startyear=2015    #Starting year to postprocess/calibrate
mycase.postproc_endyear= 2021
mycase.postproc_freq = ''
mycase.read_parm_list(os.path.join(os.environ['HOME'], 'Git', 'phenology_elm', 
                                   'calibration_files', 'parm_file_20231118_forXiaoying'))
# is in OLMT/parm_files
samples_file=os.path.join(os.environ['HOME'],'models','OLMT_SPRUCE',
                          'mcsamples_UQ_20231118_4000.txt')
mycase.samples = (np.loadtxt(samples_file,)).transpose()
mycase.nsamples=4000
mycase.np_ensemble=mycase.samples.shape[0]
mycase.npernode=40 #128
mycase.obs={}
mycase.obs_err={}
mycase.OLMTdir=os.path.join(os.environ['HOME'],'models','OLMT')
mycase.pscaler={}
mycase.yscaler={}
mycase.rundir_UQ = mycase.runroot+'/UQ/'+mycase.casename

# -------------------------------------------------------------------------------------
# Observed values
# -------------------------------------------------------------------------------------
# Put into the dictionary as numpy arrays

# -------------------------------------------------------------------------------------
# custom outputs (change units, sum variables, etc)
# -------------------------------------------------------------------------------------
# We can do additional processing of the output time series here. 
## mycase.output['NPP_correct'] = (mycase.output['FATES_NPP']-mycase.output['FATES_EXCESS_RESP'])*24*3600*365*1000
## mycase.output['NUP'] = (mycase.output['FATES_NH4UPTAKE']+mycase.output['FATES_NO3UPTAKE'])*24*3600*365*1000
## mycase.postproc_vars.append('NPP_correct')
## mycase.postproc_vars.append('NUP')

VAR_COL = ['GPP', 'NPP', 'QVEGT', 'NEE', 'TOTVEGC']
VAR_PFT = ['GPP', 'NPP', 'QVEGT']
pft_list = [2, 3, 11, 12]
nvars = len(VAR_COL) + len(pft_list) * len(VAR_PFT)

collection_all = np.full([mycase.nsamples, nvars, 11], np.nan)
for part in range(7):
  temp = np.load(os.path.join(os.environ['PROJDIR'], 'ELM_Phenology', 'output', 'extract', 
                              caseid, f'xys_uncertainty_extraction_part{part:03g}.bin'),
                  allow_pickle=True)
  collection_all[(part*200):(part*200+200), :, :] = temp
mycase.output = {}
for v, var in enumerate(VAR_COL):
  # reshape to make sure the second dimension is ensemble member
  mycase.output[var] = np.mean(collection_all[:, v, :], axis = 1).reshape(1,-1)
for v, var in enumerate(VAR_PFT):
  for t, pft in enumerate(pft_list):
    # reshape to make sure the second dimension is ensemble member
    mycase.output[var+'_pft'+str(pft)] = np.mean(collection_all[:, 
                  len(VAR_COL)+v*len(pft_list)+t, :], axis = 1).reshape(1,-1)

#mycase.create_pkl(outdir=mycase.OLMTdir+'/pklfiles/')


#------UQ -----------------------------

#Train surrogate models

##mycase.postproc_vars.append('NPP_response')
##mycase.output['NPP_response'] = mycase2.output['NPP_correct'] - mycase.output['NPP_correct']
#Create a new variable "NPP_response" that is the differece between eCO2 (mycase2) and aCO2 (mycase1) NPP
##mycase.postproc_vars.append('NPP_response')
##mycase.output['NPP_response'] = mycase2.output['NPP_correct'] - mycase.output['NPP_correct']
#Add the data for this variable to be used in calibration
##mycase.obs['NPP_response'] = [74.51825141, \
##262.4690799, \
##260.8185289, \
##320.2143628, \
##334.8027339, \
##277.0862578, \
##267.9204289, \
##332.4450538, \
##349.5163926, \
##332.9064781, \
##355.4904338, \
##191.3458051]
##mycase.obs_err['NPP_response'] = [139.369053, \
##145.1765813, \
##104.0828343, \
##39.98227572, \
##67.97362259, \
##109.2344479, \
##163.6981652, \
##92.41676215, \
##141.423613, \
##201.1801581, \
##194.9031354, \
##196.5207146]

# break out the individual PFTs here by appending _{pft} to varname
variable_list = ['GPP','GPP_pft2','GPP_pft3','GPP_pft11','GPP_pft12',
                 'NPP','NPP_pft2','NPP_pft3','NPP_pft11','NPP_pft12',
                 'QVEGT','QVEGT_pft2','QVEGT_pft3','QVEGT_pft11','QVEGT_pft12',
                 'NEE','TOTVEGC']
mycase.train_surrogate(variable_list)

#run GSA
mycase.GSA(variable_list)


def plot_GSA_treatment():
    ticklabels = ['GPP','GPP picea','GPP larix','GPP shrub','GPP moss',
                  'NPP','NPP picea','NPP larix','NPP shrub','NPP moss',
                  'QVEGT','QVEGT picea','QVEGT larix','QVEGT shrub','QVEGT moss',
                  'NEE','TOTVEGC']

    x_pos = np.cumsum(np.ones(len(variable_list)))

    #Plot main sensitivity indices
    fig, axes = plt.subplots(2, 2, figsize = (14, 11))
    for i, (pft,pftname) in enumerate(zip([2, 3, 11, 0],
                                          ['picea','larix','shrub','Column'])):
      subset = np.where(np.array(mycase.ensemble_pfts) == pft)[0]

      ax = axes.flat[i]

      bottom = np.zeros(len(x_pos))
      for s in subset:
        temp = np.array([mycase.sens_main[v][s,0] for v in variable_list])
        ax.bar(x_pos, temp, align='center', # alpha=0.5,
               bottom = bottom, label = mycase.ensemble_parms[s])
        bottom = bottom + temp

      ax.set_ylim([0, 0.99])
      ax.set_xticks(x_pos)
      ax.set_xticklabels(ticklabels, rotation=90)
      ax.set_title(f'{pftname} parameters')
      ax.legend(loc = [1.05, 0.5])
    plt.tight_layout()
    plt.savefig(f'sens_main_{caseid}.png')

    #Total sensitivity indices
    fig, axes = plt.subplots(2, 2, figsize = (14, 11))
    for i, (pft,pftname) in enumerate(zip([2, 3, 11, 0],
                                          ['Spruce','Larch','Shrub','Column'])):
      subset = np.where(np.array(mycase.ensemble_pfts) == pft)[0]

      ax = axes.flat[i]

      bottom = np.zeros(len(x_pos))
      for s in subset:
        temp = np.array([mycase.sens_tot[v][s,0] for v in variable_list])
        ax.bar(x_pos, temp, align='center', # alpha=0.5,
               bottom = bottom, label = mycase.ensemble_parms[s])
        bottom = bottom + temp

      ax.set_ylim([0, 0.99])
      ax.set_xticks(x_pos)
      ax.set_xticklabels(ticklabels, rotation=90)
      ax.set_title(f'{pftname} parameters')
      ax.legend(loc = [1.05, 0.5])
    plt.tight_layout()
    plt.savefig(f'sens_tot_{caseid}.png')


#plot GSA
plot_GSA_treatment()

###Save postprocessed output
##mycase.create_pkl(outdir=mycase.OLMTdir+'/pklfiles/')

##Set intial values for parameters
##parms=((np.array(mycase.ensemble_pmax)+np.array(mycase.ensemble_pmin))/2)

##Run MCMC for the 2 varibles of interest
##mycase.MCMC(parms, ['TLAI_ann','QRUNOFF_ann'], 100000)

