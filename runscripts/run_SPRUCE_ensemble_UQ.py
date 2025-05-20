import sys,os
sys.path.append('..')
import numpy as np
import pickle
import model_ELM
from optparse import OptionParser
import matplotlib.pyplot as plt
from string import ascii_lowercase


#Python code used to manage the ensemble simulations 
#  and perform post-processing of model output.

parser = OptionParser()

#caseid='UQ_20231118'
caseid='UQ_20240107'
#caseid='FACE_r241231_CalibrationMCMC_RD'
#caseid2='FACE_r240107_CalibrationMCMCe_RD'
compset='ICB20TRCNPRDCTCBC'
suffix=''
site='US-SPR'
sitegroup='AmeriFlux'
runroot=os.path.join(os.environ['E3SM_ROOT'], 'output')
caseroot=os.path.join(os.environ['E3SM_ROOT'],'case_dirs')
retrain_model = False # True - re-generate the "output" dict and overwrite pklfile; 
                      # False - directly read from pklfile
if len(suffix) > 0:
  casename=caseid+'_'+site+'_'+compset +'_'+suffix
else:
  casename=caseid+'_'+site+'_'+compset

temp = caseid.replace('UQ_', '')
parmfile = f'parm_file_{temp}_compact'

VAR_COL = ['GPP', 'NEE', 'HR', 'TOTVEGC', 'TOTSOMC']
VAR_PFT = ['GPP', 'AR', 'MR', 'GR', 'XR']
# variables for Xiaoying Shi
##VAR_COL = ['GPP', 'NPP', 'QVEGT', 'NEE', 'TOTVEGC']
##VAR_PFT = ['GPP', 'NPP', 'QVEGT']
pft_list = [2, 3, 11, 12]
nvars = len(VAR_COL) + len(pft_list) * len(VAR_PFT)


# break out the individual PFTs here by appending _{pft} to varname
variable_list = []
for var in VAR_COL:
    variable_list.append(var)
    if var in VAR_PFT:
        variable_list.extend([var+'_pft'+str(pft) for pft in pft_list])
for var in VAR_PFT:
    if not var in VAR_COL:
        variable_list.extend([var+'_pft'+str(pft) for pft in pft_list])


#Make sure to back up the old pkl files before this step!
#Create case object
if retrain_model:
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
                                    'calibration_files', parmfile))
  # is in OLMT/parm_files
  samples_file=os.path.join(os.environ['HOME'],'models','OLMT_SPRUCE',
                            f'mcsamples_{caseid}_4000.txt')
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

  BLOCK = 200
  collection_all = np.full([mycase.nsamples, nvars, 11], np.nan)
  for part in range(mycase.nsamples // BLOCK):
    temp = np.load(os.path.join(os.environ['PROJDIR'], 'ELM_Phenology', 'output', 'extract',
                                caseid, f'uncertainty_extraction_part{part:03g}.bin'), # f'xys_uncertainty_extraction_part{part:03g}.bin'),
                    allow_pickle=True)
    collection_all[(part*BLOCK):(part*BLOCK+BLOCK), :, :] = temp
  mycase.output = {}
  for v, var in enumerate(VAR_COL):
    # reshape to make sure the second dimension is ensemble member
    mycase.output[var] = np.mean(collection_all[:, v, :], axis = 1).reshape(1,-1)
  for v, var in enumerate(VAR_PFT):
    for t, pft in enumerate(pft_list):
      # reshape to make sure the second dimension is ensemble member
      mycase.output[var+'_pft'+str(pft)] = np.mean(collection_all[:, 
                    len(VAR_COL)+v*len(pft_list)+t, :], axis = 1).reshape(1,-1)

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

  #mycase.create_pkl(outdir=mycase.OLMTdir+'/pklfiles/')

  mycase.train_surrogate(variable_list)

  #run GSA
  mycase.GSA(variable_list)

  #Save postprocessed output
  mycase.create_pkl(outdir=mycase.OLMTdir+'/pklfiles/')

else:

  myfile=open(os.path.join(os.environ['HOME'],'models','OLMT','pklfiles', casename+'.pkl'),'rb')
  mycase=pickle.load(myfile)
  myfile.close()


def plot_GSA_treatment():
    pft_names = ['Spruce','Tamarack','Shrub','Moss']
    ticklabels = []
    for var in VAR_COL:
        ticklabels.append(var)
        if var in VAR_PFT:
            ticklabels.extend([var+' '+pname for pname in pft_names])
    for var in VAR_PFT:
        if not var in VAR_COL:
            ticklabels.extend([var+' '+pname for pname in pft_names])

    x_pos = np.cumsum(np.ones(len(variable_list)))

    #Plot main sensitivity indices
    fig, axes = plt.subplots(2, 2, figsize = (14, 11), sharex = True, sharey = True)
    for i, (pft,pftname) in enumerate(zip([2, 3, 11, 0], pft_names[:-1] + ['Column'])):
      subset = np.where(np.array(mycase.ensemble_pfts) == pft)[0]

      ax = axes.flat[i]

      bottom = np.zeros(len(x_pos))
      for s in subset:
        temp = np.array([mycase.sens_main[v][s,0] for v in variable_list])
        ax.bar(x_pos, temp, align='center', # alpha=0.5,
               bottom = bottom, label = mycase.ensemble_parms[s])
        bottom = bottom + temp

      # add a line for total
      total = np.array([mycase.sens_main[v][:,0].sum() for v in variable_list])
      ax.plot(x_pos, total, '-k', label = 'Total')

      ax.set_xticks(x_pos)
      ax.set_xticklabels(ticklabels, rotation=90)
      ax.set_title(f'{pftname} parameters')
      ax.legend(loc = [1.05, 0.5])
    for ax, lab in zip(np.ravel(axes), ascii_lowercase):
        ax.text(-0.15, 1.05, lab, transform=ax.transAxes, fontweight = 'bold')
    plt.tight_layout()
    plt.savefig(f'sens_main_{caseid}.png')

    #Total sensitivity indices
    fig, axes = plt.subplots(2, 2, figsize = (14, 11), sharex = True, sharey = True)
    for i, (pft,pftname) in enumerate(zip([2, 3, 11, 0], pft_names[:-1] + ['Column'])):
      subset = np.where(np.array(mycase.ensemble_pfts) == pft)[0]

      ax = axes.flat[i]

      bottom = np.zeros(len(x_pos))
      for s in subset:
        temp = np.array([mycase.sens_tot[v][s,0] for v in variable_list])
        ax.bar(x_pos, temp, align='center', # alpha=0.5,
               bottom = bottom, label = mycase.ensemble_parms[s])
        bottom = bottom + temp

      # add a line for total
      total = np.array([mycase.sens_tot[v][:,0].sum() for v in variable_list])
      ax.plot(x_pos, total, '-k', label = 'Total')

      ax.set_xticks(x_pos)
      ax.set_xticklabels(ticklabels, rotation=90)
      ax.set_title(f'{pftname} parameters')
      ax.legend(loc = [1.05, 0.5])
    for ax, lab in zip(np.ravel(axes), ascii_lowercase):
        ax.text(-0.15, 1.05, lab, transform=ax.transAxes, fontweight = 'bold')
    plt.tight_layout()
    plt.savefig(f'sens_tot_{caseid}.png')


#plot GSA
plot_GSA_treatment()

##Set intial values for parameters
##parms=((np.array(mycase.ensemble_pmax)+np.array(mycase.ensemble_pmin))/2)

##Run MCMC for the 2 varibles of interest
##mycase.MCMC(parms, ['TLAI_ann','QRUNOFF_ann'], 100000)
