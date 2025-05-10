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
UQ_only = True # False - re-generate the "output" dict and overwrite pklfile; 
                # True - directly read from pklfile
if len(suffix) > 0:
  casename=caseid+'_'+site+'_'+compset +'_'+suffix
else:
  casename=caseid+'_'+site+'_'+compset

#Make sure to back up the old pkl files before this step!
#Create case object
if (not UQ_only):
 mycase = model_ELM.ELMcase(caseid=caseid,compset=compset,site=site, \
         sitegroup=sitegroup, machine='cases-baseline', suffix=suffix, \
         runroot=runroot, caseroot=caseroot)

 mycase.casename=casename

 print(mycase.casename)

 mycase.startyear=2015   #Starting year of run; doesn't have to be same as actual run
 mycase.run_n=7         #Number of years for the run (to post process)
 mycase.postproc_pfts=[2,3,11,12]  #PFTs to postprocess
 mycase.postproc_vars=['GPP','GPP_pft','NPP','NPP_pft','QVEGT','QVEGT_pft','NEE','TOTVEGC']
 mycase.postproc_startyear=2015    #Starting year to postprocess/calibrate
 mycase.postproc_endyear= 2021
 mycase.postproc_freq = 'annual'
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
else:
 myfile=open(os.path.join(os.environ['HOME'],'models','OLMT','pklfiles', casename+'.pkl'),'rb')
 mycase=pickle.load(myfile)
 #myfile2=open('./pklfiles/'+casename2+'.pkl','rb')
 #mycase2=pickle.load(myfile2)

# -------------------------------------------------------------------------------------
# Observed values
# -------------------------------------------------------------------------------------
# Put into the dictionary as numpy arrays

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
postproc_only = True
#get the node file and parse
def get_nodelist():
  mynodes=[]
  nodelist=os.environ['SLURM_JOB_NODELIST'].split('xxx')
  print(nodelist)
  for n in nodelist:
    if ('[' in n):
        node_prefix=n.split('[')[0]
        nodelist2=n.split('[')[1].split(',')
        for n2 in nodelist2:
          if ('-' in n2):
            firstnode=n2.split('-')[0]
            lastnode=n2.split('-')[1].strip(']')
            for nn in range(int(firstnode),int(lastnode)+1):
              if ('baseline' in mycase.machine):
                nstr = str(nn)
              else:
                nstr = str(10000+nn)[1:]
              mynodes.append(node_prefix+nstr)
          else:
              if ('baseline' in mycase.machine):
                nstr=str(n2)
              else:
                nstr=str(10000+n2)[1:]
              mynodes.append(node_prefix+nstr)
    else:
        mynodes.append(n)
  return mynodes

def get_node_submit(pactive,process_nodes,mynodes):
    node_submit=0
    for n in range(0,len(mynodes)):
         ctn=0    #Counter for active processes on each node
         for p in range(0,len(processes)):
                if pactive[p] == 1 and process_nodes[p] == n:
                    ctn=ctn+1
         if (ctn < mycase.npernode/mycase.np):
             #If this node is not full, submit
             node_submit=n
    return(node_submit)

def check_run_success(n):
    success=False
    jobst = str(100000+n)
    rundir = mycase.runroot+'/UQ/'+mycase.casename+'/g'+jobst[1:]
    yst = str(10000+mycase.startyear+mycase.run_n)[1:]
    if (os.path.isfile(rundir+'/'+mycase.casename+'.elm.r.'+yst+'-01-01-00000.nc')):
        success=True
        print(success)
    return success

def active_processes(processes,process_jobnum,process_hang):
    """Returns the number of processes that are still running."""
    pactive=[]
    n=0
    for process in processes:
        if process.poll() is None:  # None means the process is still running
            #Check if final restart file created
            pactive.append(1)
            if (check_run_success(process_jobnum[n])):
                process_hang[n] = process_hang[n]+1
            if (process_hang[n] > 30):
                process.kill()  # Force kill the process
        else:
            pactive.append(0)
            #Post-process ensemble member if it hasn't yet been done
            if (mycase.postprocessed[n] == 0):
                print(n, check_run_success(process_jobnum[n]))
                if (check_run_success(process_jobnum[n])):
                    ierr = postprocess_ensemble(process_jobnum[n])
                else:
                    print('Ensemble member '+str(process_jobnum[n])+ \
                            'Failed to complete')
                mycase.postprocessed[n] = 1
        n=n+1
    return pactive


def postprocess_treatment(var, ens_num, startyear, endyear, index, hnum):
  gst = str(100000+ens_num)[1:]
  rundir = mycase.rundir_UQ+'/g'+gst
  os.chdir(rundir)
  lnd_in = open('./lnd_in')
  #Get history file info from lnd_in
  for s in lnd_in:
      if (s.split('=')[0].strip() == 'hist_mfilt'):
          hist_mfilt = int((s.split('=')[1].strip()).split(',')[hnum].strip())
      if (s.split('=')[0].strip() == 'hist_nhtfrq'):
          hist_nhtfrq = int((s.split('=')[1].strip()).split(',')[hnum].strip())
  lnd_in.close()
  if (hist_nhtfrq == 0):
    nperyear=12
  else:
    nperyear = abs(8760/hist_nhtfrq)

  plot_list = ['TAMB','T0.00','T0.00','T0.00CO2','T2.25','T2.25CO2','T4.50','T4.50CO2',
               'T6.75','T6.75CO2','T9.00','T9.00CO2']

  var_out = var
  if ('_pft' in var):
      var_out = var_out+str(index)
  if (not var_out in mycase.output):
      mycase.output[var_out] = np.zeros([len(plot_list), mycase.nsamples], float)

  for pind, plot in enumerate(plot_list):
    file_list = []
    file_list_all = np.sort(glob.glob(plot+'/'+mycase.casename+'.elm.h'+str(hnum)+'.*.nc'))
    #Filter the requested years
    for f in file_list_all:
        if (hist_nhtfrq == 0):
            yr = int(f.split('-')[-2][-4:])
        else:
            yr = int(f.split('-')[-4][-4:])
        if (startyear < 0):
            startyear = yr
        if (endyear >= 9999):
            lastyr = yr
        if (yr >= startyear and yr <= endyear):
            file_list.append(f)
    if (endyear >= 9999):
        endyear=lastyr
        if (nperyear != 12):
          #If not monthly files, ignore the last file (it only represents a single timestep)
          file_list = file_list[:-1]

    os.system('ncrcat -O -v '+var.split('_pft')[0]+' '+' '.join(file_list)+' '+var+'_'+plot+'.nc')
    myoutput = Dataset(var+'_'+plot+'.nc','r')
    if (myoutput[var.split('_pft')[0]][:].ndim == 4):
      #2D output with vertical structure
      values = myoutput[var.split('_pft')[0]][:,index,0,0]
    elif (myoutput[var.split('_pft')[0]][:].ndim == 3):
      #2D output or 1D output with vertical structure (currently assumes 1D)
      values = myoutput[var.split('_pft')[0]][:,index,0]
    else:
      #1D output (unstructured grid)
      if ('_pft' in var):  #PFT-level output
          values = myoutput[var.split('_pft')[0]][:,index]
      else:
          values = myoutput[var.split('_pft')[0]][:,0]

    # all time average
    values_out = np.array([np.mean(values[:])])
    mycase.output[var_out][pind,ens_num-1] = values_out


def postprocess_ensemble(n):
  #Postprocess
  if (mycase.postproc_vars != []):
      for v in mycase.postproc_vars:
        hnum=1 # column leve data ".h1."
        mypfts=[0]
        if ('_pft' in v):
          #PFT level outputs requested
          hnum=2 # PFT level ".h2."
          mypfts=mycase.postproc_pfts
        for p in mypfts:
          postprocess_treatment(v, ens_num=n, startyear=mycase.postproc_startyear, \
                                endyear=mycase.postproc_endyear,index=p,hnum=hnum)
  return 0

if (not UQ_only):
 workdir = os.getcwd()

 processes=[]
 process_jobnum=[]
 process_hang=[]    #Keep track of how long process has been hanging
 mycase.postprocessed=np.zeros([mycase.nsamples],int)
 n_job = 1
 if (mycase.noslurm == False):
    process_nodes = []
    mynodes = get_nodelist()

 #Run the simulations

 while (n_job <= mycase.nsamples):
  pactive = active_processes(processes,process_jobnum,process_hang)
  if (sum(pactive) < int(mycase.np_ensemble)):
    jobst = str(100000+n_job)
    rundir = mycase.runroot+'/UQ/'+mycase.casename+'/g'+jobst[1:]+'/'
    log_file_path = f"{rundir}e3sm_log.txt"
    #log_file_path = '/gpfs/wolf2/cades/cli185/proj-shared/zdr/OLMT/e3sm_log.txt'
    #Copy relevant files
    if not postproc_only:
        mycase.ensemble_copy(n_job)
    with open(log_file_path, "w") as log_file:
       if (mycase.noslurm == False):
         node_submit=get_node_submit(pactive,process_nodes,mynodes)
         command = ['srun -n '+str(mycase.np)+' -c 1 -w '+mynodes[node_submit]+' '+mycase.exeroot+'/e3sm.exe']
         if postproc_only:
             command = 'ls'
         process_nodes.append(node_submit)
       else:
         command = [mycase.exeroot+'/e3sm.exe']
       process = subprocess.Popen(command, shell=True, stderr=subprocess.STDOUT, cwd=rundir, stdout=log_file)
       processes.append(process)
       process_jobnum.append(n_job)
       process_hang.append(0)
    n_job=n_job+1
  else:
    time.sleep(0.1)

 while (sum(pactive) > 0):
    pactive = active_processes(processes,process_jobnum,process_hang)
    time.sleep(0.1)

 mycase.create_pkl(outdir=mycase.OLMTdir+'/pklfiles/')

# -------------------------------------------------------------------------------------
# custom outputs (change units, sum variables, etc)
# -------------------------------------------------------------------------------------
# We can do additional processing of the output time series here. 
## mycase.output['NPP_correct'] = (mycase.output['FATES_NPP']-mycase.output['FATES_EXCESS_RESP'])*24*3600*365*1000
## mycase.output['NUP'] = (mycase.output['FATES_NH4UPTAKE']+mycase.output['FATES_NO3UPTAKE'])*24*3600*365*1000
## mycase.postproc_vars.append('NPP_correct')
## mycase.postproc_vars.append('NUP')


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


def plot_GSA_treatment(myvars):
    plot_list = ['TAMB','T0.00','T0.00','T0.00CO2','T2.25','T2.25CO2','T4.50','T4.50CO2',
                'T6.75','T6.75CO2','T9.00','T9.00CO2']

    for v in myvars:
      #Plot main sensitivity indices
      fig,ax = plt.subplots()
      nvar = mycase.sens_main[v].shape[1]
      x_pos = np.cumsum(np.ones(nvar))
      ax.bar(x_pos, mycase.sens_main[v][0,:], align='center', alpha=0.5)
      ax.set_xticks(x_pos)
      ax.set_xticklabels(plot_list, rotation=45)

      bottom=mycase.sens_main[v][0,:]
      for p in range(1,mycase.nparms_ensemble):
       ax.bar(x_pos, mycase.sens_main[v][p,:], bottom=bottom)
       bottom=bottom+mycase.sens_main[v][p,:]
      plt.legend(mycase.ensemble_parms)
      plt.savefig('sens_main_'+v+'.png')

      #Total sensitivity indices
      fig,ax = plt.subplots()
      ax.bar(x_pos, mycase.sens_tot[v][0,:], align='center', alpha=0.5)
      ax.set_xticks(x_pos)
      ax.set_xticklabels(plot_list, rotation=45)
      bottom=mycase.sens_tot[v][0,:]
      for p in range(1,mycase.nparms_ensemble):
       ax.bar(x_pos, mycase.sens_tot[v][p,:], bottom=bottom)
       bottom=bottom+mycase.sens_tot[v][p,:]
      plt.legend(mycase.ensemble_parms)
      plt.savefig('sens_tot_'+v+'.png')


#plot GSA
mycase.plot_GSA(variable_list)

#Save postprocessed output
mycase.create_pkl(outdir=mycase.OLMTdir+'/pklfiles/')

##Set intial values for parameters
##parms=((np.array(mycase.ensemble_pmax)+np.array(mycase.ensemble_pmin))/2)

##Run MCMC for the 2 varibles of interest
##mycase.MCMC(parms, ['TLAI_ann','QRUNOFF_ann'], 100000)

