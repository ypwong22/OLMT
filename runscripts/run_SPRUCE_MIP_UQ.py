import sys,os
sys.path.append('..')
import numpy as np
import pickle
import model_ELM
from optparse import OptionParser
import matplotlib.pyplot as plt
from string import ascii_lowercase


def plot_GSA_treatment(plot):
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
    fig, axes = plt.subplots(2, 3, figsize = (20, 11), sharex = True, sharey = True)
    fig.subplots_adjust(wspace = 0.5)
    for i, (pft,pftname) in enumerate(zip([2, 3, 11, 12, 0], pft_names + ['Column'])):
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
    for ax, lab in zip(np.ravel(axes)[:-1], ascii_lowercase):
        ax.text(-0.15, 1.05, f'{lab})', transform=ax.transAxes, fontsize = 16)
    axes.flat[-1].axis('off')
    plt.savefig(f'sens_main_20240723_SPRUCE_MIP_{plot}.png')

    #Total sensitivity indices
    fig, axes = plt.subplots(2, 3, figsize = (20, 11), sharex = True, sharey = True)
    fig.subplots_adjust(wspace = 0.5)
    for i, (pft,pftname) in enumerate(zip([2, 3, 11, 12, 0], pft_names + ['Column'])):
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
    for ax, lab in zip(np.ravel(axes)[:-1], ascii_lowercase):
        ax.text(-0.15, 1.05, f'{lab})', transform=ax.transAxes, fontsize = 16)
    axes.flat[-1].axis('off')
    plt.savefig(f'sens_tot_20240723_SPRUCE_MIP_{plot}.png')


#Python code used to manage the ensemble simulations 
#  and perform post-processing of model output.

parser = OptionParser()

runroot='/gpfs/wolf2/cades/cli185/proj-shared/zdr/SPRUCE/e3sm_run/UQ'

VAR_COL = ['NEE','NPP','CH4PROD','FCH4','HR']
VAR_PFT = ['NPP', 'TLAI', 'TOTVEGC']
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
regen_pkl = False
#PLOT_LIST = ['TAMB','T0.00','T0.00eCO2','T2.25','T2.25eCO2','T4.50','T4.50eCO2',
#             'T6.75','T6.75eCO2','T9.00','T9.00eCO2']
PLOT_LIST = ['T0.00','T0.00eCO2']
#PLOT_LIST = ['T2.25','T2.25eCO2','T4.50','T4.50eCO2']
#PLOT_LIST = ['T6.75','T6.75eCO2','T9.00','T9.00eCO2'] 
for i, plot in enumerate(PLOT_LIST):
  
    if regen_pkl:

        myfile=open('/ccsopen/home/zdr/models/OLMT/pklfiles/20240723_US-SPR_ICB20TRCNPRDCTCBC' + \
                    f'_{plot}.pkl','rb')
        mycase=pickle.load(myfile)
        myfile.close()

        print(mycase.casename)

        # mycase.samples & mycase.ensemble_parms & mycase.nsamples are already in 

        mycase.postproc_pfts=[]
        mycase.postproc_vars=[]
        mycase.postproc_startyear=2015    #Starting year to postprocess/calibrate
        mycase.postproc_endyear= 2021
        mycase.postproc_freq = 'annual'

        mycase.obs={}
        mycase.obs_err={}
        mycase.pscaler={}
        mycase.yscaler={}

        mycase.OLMTdir=os.path.join(os.environ['HOME'],'models','OLMT')
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

        BLOCK = 160
        collection_all = np.full([mycase.nsamples, nvars], np.nan)
        for part in range(mycase.nsamples // BLOCK):
            temp = np.load(os.path.join(os.environ['PROJDIR'], 'ELM_Phenology', 'output', 'extract',
                                        '20240723', f'uncertainty_extraction_part{part:03g}.bin'),
                            allow_pickle=True)
            # Only take the values at this plot
            collection_all[(part*BLOCK):(part*BLOCK+BLOCK), :] = temp[:, :, i]

        mycase.output = {}
        for v, var in enumerate(VAR_COL):
            # reshape to make sure the second dimension is ensemble member
            mycase.output[var] = collection_all[:, v].reshape(1,-1)
        for v, var in enumerate(VAR_PFT):
            for t, pft in enumerate(pft_list):
                # reshape to make sure the second dimension is ensemble member
                mycase.output[var+'_pft'+str(pft)] = collection_all[:,
                                len(VAR_COL)+v*len(pft_list)+t].reshape(1,-1)

        #------UQ -----------------------------

        #Train surrogate models

        mycase.train_surrogate(variable_list)

        #run GSA
        mycase.GSA(variable_list)

        #Save postprocessed output
        mycase.create_pkl(outdir=mycase.OLMTdir+'/pklfiles/')

    else:
        myfile=open(os.path.join(os.environ['HOME'],'models','OLMT', 'pklfiles',
                                 '20240723_US-SPR_ICB20TRCNPRDCTCBC' + f'_{plot}.pkl'),'rb')
        mycase=pickle.load(myfile)
        myfile.close()

        print(mycase.casename)

    #plot GSA
    plot_GSA_treatment(plot)
