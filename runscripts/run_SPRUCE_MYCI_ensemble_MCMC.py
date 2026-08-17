"""Python code used to perform post-processing of SPRUCE treatment output and parameter inversion.
   Assemble all the treatment ensemble runs. Cannot run the ensembles. 
"""
import sys,os, time
sys.path.append('..')
import numpy as np
import subprocess
import pickle
import model_ELM
from optparse import OptionParser
import pandas as pd

parser = OptionParser()

case_prefix = '20260809'
case_suffix = 'MYCI_full'
mode = 'MYCI' # 'RD' or 'MYCI'
case_base = f'{case_prefix}_US-SPR_ICB20TRCNP{mode}CTCBC_{case_suffix}'
case_treatments = [f'{case_prefix}_US-SPR_ICB20TRCNP{mode}CTCBC_{trt}_{case_suffix}' for trt in \
    ['TAMB','T0.00','T2.25','T4.50','T6.75','T9.00','T0.00eCO2','T2.25eCO2', \
     'T4.50eCO2','T6.75eCO2','T9.00eCO2']]
plot_treatments = ['plot07', 'plot06', 'plot20', 'plot13', 'plot08', 'plot17', 'plot19', 
                   'plot11', 'plot04', 'plot16', 'plot10']
case_plots = []

postproc_startyear = 2015
postproc_endyear   = 2023
postproc_freq      = 'annual'   #Can be daily, monthly, annual
postproc_sim       = False # True - postprocess each treatment run to obtain the necessary outputs
                           # False - directly read the postprocessed model output from each treatment run

# Load the 20th century case's config. The treatment cases will be used to
# re-populate the output
myfile=open(os.path.join(os.environ['HOME'],'models','OLMT','pklfiles', case_base+'.pkl'),'rb')
mycase=pickle.load(myfile)

# -------------------------------------------------------------------------------------
# Observed values
# -------------------------------------------------------------------------------------
chamber_list_complete = [7, 6, 20, 13, 8, 17, 19, 11, 4, 16, 10] # 21,
chamber_list_names_complete = ["TAMB", "T0.00", "T2.25", "T4.50", "T6.75", "T9.00",
                               "T0.00eCO2", "T2.25eCO2", "T4.50eCO2", "T6.75eCO2", "T9.00eCO2"]
chamber_list_complete_dict = dict({f"P{a:02d}": b for  a,b in zip(chamber_list_complete, chamber_list_names_complete)})
chamber_list_complete_dict2 = dict({a: b for  a,b in zip(chamber_list_complete, chamber_list_names_complete)})

myobsdir='/projects/hpcl-cli185/proj-shared/ywo/ELM_Allocation/output'

cobs = pd.read_csv(myobsdir+'/extract_obs_productivity.csv', index_col=[0,1]).drop(['eCO2','NPP','NEE'],axis=1) # skip parameters not useful in calibration

# drop BGNPP except 2016-2017, because other years were not real measurements
cobs.loc[(cobs.index.get_level_values('Year') < 2016) | (cobs.index.get_level_values('Year') > 2017), 'BGNPP_TreeShrub'] = np.nan

tlai = pd.read_csv(myobsdir+'/extract_tlai_tls.csv', index_col=[0,1])
myobs = pd.concat([cobs, tlai], axis=1).rename(lambda x: chamber_list_complete_dict[x], level='Plot').sort_index()


## Based on discussion with Claude, because I already removed pre-treatment
## unevenness from the biomass observations, while pre-treatment unevenness
## in the flux observations are already aliased in the inter-plot availability,
## no extra pre-treatment uncertainty propagation is needed. 
## 
## inherent observation uncertainty based on pre-treatment levels (sd / mean)
# calibrate on mean & slope ##, with uncertainty propagated
rel_sd = {'Tair': 0.5/6.25, 'AGNPP_Spruce': 32/73, 'AGNPP_Tamarack': 32/73,
            'AGNPP_Shrub': 34/104, 'NPP_moss': 67/208, 'BGNPP_TreeShrub': 4.8/3.4,
            'HR': 53/283, 'AGBiomass_Spruce': 439.694598/748.5575,
            'AGBiomass_Tamarack': 143.243042/215.0901, 'AGBiomass_Shrub': 99.036822/243.5690,
            'LAImax_Spruce': 0.671678/1.472463, 'LAImax_Tamarack': 0.268715/0.469105,
            'LAImax_Shrub': 0.408444/1.825093}

def regress_uncert(df, xvar='Tair', rel_sd=rel_sd, n=5000, seed=0):
    rng = np.random.default_rng(seed)

    rows = {}

    for v in df.columns.drop(xvar):
        d = df[[xvar, v]].dropna()
        if len(d) < 3:
            continue

        x = d[xvar].to_numpy()
        y = d[v].to_numpy()

        sx = rel_sd[xvar] * abs(x.mean())
        sy = rel_sd[v] * abs(y.mean())
        xs = x + rng.normal(0, sx, (n, x.size))
        ys = y + rng.normal(0, sy, (n, y.size))
        
        xc = xs - xs.mean(axis=1, keepdims=True)
        yc = ys - ys.mean(axis=1, keepdims=True)
        b_mc = (xc * yc).sum(axis=1) / (xc**2).sum(axis=1)
        a_mc = ys.mean(axis=1) - b_mc * xs.mean(axis=1)

        x_mean = x.mean(); y_mean = y.mean()
        xc0 = x - x_mean; yc0 = y - y_mean
        sxx = (xc0**2).sum()

        # Skip regressions for which the predictor is constant.
        if sxx == 0:
            continue

        b = (xc0 * yc0).sum() / sxx
        a = y_mean - b * x_mean
        resid = y - (a + b * x)

        residual_variance = (resid**2).sum() / (len(d) - 2)

        rows[v] = dict(
            n=len(d),
            slope=b, slope_sd_meas=b_mc.std(ddof=1),
            slope_se_fit=np.sqrt(residual_variance / sxx),
            intercept=a, intercept_sd_meas=a_mc.std(ddof=1),
            intercept_se_fit=np.sqrt(
                residual_variance * (1 / len(d) + x_mean**2 / sxx)
            ),
        )

    return pd.DataFrame(rows).T


res = regress_uncert(myobs.loc[['TAMB','T0.00','T2.25','T4.50','T6.75','T9.00'],:])
res_co2 = regress_uncert(myobs.loc[['T0.00eCO2','T2.25eCO2','T4.50eCO2','T6.75eCO2','T9.00eCO2'],:])

# populate the obs object
# this round of optimization does not focus on biomass
for vv in ['AGNPP_Spruce','AGNPP_Tamarack','AGNPP_Shrub','NPP_moss','BGNPP_TreeShrub','HR']:
   mycase.obs[f'{vv}_slope'] = np.array([float(res.loc[vv, 'slope']), float(res_co2.loc[vv, 'slope'])])
   mycase.obs_err[f'{vv}_slope'] = np.array([float(res.loc[vv, 'slope_se_fit']), 
                                             float(res_co2.loc[vv, 'slope_se_fit'])])

   mycase.obs[f'{vv}_intercept'] = np.array([float(res.loc[vv, 'intercept']), float(res_co2.loc[vv, 'intercept'])])
   mycase.obs_err[f'{vv}_intercept'] = np.array([float(res.loc[vv, 'intercept_se_fit']), 
                                                 float(res_co2.loc[vv, 'intercept_se_fit'])])
   ##mycase.obs[f'{vv}_T0.00'] = myobs.loc['T0.00', vv].mean() # not enough measurement in TAMB
   ##mycase.obs_err[f'{vv}_T0.00'] = rel_sd[vv] * myobs.loc['T0.00', vv].mean() # not enough measurement in TAMB

# -------------------------------------------------------------------------------------
# Loop through the treatment pklfiles to obtain the custom output
# -------------------------------------------------------------------------------------
def regress_vectorized(x, y):
    """
    Apply ordinary least-squares linear regression to every column of y.

    Model:
        y[:, j] = intercept[j] + slope[j] * x

    Parameters
    ----------
    x : array_like, shape (n_samples,)
    y : array_like, shape (n_samples, n_regressions)

    Returns
    -------
    slope : ndarray, shape (n_regressions,)
    intercept : ndarray, shape (n_regressions,)
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    if x.ndim != 1:
        raise ValueError("x must be a 1D array")
    if y.ndim != 2:
        raise ValueError("y must be a 2D array")
    if y.shape[0] != x.size:
        raise ValueError("x and y must have the same first dimension")
    if x.size < 2:
        raise ValueError("At least two observations are required")

    x_centered = x - x.mean()
    denominator = x_centered @ x_centered

    if denominator == 0:
        raise ValueError("Regression is undefined because x has zero variance")

    # OLS equation: slope = Σ((x - x̄)(y - ȳ)) / Σ((x - x̄)²)
    slope = (x_centered @ y) / denominator # centering y is unecessary because (x - x̄) sums to zero
    intercept = y.mean(axis=0) - slope * x.mean()

    return slope.reshape(1,-1), intercept.reshape(1,-1)


sphagnum_fraction = pd.read_excel(os.path.join(os.environ['SHARDIR'], 'ELM_Allocation', 
                                               'data', 'Sphagnum_fraction.xlsx'), skiprows=1, index_col=[1])
sphagnum_fraction.index = [f'plot{p:02d}' for p in sphagnum_fraction.index]
sphagnum_fraction[2015] = sphagnum_fraction[2016].copy()
sphagnum_fraction = sphagnum_fraction.drop(['filename','Temp','CO2'],axis=1)
sphagnum_fraction = sphagnum_fraction.T.sort_index()

tair = []
tair_co2 = []
output = {}
for vv in ['AGNPP_Spruce','AGNPP_Tamarack','AGNPP_Shrub','NPP_moss','BGNPP_TreeShrub','HR']:
   output[vv] = []
   output[f'{vv}_CO2'] = []
for case, plot in zip(case_treatments, plot_treatments):

  with open(os.path.join(os.environ['HOME'],'models','OLMT','pklfiles', case+'.pkl'),'rb') as f:
    temp=pickle.load(f)

  metfile = temp.inputdata_path + '/atm/datm7/CLM1PT_data/SPRUCE_data/' + plot + '/all_hourly.nc'
  temp_tair = temp.getncvar(metfile, 'TBOT').reshape(-1, 365*48).mean(axis=1)[:7].data - 273.15 # 2015-2021
  if 'CO2' in case:
    tair_co2.append(temp_tair)
  else:
    tair.append(temp_tair)

  if 'CO2' in case:
    suffix='_CO2'
  else:
    suffix=''

  output['AGNPP_Spruce'+suffix].append(temp.output['AGNPP_pft2'] * 0.36 * 86400 * 365)
  output['AGNPP_Tamarack'+suffix].append(temp.output['AGNPP_pft3'] * 0.14 * 86400 * 365)
  output['AGNPP_Shrub'+suffix].append(temp.output['AGNPP_pft11'] * 0.25 * 86400 * 365)

  output['NPP_moss'+suffix].append(temp.output['NPP_pft12'] * sphagnum_fraction.loc[:, [plot]].values / 100 * 86400 * 365)

  output['BGNPP_TreeShrub'+suffix].append(temp.output['FROOTC_ALLOC_pft2'] * 0.36 * 86400 * 365 + \
                                   temp.output['FROOTC_ALLOC_pft3'] * 0.14 * 86400 * 365 + \
                                   temp.output['FROOTC_ALLOC_pft11'] * 0.25 * 86400 * 365)
  output['HR'+suffix].append(temp.output['HR'] * 86400 * 365)
tair = np.concatenate(tair)
tair_co2 = np.concatenate(tair_co2)
for vv in output.keys():
   output[vv] = np.concatenate(output[vv], axis=0)


# populate the output object
mycase.output = {}
for vv in ['AGNPP_Spruce','AGNPP_Tamarack','AGNPP_Shrub','NPP_moss','BGNPP_TreeShrub','HR']:
   slope, intercept = regress_vectorized(tair, output[vv])
   slope_co2, intercept_co2 = regress_vectorized(tair_co2, output[vv+'_CO2'])

   mycase.output[f'{vv}_slope'] = np.vstack([slope,slope_co2])
   mycase.output[f'{vv}_intercept'] = np.vstack([intercept,intercept_co2])

#------UQ -----------------------------

#Train surrogate models
myvars = list(mycase.output.keys())
mycase.train_surrogate(myvars)

#Save surrogate models because this takes a long time
mycase.create_pkl(outdir=mycase.OLMTdir+'/pklfiles/')

#run GSA
mycase.GSA(myvars)
#plot GSA
mycase.plot_GSA(myvars)

#Set intial values for parameters
parms=((np.array(mycase.ensemble_pmax)+np.array(mycase.ensemble_pmin))/2)

#Run MCMC for the varibles of interest
mycase.MCMC(parms, myvars, 100000)

#Save postprocessed output
mycase.create_pkl(outdir=mycase.OLMTdir+'/pklfiles/')
