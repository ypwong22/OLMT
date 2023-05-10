from netCDF4 import Dataset
from sklearn.metrics import f1_score
from sklearn.utils import shuffle
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier # RandomForestRegressor doesn't work well compared to MLPRegressor
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.neural_network import MLPRegressor, MLPClassifier # doesn't work well compared to RandomForestClassifier
from scipy.stats.qmc import LatinHypercube
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os, math, sys
import numpy as np
#from mpi4py import MPI
import pickle
from optparse import OptionParser
from xgboost import XGBClassifier, XGBRegressor
from sklearn.multioutput import MultiOutputRegressor


parser = OptionParser()
parser.add_option("--case", dest="casename", default="", \
                  help="Name of case")
(options, args) = parser.parse_args()

UQ_output = 'UQ_output/'+options.casename
datapath = UQ_output+'/data/'
os.system('mkdir -p '+UQ_output+'/NN_surrogate')

print(datapath+'/ptrain.dat')
ptrain = np.loadtxt(datapath+'/ptrain.dat')
ytrain = np.loadtxt(datapath+'/ytrain.dat')
pval   = np.loadtxt(datapath+'/pval.dat')
yval   = np.loadtxt(datapath+'/yval.dat')
varnames_file = open(datapath+'/outnames.txt')
outnames=[]
for s in varnames_file:
  outnames.append(s)
varnames_file.close()

nparms = ptrain.shape[1]
nqoi   = ytrain.shape[1]

good = np.where(ytrain[:,1].squeeze() > -9999)[0]
ytrain = ytrain[good,:].copy()
ptrain = ptrain[good,:].copy()
ntrain = ptrain.shape[0]
print(ntrain, 'Training points')
good = np.where(yval[:,1].squeeze() > -9999)[0]
yval = yval[good,:].copy()
pval = pval[good,:].copy()
nval = pval.shape[0]
print(nval, 'Validation points')


#Normalize parameters
prange = np.zeros([2,nparms],float)
for i in range(0,nparms):
  prange[0,i] = min(ptrain[:,i])
  prange[1,i] = max(ptrain[:,i])

ptrain_norm = ptrain.copy()
pval_norm   = pval.copy()
for i in range(0,nparms):
  ptrain_norm[:,i] = (ptrain[:,i] - prange[0,i])/(prange[1,i] - prange[0,i])
  pval_norm[:,i]   = (pval[:,i] - prange[0,i])/(prange[1,i] - prange[0,i])
  for j in range(0,nval):
    pval_norm[j,i] = min(max(pval_norm[j,i], 0.0), 1.0)


#Label: true = all PFTs grow
ytrain_norm_lab = ytrain[:,1].squeeze() < 9998
ptrain_norm_lab = ptrain_norm
yval_norm_lab = yval[:,1].squeeze() < 9998
pval_norm_lab = pval_norm

#Values: only data points where all PFTs grow
ytrain_val = ytrain[ytrain_norm_lab, :].copy()
ptrain_norm_val = ptrain_norm[ytrain_norm_lab, :].copy()

yval_val = yval[yval_norm_lab, :].copy()
pval_norm_val = pval_norm[yval_norm_lab, :].copy()


#Normalize variables
yrange_val = np.zeros([2,nqoi],float)
qoi_good = []
for i in range(0,nqoi):
  yrange_val[0,i] = min(ytrain_val[:,i])
  yrange_val[1,i] = max(ytrain_val[:,i])
  if (yrange_val[0,i] != yrange_val[1,i]):
    qoi_good.append(i)
np.savetxt(UQ_output+'/NN_surrogate/qoi_good.txt',np.array(qoi_good))
np.savetxt(UQ_output+'/NN_surrogate/yrange_val.txt',np.array(yrange_val))


ytrain_norm_val = ytrain_val.copy()
yval_norm_val   = yval_val.copy()
for i in range(0,nqoi):
  if (yrange_val[0,i] != yrange_val[1,i]):
    ytrain_norm_val[:,i] = (ytrain_val[:,i] - yrange_val[0,i])/(yrange_val[1,i]-yrange_val[0,i])  
    yval_norm_val[:,i]   = (yval_val[:,i]  -  yrange_val[0,i])/(yrange_val[1,i]-yrange_val[0,i])
    for j in range(0,sum(yval_norm_lab)):
      yval_norm_val[j,i] = min(max(yval_norm_val[j,i], 0.0), 1.0)

print(f'All PFTs grow points, training: {sum(ytrain_norm_lab)}/{len(ytrain_norm_lab)}, validation {sum(yval_norm_lab)}/{len(yval_norm_lab)}')


"""
if (sum(ytrain_norm_lab) > 0) & (sum(yval_norm_lab) > 0) & (sum(~ytrain_norm_lab) > 0) & (sum(~yval_norm_lab) > 0):
  #Latin hypercube sampling for parameters
  n_samples = 100
  samples = LatinHypercube(5, seed = 200).random(n_samples)

  params_names = ['n_estimators', 'max_depth', 'min_samples_split', 'min_samples_leaf', 'max_features']
  params_samples = np.full((n_samples, 5), 0) # int
  for i in range(n_samples):
      params_samples[i, 0] = 100 + int(samples[i, 0] * 1500)
      params_samples[i, 1] = 3 + int(samples[i, 1] * (nparms - 3))
      params_samples[i, 2] = 2 + int(samples[i, 2] * 10)
      params_samples[i, 3] = 2 + int(samples[i, 3] * 5)
      params_samples[i, 4] = 3 + int(samples[i, 4] * (nparms - 3))

  brier_best = 99999

  for i in range(n_samples):
    params = dict(zip(params_names, params_samples[i, :]))
    clc = RandomForestClassifier(**params, n_jobs = -1)
    clc.fit(ptrain_norm_lab, ytrain_norm_lab)

    ypredict_train_lab = clc.predict_proba(ptrain_norm_lab)[:, 1]
    ypredict_val_lab = clc.predict_proba(pval_norm_lab)[:, 1]

    brier_train = np.sqrt(np.mean(np.power(ypredict_train_lab - ytrain_norm_lab, 2)))
    brier_val = np.sqrt(np.mean(np.power(ypredict_val_lab - yval_norm_lab, 2)))

    print(f'validation {i}, brier score: {brier_val}')

    if (brier_val < brier_best):
      myfile = open(UQ_output+'/NN_surrogate/fitstats_lab.txt','w')
      myfile.write('Number of parameters:                        '+str(nparms)+'\n')
      myfile.write('Number of outputs:                           1\n')
      myfile.write('Number of training samples:                  '+str(ntrain)+'\n')
      myfile.write('Number of training samples, all PFTs grow:   '+str(sum(ytrain_norm_lab))+'\n')
      myfile.write('Number of validation samples:                '+str(nval)+'\n\n')
      myfile.write('Number of validation samples, all PFTs grow: '+str(sum(yval_norm_lab))+'\n')
      myfile.write('Best random forest classifier:\n')
      for j, par in enumerate(params_names):
        myfile.write(par+': '+str(params_samples[i,j])+'\n')
      myfile.write('Brier score: '+str(brier_val)+'\n')
      brier_best = brier_val
      pkl_filename = UQ_output+'/NN_surrogate/RFmodel.pkl'
      with open(pkl_filename,'wb') as file:
        pickle.dump(clc, file)
      myfile.close()

      ypredict_val_best = ypredict_val_lab
      plt.clf()
      plt.boxplot([ypredict_val_best[~yval_norm_lab], ypredict_val_best[yval_norm_lab]], positions = [0, 1])
      plt.gca().set_xticks([0, 1])
      plt.xlabel('Model all PFTs grow (N/Y)')
      plt.ylabel('Surrogate probability')
      plt.savefig(UQ_output+'/NN_surrogate/rffit.pdf')

    if (brier_best < 0.01):
      print('Brier score < 0.01')
      break
"""

#Train the classification for non-zero PFTs
if (sum(ytrain_norm_lab) > 0) & (sum(yval_norm_lab) > 0) & (sum(~ytrain_norm_lab) > 0) & (sum(~yval_norm_lab) > 0):
  brier_best = 99999

  for n in range(1):
    nmin = 10*np.sqrt(n+1) #max(10, ntrain/20)
    nmax = 20*np.sqrt(n+1) #min(ntrain/4, 100)
    nl = int(np.random.uniform(nmin,nmax))
    nl2 = int(np.random.uniform(nmin,nmax))*2
    do3 = 0 #np.random.uniform(0,1)
    nl3 = int(np.random.uniform(nmin,nmax))
    if (do3 > 0.5):
      clc = MLPClassifier(solver='adam', activation = 'relu', early_stopping=True, tol=1e-7, hidden_layer_sizes=(nl,nl2,nl3,), max_iter=500, validation_fraction=0.2)
    else:
      clc = MLPClassifier(solver='adam', activation = 'relu', early_stopping=True, tol=1e-7, hidden_layer_sizes=(nl,nl2,), max_iter=500, validation_fraction=0.2)

    #clc = GaussianProcessClassifier(kernel = 1.0 * RBF(0.5), random_state = 0)
    #clc = RandomForestClassifier()
    #clc.fit(ptrain_norm_lab, ytrain_norm_lab)
    clc = XGBClassifier(n_estimators=200, max_depth=int(nparms*1/4), learning_rate=0.0001, objective='binary:logitraw')

    temp = np.concatenate([ptrain_norm_lab, np.tile(ptrain_norm_lab[ytrain_norm_lab, :], (6,1))], axis = 0)
    temp2 = np.concatenate([ytrain_norm_lab, np.ones(sum(ytrain_norm_lab) * 6)], axis = 0)
    temp, temp2 = shuffle(temp, temp2, random_state = 0)
    clc.fit(temp, temp2)

    ypredict_train_lab = clc.predict(ptrain_norm_lab)
    ypredict_val_lab   = clc.predict(pval_norm_lab)
    print(f1_score(ypredict_train_lab, ytrain_norm_lab))
    print(f1_score(ypredict_val_lab, yval_norm_lab))

    #ypredict_train_lab = clc.predict_proba(ptrain_norm_lab)[:, 1]
    #ypredict_val_lab   = clc.predict_proba(pval_norm_lab)[:, 1]
    brier_train = np.sqrt(np.mean(np.power(ypredict_train_lab - ytrain_norm_lab, 2)))
    brier_val = np.sqrt(np.mean(np.power(ypredict_val_lab - yval_norm_lab, 2)))
    print(f'validation {n}, brier score: {brier_train} {brier_val}')

    if (brier_val < brier_best):
      myfile = open(UQ_output+'/NN_surrogate/fitstats_lab.txt','w')
      myfile.write('Number of parameters:                        '+str(nparms)+'\n')
      myfile.write('Number of outputs:                           1\n')
      myfile.write('Number of training samples:                  '+str(ntrain)+'\n')
      myfile.write('Number of training samples, all PFTs grow:   '+str(sum(ytrain_norm_lab))+'\n')
      myfile.write('Number of validation samples:                '+str(nval)+'\n\n')
      myfile.write('Number of validation samples, all PFTs grow: '+str(sum(yval_norm_lab))+'\n')
      myfile.write('Best neural network:\n')
      myfile.write('Size of NN layer 1: '+str(nl)+'\n')
      myfile.write('Size of NN layer 2: '+str(nl2)+'\n')
      if (do3 == 1):
        myfile.write('Size of NN layer 3: '+str(nl3)+'\n')
      myfile.write('Brier score: '+str(brier_val)+'\n')
      brier_best = brier_val
      pkl_filename = UQ_output+'/NN_surrogate/RFmodel.pkl'
      with open(pkl_filename,'wb') as file:
        pickle.dump(clc, file)
      myfile.close()

      ypredict_val_best = ypredict_val_lab
      plt.clf()
      plt.boxplot([ypredict_val_best[~yval_norm_lab], ypredict_val_best[yval_norm_lab]], positions = [0, 1])
      plt.gca().set_xticks([0, 1])
      plt.xlabel('Model all PFTs grow (N/Y)')
      plt.ylabel('Surrogate probability')
      plt.savefig(UQ_output+'/NN_surrogate/rffit.pdf')

    if (brier_best < 0.01):
      print('Brier score < 0.01')
      break


def gen_circle():
    "Generate a sample dataset that y is a 2 dim circle."
    rng = np.random.RandomState(1994)
    X = np.sort(200 * rng.rand(100, 1) - 100, axis=0)
    y = np.array([np.pi * np.sin(X).ravel(), np.pi * np.cos(X).ravel()]).T
    y[::5, :] += 0.5 - rng.rand(20, 2)
    y = y - y.min()
    y = y / y.max()
    return X, y


#Train the non-zero values
rmse_best = 9999
corr_best = 0
for n in range(100):
  """
  nmin = 10*np.sqrt(n+1) #max(10, ntrain/20)
  nmax = 20*np.sqrt(n+1) #min(ntrain/4, 100)
  nl = int(np.random.uniform(nmin,nmax))
  nl2 = int(np.random.uniform(nmin,nmax))*2
  do3 = 0 #np.random.uniform(0,1)
  nl3 = int(np.random.uniform(nmin,nmax))
  if (do3 > 0.5):
    clf = MLPRegressor(solver='adam', early_stopping=True, tol=1e-7, hidden_layer_sizes=(nl,nl2,nl3,), max_iter=500, validation_fraction=0.2)
  else:
    clf = MLPRegressor(solver='adam', early_stopping=True, tol=1e-7, hidden_layer_sizes=(nl,nl2,), max_iter=500, validation_fraction=0.2)
  """
  clf = MultiOutputRegressor(XGBRegressor(n_estimators=int(30*np.sqrt(n+1)), max_depth=int(nparms*1/4), learning_rate=1e-5, objective='reg:squaredlogerror'))
  clf.fit(ptrain_norm_val, ytrain_norm_val[:,qoi_good])

  ypredict_train_temp = clf.predict(ptrain_norm_val)
  ypredict_val_temp   = clf.predict(pval_norm_val)
  ypredict_train_val = ytrain_norm_val.copy()
  ypredict_val_val = yval_norm_val.copy()
  ypredict_train_val[:,qoi_good] = ypredict_train_temp
  ypredict_val_val[:,qoi_good] = ypredict_val_temp

  corr_train=[]
  rmse_train = []
  corr_val=[]
  rmse_val=[]
  for qoi in qoi_good:
    corr_train.append((np.corrcoef(ytrain_norm_val.astype(float)[:,qoi], ypredict_train_val.astype(float)[:,qoi])[0,1])**2)
    rmse_train.append((sum((ypredict_train_val[:,qoi]-ytrain_norm_val[:,qoi])**2)/ntrain)**0.5)
    corr_val.append((np.corrcoef(yval_norm_val.astype(float)[:,qoi], ypredict_val_val.astype(float)[:,qoi])[0,1])**2)
    rmse_val.append((sum((ypredict_val_val[:,qoi]-yval_norm_val[:,qoi])**2)/nval)**0.5)
  print(f'validation {n}, mean corr: {np.mean(corr_val)}, mean rmse: {np.mean(rmse_val)}')

  if (sum(corr_val) > corr_best):
    myfile = open(UQ_output+'/NN_surrogate/fitstats_val.txt','w')
    myfile.write('Number of parameters:         '+str(nparms)+'\n')
    myfile.write('Number of outputs:            '+str(nqoi)+'\n')
    myfile.write('Number of good outputs:       '+str(len(qoi_good))+'\n')
    myfile.write('Number of training samples:   '+str(ntrain)+'\n')
    myfile.write('Number of validation samples: '+str(nval)+'\n\n')
    #myfile.write('Best neural network:\n')
    #myfile.write('Size of NN layer 1: '+str(nl)+'\n')
    #myfile.write('Size of NN layer 2: '+str(nl2)+'\n')
    #if (do3 == 1):
    #  myfile.write('Size of NN layer 3: '+str(nl3)+'\n')
    for q in range(0, len(qoi_good)):
      myfile.write('QOI validation '+str(qoi_good[q])+' (R2,rmse): '+str(corr_val[q])+'  '+ \
             str(rmse_val[q]**2)+'\n')
    corr_best = sum(corr_val)
    pkl_filename = UQ_output+'/NN_surrogate/NNmodel.pkl'
    ypredict_val_best = ypredict_val_val
    with open(pkl_filename,'wb') as file:
      pickle.dump(clf, file)
    for q in qoi_good:
      plt.clf()
      plt.scatter(yval_norm_val[:,q]*(yrange_val[1,q]-yrange_val[0,q])+yrange_val[0,q], ypredict_val_best[:,q]*(yrange_val[1,q]-yrange_val[0,q])+yrange_val[0,q])
      plt.xlabel('Model '+outnames[q])
      plt.ylabel('Surrogate '+outnames[q])
      plt.savefig(UQ_output+'/NN_surrogate/nnfit_qoi'+str(q)+'.pdf')
    myfile.close()
  if (min(corr_val) > 0.99):
    print('All QOIs have R2 > 0.99')
    break
