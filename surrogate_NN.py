from netCDF4 import Dataset
from sklearn.metrics import recall_score, precision_score
from sklearn.neural_network import MLPRegressor
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os, math, sys
import numpy as np
#from mpi4py import MPI
import pickle
from optparse import OptionParser
from xgboost import XGBClassifier # MLPClassifier, RandomForestClassifier do not work well on shrubs


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


#Normalize variables
yrange = np.zeros([2,nqoi],float)
qoi_good = [] # need to label if all PFTs grew or not
for i in range(0,nqoi):
  yrange[0,i] = min(ytrain[:,i])
  yrange[1,i] = max(ytrain[:,i])
  if (yrange[0,i] != yrange[1,i]):
    qoi_good.append(i)
np.savetxt(UQ_output+'/NN_surrogate/qoi_good.txt',np.array(qoi_good))
np.savetxt(UQ_output+'/NN_surrogate/yrange.txt',np.array(yrange))

ytrain_norm = ytrain.copy()
yval_norm   = yval.copy()
for i in range(0,nqoi):
  if (yrange[0,i] != yrange[1,i]):
    ytrain_norm[:,i] = (ytrain[:,i] - yrange[0,i])/(yrange[1,i]-yrange[0,i])
    yval_norm[:,i]   = (yval[:,i]  -  yrange[0,i])/(yrange[1,i]-yrange[0,i])
    for j in range(0,yval_norm.shape[0]):
      yval_norm[j,i] = min(max(yval_norm[j,i], 0.0), 1.0)

#For the variables that are zero if the PFT did not grow:
#label: 0 = value is non-zero (PFT grew), 1 = value is zero (PFT did not grow)
ytrain_norm_label = np.zeros(ytrain.shape)
yval_norm_label = np.zeros(yval.shape)
qoi_need_lab = []
for i in range(0,nqoi):
  ymin = ytrain[:,i].min()
  ymax = ytrain[:,i].max()
  if ~((ymin < 0) & (ymax > 0)):
    near_zero = (np.abs(ytrain[:,i]) / max(abs(ymin), abs(ymax))) < 1e-8
    # more than 10% is zero
    if sum(near_zero) > (0.1 * nval):
      qoi_need_lab.append(i) # values do not naturally span zero
      ytrain_norm_label[near_zero, i] = 1
      near_zero2 = (np.abs(yval[:,i]) / max(abs(ymin), abs(ymax))) < 1e-8
      yval_norm_label[near_zero2, i] = 1
np.savetxt(UQ_output+'/NN_surrogate/qoi_need_lab.txt',np.array(qoi_need_lab))


#Train XGBoost to distinguish between zero and non-zero values
np.random.seed(10)
for q in qoi_need_lab:
  f1_best = -9999
  for i in range(9,100,3):
    clc = XGBClassifier(n_estimators = i, booster = 'gbtree', max_depth = int(np.sqrt(i*2.5)), objective='binary:logistic', verbosity = 0, use_label_encoder = False)

    clc.fit(ptrain_norm, ytrain_norm_label[:, q])
    ypredict_train_label = clc.predict(ptrain_norm)
    ypredict_val_label = clc.predict(pval_norm)

    preci_train = precision_score(ypredict_train_label, ytrain_norm_label[:, q])
    preci_val = precision_score(ypredict_val_label, yval_norm_label[:, q])
    recall_train = recall_score(ypredict_train_label, ytrain_norm_label[:, q])
    recall_val = recall_score(ypredict_val_label, yval_norm_label[:, q])

    print(f'qoi {q}, round {i}, validation, precision: {preci_val}, recall: {recall_val}')

    f1_val = preci_val * recall_val / (preci_val + recall_val) * 2

    if (f1_val > f1_best):
      myfile = open(UQ_output+'/NN_surrogate/fitstats_lab_'+str(q)+'.txt','w')
      myfile.write('Number of parameters:                        '+str(nparms)+'\n')
      myfile.write('Number of outputs:                           1\n')
      myfile.write('Number of training samples:                  '+str(ntrain)+'\n')
      myfile.write('Number of training samples, PFT grow:        '+str(sum(ytrain_norm_label[:,q]))+'\n')
      myfile.write('Number of validation samples:                '+str(nval)+'\n\n')
      myfile.write('Number of validation samples, PFT grow:      '+str(sum(yval_norm_label[:,q]))+'\n')
      myfile.write('Best XGBoost classifier:\n')
      myfile.write('    n_estimators:                            '+str(i)+'\n')
      myfile.write('    max_depth:                               '+str(int(np.sqrt(i*2.5)))+'\n')
      myfile.write('QOI validation '+str(q)+' (precision,recall): '+str(preci_val)+'  '+str(recall_val)+'\n')
      f1_best = f1_val
      pkl_filename = UQ_output+'/NN_surrogate/classify_'+str(q)+'.pkl'
      with open(pkl_filename,'wb') as file:
        pickle.dump(clc, file)
      myfile.close()

    if (f1_best > 0.95):
      print('F1 score > 0.95')
      break


#Train MLP on both zero and non-zero values values
corr_best = -99999
for n in range(0,100):
  nmin = 10*np.sqrt(n+1) #max(10, ntrain/20)
  nmax = 20*np.sqrt(n+1) #min(ntrain/4, 100)
  nl = int(np.random.uniform(nmin,nmax))
  nl2 = int(np.random.uniform(nmin,nmax))*2
  do3 = 0 #np.random.uniform(0,1)
  nl3 = int(np.random.uniform(nmin,nmax))
  if (do3 > 0.5):
    clf = MLPRegressor(solver='adam', early_stopping=True, tol=1e-7, hidden_layer_sizes=(nl,nl2,nl3,), max_iter=200, validation_fraction=0.2)
  else: 
    clf = MLPRegressor(solver='adam', early_stopping=True, tol=1e-7, hidden_layer_sizes=(nl,nl2,), max_iter=200, validation_fraction=0.2)
  clf.fit(ptrain_norm, ytrain_norm[:,qoi_good]) 

  ypredict_train_temp = clf.predict(ptrain_norm)
  ypredict_val_temp   = clf.predict(pval_norm)
  ypredict_train = ytrain_norm.copy()
  ypredict_val = yval_norm.copy()
  ypredict_train[:,qoi_good] = ypredict_train_temp
  ypredict_val[:,qoi_good] = ypredict_val_temp  

  # Try censoring the values
  for qoi in qoi_need_lab:
    pkl_filename = UQ_output+'/NN_surrogate/classify_'+str(qoi)+'.pkl'
    if os.path.exists(pkl_filename):
      with open(pkl_filename, 'rb') as file:
        clc = pickle.load(file)
      ypredict_train_label = clc.predict(ptrain_norm)
      ypredict_val_label = clc.predict(pval_norm)
      ypredict_train[ypredict_train_label,qoi] = 0.
      ypredict_val[ypredict_val_label,qoi] = 0.

  corr_train=[]
  rmse_train = []
  corr_val=[]
  rmse_val=[]
  for qoi in qoi_good:
    corr_train.append((np.corrcoef(ytrain_norm.astype(float)[:,qoi], ypredict_train.astype(float)[:,qoi])[0,1])**2)
    rmse_train.append((sum((ypredict_train[:,qoi]-ytrain_norm[:,qoi])**2)/ntrain)**0.5)
    corr_val.append((np.corrcoef(yval.astype(float)[:,qoi], ypredict_val.astype(float)[:,qoi])[0,1])**2)
    rmse_val.append((sum((ypredict_val[:,qoi]-yval_norm[:,qoi])**2)/nval)**0.5)
  print(corr_val)
  print(rmse_val)

  if (sum(corr_val) > corr_best):
    myfile = open(UQ_output+'/NN_surrogate/fitstats.txt','w')
    myfile.write('Number of parameters:         '+str(nparms)+'\n')
    myfile.write('Number of outputs:            '+str(nqoi)+'\n')
    myfile.write('Number of good outputs:       '+str(len(qoi_good))+'\n')
    myfile.write('Number of training samples:   '+str(ntrain)+'\n')
    myfile.write('Number of validation samples: '+str(nval)+'\n\n')
    myfile.write('Best neural network:\n')
    myfile.write('Size of NN layer 1: '+str(nl)+'\n')
    myfile.write('Size of NN layer 2: '+str(nl2)+'\n')
    if (do3 == 1):
      myfile.write('Size of NN layer 3: '+str(nl3)+'\n')
    for q in range(0, len(qoi_good)):
      myfile.write('QOI validation '+str(qoi_good[q])+' (R2,rmse): '+str(corr_val[q])+'  '+ \
             str(rmse_val[q]**2)+'\n')
    corr_best = sum(corr_val)
    pkl_filename = UQ_output+'/NN_surrogate/NNmodel.pkl'
    ypredict_val_best = ypredict_val
    with open(pkl_filename,'wb') as file:
      pickle.dump(clf, file)
    for q in qoi_good:
      plt.clf()
      plt.scatter(yval[:,q], ypredict_val_best[:,q]*(yrange[1,q]-yrange[0,q])+yrange[0,q])
      plt.xlabel('Model '+outnames[q])
      plt.ylabel('Surrogate '+outnames[q])
      plt.savefig(UQ_output+'/NN_surrogate/nnfit_qoi'+str(q)+'.pdf')
    myfile.close()
  if (min(corr_val) > 0.99):
    print('All QOIs have R2 > 0.99')
    break