import numpy as np
from xgboost import XGBClassifier # MLPClassifier, RandomForestClassifier do not work well on shrubs
from sklearn.neural_network import MLPRegressor
import os
import pickle

class MyModel(object):

    def __init__(self,case=''):

        UQdir = './UQ_output/'+case
        self.ptrain = np.loadtxt(UQdir+'/data/ptrain.dat')
        self.ytrain = np.loadtxt(UQdir+'/data/ytrain.dat')
        self.nparms = self.ptrain.shape[1]
        self.nobs   = self.ytrain.shape[1]
        self.ntrain=self.ptrain.shape[0]
        self.yrange = np.loadtxt(UQdir+'/NN_surrogate/yrange.txt').astype(int)

        self.pmin = np.zeros([self.nparms], float)
        self.pmax= np.zeros([self.nparms], float)
        for i in range(0,self.nparms):
          self.pmin[i] = min(self.ptrain[:,i])
          self.pmax[i] = max(self.ptrain[:,i])

        pnamefile = open(UQdir+'/data/pnames.txt','r')
        self.parm_names = []
        for s in pnamefile:
          self.parm_names.append(s[:-1])
        pnamefile.close()
        print(self.parm_names, self.nobs, self.nparms)
        self.pdef = np.zeros([self.nparms], float)
        for i in range(0,self.nparms):
          self.pdef[i] = (self.pmin[i]+self.pmax[i])/2
        self.obs = np.zeros([self.nobs],float)
        self.obs_err = np.zeros([self.nobs],float)
        self.obs_name = []
        obsfile = open(UQdir+'/data/obs.dat')
        i=0
        for s in obsfile:
          self.obs[i] = float(s.split()[0])
          self.obs_err[i] = float(s.split()[1])
          i=i+1
        obsfile.close()
        outnamesfile = open(UQdir+'/data/outnames.txt')
        for s in outnamesfile:
          self.obs_name.append(s[:-1])
        outnamesfile.close()
        pkl_filename = UQdir+'/NN_surrogate/NNmodel.pkl'
        with open(pkl_filename, 'rb') as file:
          self.nnmodel = pickle.load(file)
        self.qoi_good = np.loadtxt(UQdir+'/NN_surrogate/qoi_good.txt').astype(int)

        self.qoi_need_lab = np.loadtxt(UQdir+'/NN_surrogate/qoi_need_lab.txt').astype(int)
        self.classify_model = {}
        for qoi in self.qoi_need_lab:
          pkl_filename = UQdir+'/NN_surrogate/classify_'+str(qoi)+'.pkl'
          with open(pkl_filename, 'rb') as file:
            self.classify_model[qoi] = pickle.load(file)

    def run(self,parms):
        if ((parms).ndim == 1):
            nsamples=1
            theseparms=parms
        else:
            nsamples=parms.shape[0]
        parms_nn = np.zeros([nsamples,self.nparms],float)
        for n in range(0,nsamples):
          if nsamples > 1:
              theseparms = parms[n,:]
          for p in range(0,self.nparms):
            parms_nn[n,p] = (theseparms[p]-self.pmin[p])/(self.pmax[p]-self.pmin[p])
        self.output = np.zeros([nsamples,self.nobs])
        output_val_temp = self.nnmodel.predict(parms_nn)

        qgood=0
        for q in range(0,self.nobs):
          if (q in self.qoi_good):
            self.output[:,q] = output_val_temp[:,qgood]*(self.yrange[1,q]-self.yrange[0,q])+self.yrange[0,q]
            qgood=qgood+1
            if (q in self.qoi_need_lab):
              output_iszero = self.classify_model[q].predict(parms_nn)
              output_val_temp[output_iszero, q] = -9999 # put an invalid value to make sure the PFT grows
          else:
            self.output[:,q] = self.yrange[1,q]
