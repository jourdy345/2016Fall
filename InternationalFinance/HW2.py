import numpy as np
from scipy.stats import multivariate_normal
import matplotlib.pyplot as plt

def genBeta(Y,X,b0,B0inv,sig2):
  XX = X.T.dot(X)
  XY = X.T.dot(Y)
  B1inv = (1./sig2)*XX+B0inv
  B1    = np.linalg.inv(B1inv)
  A     = (1./sig2)*XY+B0inv.dot(b0)
  beta1 = B1.dot(A)
  return beta1+np.linalg.cholesky(B1).T.dot(np.random.randn(2))

def genSig2(Y,X,a0,d0,beta):
  T = Y.shape[0]
  a1 = a0+T
  ehat = Y-X.dot(beta)
  d1 = (ehat*ehat).sum()+d0
  return np.random.gamma(shape=a1/2.,scale=2./d1,size=1)[0]

def gibbs(burnin,nSample):
  b1 = np.array([1.,2.])
  b2 = np.array([1.,4.])
  T = 100
  tru_tau = 49
  T1 = tru_tau-1
  T2 = T-tru_tau+1
  X1 = np.column_stack((np.ones(T1),np.random.rand(T1)*10.))
  X2 = np.column_stack((np.ones(T2),np.random.rand(T2)*10.))
  #X1 = np.array([np.ones(T1),np.random.rand(T1)*10.]).reshape(2,T1).T
  #X2 = np.array([np.ones(T2),np.random.rand(T2)*10.]).reshape(2,T2).T
  sig21 = 1.5
  sig22 = sig21
  Y1 = X1.dot(b1)+np.random.randn(T1)*np.sqrt(sig21)
  Y2 = X2.dot(b2)+np.random.randn(T2)*np.sqrt(sig22)

  #Y  = np.column_stack((Y1,Y2))
  #X  = np.column_stack((X1,X2))
  Y = np.append(Y1,Y2)
  X = np.vstack((X1,X2))

  k  = X.shape[1]
  b0 = np.zeros(k)
  B0 = 10.*np.eye(k)

  a0 = 10.
  d0 = 10.

  sig21j = d0/a0
  sig22j = sig21j
  tauj = round(0.5*T)

  # n0 = 200
  # n1 = 1000
  n = burnin+nSample
  lb = round(0.2*T)
  ub = round(0.8*T)
  b1m = np.zeros((n,k))
  b2m = np.zeros((n,k))
  sig21m = np.zeros(n)
  sig22m = np.zeros(n)
  taum = np.zeros(n)
  tau_postm = np.zeros(T)

  B0inv = np.linalg.inv(B0)
  for it in range(n):
    Y1 = Y[:tauj-2]
    X1 = X[:tauj-2,:]
    Y2 = Y[tauj-1:]
    X2 = X[tauj-1:,:]
    beta1j = genBeta(Y1,X1,b0,B0inv,sig21j)
    b1m[it,:] = beta1j
    sig21j = genSig2(Y1,X1,a0,d0,beta1j)
    sig21m[it] = sig21j
    beta2j = genBeta(Y2,X2,b0,B0inv,sig22j)
    b2m[it,:] = beta2j
    sig22j = genSig2(Y2,X2,a0,d0,beta2j)
    sig22m[it] = sig22j
    likm = np.zeros(ub-lb-1)

    for t in range(lb+1,ub):
      lik1 = multivariate_normal.logpdf(Y[:(t-1)],X[:(t-1),:].dot(beta1j),sig21j*np.eye(t-1))
      lik2 = multivariate_normal.logpdf(Y[(t-1):],X[(t-1):,:].dot(beta2j),sig22j*np.eye(T-t+1))
      likm[t-lb-1] = np.exp(lik1+lik2)
    C = 1./likm.sum()
    tau_probm = C*likm
    cumsumprob = np.cumsum(tau_probm)
    u = np.random.rand(1)[0]
    for i in range(cumsumprob.size):
      if i == 0:
        if u>0. and u < cumsumprob[0]:
          tauj = i+lb
      else:
        if u > cumsumprob[i-1] and u < cumsumprob[i]:
          tauj = i+lb
    taum[it] = tauj
    if it > n0:
      tau_postm[tauj] = tau_postm[tauj]+1
    print('current iteration is',it+1)
  return {'beta1':b1m,'beta2':b2m,'sigma2_1':sig21m,'sigma2_2':sig22m,'taum':taum,'tau_postm':tau_postm}

res = gibbs(200,1000)

plt.plot(res['tau_postm'],linewidth=2.)
plt.title('tau probability')
plt.ylabel('Mass')
plt.xlabel('t')
plt.show()