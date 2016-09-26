import matplotlib.pyplot as plt
import numpy as np
def simulateLinearRegression(T,p,burnin,thinin,nSample,drawPlot=True):
  if T < p: raise Exception('Singular design matrix: p should be smaller than T')
  X = 5.*np.random.rand(T,p)
  b = np.arange(1.,float(p)+1.)
  sig2 = 1.5
  Y = X.dot(b)+np.sqrt(2.)*np.random.randn(T)

  b0 = np.zeros(p)
  B0 = np.eye(p)
  a0 = 10.
  d0 = 10.

  sig2j = d0/a0
  bm = np.zeros((nSample,p))
  sig2m = np.zeros(nSample)
  XX = X.transpose().dot(X)
  XY = X.transpose().dot(Y)
  B0inv = np.linalg.inv(B0)
  a1 = a0+float(T) # a0 and T does not change, inefficient if recomputed in every iteration

  for _ in range(burnin):
    B1inv = (1./sig2j)*XX+B0inv
    B1    = np.linalg.inv(B1inv)
    A     = (1./sig2j)*XY+B0inv.dot(b0)
    beta1 = B1.dot(A)
    betaj = beta1+np.linalg.cholesky(B1).transpose().dot(np.random.randn(p))

    ehat  = Y-X.dot(betaj)
    d1    = d0+(ehat*ehat).sum()
    sig2j = 1./np.random.gamma(shape=a1/2.,scale=2./d1,size=1)[0]
  for iter in range(nSample):
    for _ in range(thinin):
      B1inv = (1./sig2j)*XX+B0inv
      B1    = np.linalg.inv(B1inv)
      A     = (1./sig2j)*XY+B0inv.dot(b0)
      beta1 = B1.dot(A)
      betaj = beta1+np.linalg.cholesky(B1).transpose().dot(np.random.randn(p))
      
      ehat = Y-X.dot(betaj)
      d1 = d0+(ehat*ehat).sum()
      sig2j = 1./np.random.gamma(shape=a1/2.,scale=2./d1,size=1)[0]

    bm[iter,:] = betaj
    sig2m[iter] = sig2j
  if drawPlot:
    plt.close('all')
    fig = plt.figure(1)
    fig.suptitle('Histograms of posterior distributions')
    ax = plt.subplot(131)
    ax.set_title(r'Distribution of $\beta_{1}$')
    ax.hist(bm[:,0])
    ax = plt.subplot(132)
    ax.set_title(r'Distribution of $\beta_{2}$')
    ax.hist(bm[:,1])
    ax = plt.subplot(133)
    ax.set_title(r'Distribution of $\sigma^{2}$')
    ax.hist(sig2m)
    plt.show()
    print('We only provide plots up to beta2')
    return {'beta':bm,'sigma':sig2m,'betaPosteriorMean':bm.mean(axis=0),'sigmaPosteriorMean':sig2m.mean()}
  else:
    return {'beta':bm,'sigma':sig2m,'betaPosteriorMean':bm.mean(axis=0),'sigmaPosteriorMean':sig2m.mean()}

res = simulateLinearRegression(500,30,10000,20,5000)