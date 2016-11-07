import numpy as np
import matplotlib.pyplot as plt


def ppl_ppd(yfull,Xfull,mu_beta0,Sigma_beta0,r0,s0,burnin,thinin,nSample,H):
  invSigma_beta0 = np.linalg.inv(Sigma_beta0)
  invSigma_muBeta = invSigma_beta0.dot(mu_beta0)
  (T,p) = Xfull.shape
  r_q = 0.5*(float(T)+r0)
  # initialize sigma
  sigma = 0.5*r0/(0.5*s0-1.)
  if H == 0.:
    X = Xfull
    y = yfull
    XtX = X.T.dot(X)
    Xty = X.T.dot(y)

    for _ in range(burnin):
      # sample beta
      Sigma_q = np.linalg.inv(1./sigma*XtX+invSigma_beta0)
      mu_q    = Sigma_q.dot(1./sigma*Xty+invSigma_muBeta)
      beta = mu_q+np.linalg.cholesky(Sigma_q).T.dot(np.random.randn(p))


      # sample sigma
      y_Xbeta = y-X.dot(beta)
      s_q = 2./(s0+np.sum(y_Xbeta*y_Xbeta))
      sigma = 1./np.random.gamma(shape=r_q,scale=s_q,size=1)[0]
    beta_container = np.zeros((nSample,p))
    sigma_container = np.zeros(nSample)
    for i in range(nSample):
      for _ in range(thinin):
        # sample beta
        Sigma_q = np.linalg.inv(1./sigma*XtX+invSigma_beta0)
        mu_q    = Sigma_q.dot(1./sigma*Xty+invSigma_muBeta)
        beta = mu_q+np.linalg.cholesky(Sigma_q).T.dot(np.random.randn(p))

        # sample sigma
        y_Xbeta = y-X.dot(beta)
        s_q = 2./(s0+np.sum(y_Xbeta*y_Xbeta))
        sigma = 1./np.random.gamma(shape=r_q,scale=s_q,size=1)[0]
      beta_container[i,:] = beta
      sigma_container[i] = sigma
    return {'beta':beta_container,'sigma':sigma_container,'PPL':0.,'lnPPL':0.}
  else:
    PPLm = np.zeros(H)
    T0   = yfull.shape[0]-H
    n = burnin+nSample
    beta_container = np.zeros((n,p))
    sigma_container = np.zeros(n)
    for indH in range(H):
      y = yfull[:(T0+indH-2)]
      X = Xfull[:(T0+indH-2)]
      yf = yfull[T0+indH-1]
      xf = Xfull[(T0+indH-1),:]
      PPL_Hm = np.zeros(n)
      for it in range(n):
        k = X.shape[1]
        XtX = X.T.dot(X)
        Xty = X.T.dot(y)
        Sigma_q = np.linalg.inv(1./sigma*XtX+invSigma_beta0)
        mu_q    = Sigma_q.dot(1./sigma*Xty+invSigma_muBeta)
        beta = mu_q+np.linalg.cholesky(Sigma_q).T.dot(np.random.randn(k))
        beta_container[it,:] = beta

        y_Xbeta = y-X.dot(beta)
        s_q = 2./(s0+np.sum(y_Xbeta*y_Xbeta))
        sigma = 1./np.random.gamma(shape=r_q,scale=s_q,size=1)[0]
        sigma_container[it] = sigma

        c = -0.5*np.log(2.*sigma*np.pi)
        e = yf-xf.dot(beta)
        e2 = e*e
        PPL_Hm[it] = c-0.5*e2/sigma
      PPLm[indH] = np.log(np.mean(np.exp(PPL_Hm[burnin:])))
    lnPPL = np.sum(PPLm)
    return {'beta':beta_container,'sigma':sigma_container,'PPLm':PPLm,'lnPPL':lnPPL}


b = np.array([1.,2.,3.,4.])
X = 5.*np.random.randn(100,4)
y = X.dot(b)+np.random.randn(100)*np.sqrt(2.)
mu_beta0 = np.zeros(4)
Sigma_beta0 = np.eye(4)
r0 = 10.
s0 = 10.
burnin = 10000
thinin = 20
nSample = 5000
lin2 = ppl_ppd(y,X,mu_beta0,Sigma_beta0,r0,s0,burnin,thinin,nSample,5)


plt.hist(lin['sigma'])
plt.show()