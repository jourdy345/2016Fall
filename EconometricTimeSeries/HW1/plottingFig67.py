import numpy as np
import matplotlib.pyplot as plt
import matplotlib.text as txt
import types

def rate(y,lag):
  T = len(y)
  g = (np.log(y[lag:T])-np.log(y[0:(T-lag)]))*100.
  return g

data = np.loadtxt('/Users/daeyounglim/Downloads/data.txt',skiprows=1)
date = data[:,0]
base = data[:,1]
cpi  = data[:,2]
inf  = rate(cpi,12)
imp_price = data[:,3]
dimp_price = rate(imp_price,12)
ip = data[:,5]
dip = rate(ip,12)

labels = ['01.Jan','02.Jan','03.Jan','04.Jan','05.Jan','06.Jan','07.Jan','08.Jan','09.Jan','10.Jan','11.Jan','12.Jan','13.Jan']
plt.figure(figsize=(2,2))
ax1 = plt.subplot(221)
ax2 = plt.subplot(222)
ax3 = plt.subplot(223)
ax4 = plt.subplot(224)
# fig,(ax1,ax2,ax3,ax4) = plt.subplots(2,2,sharex=False)
line1, = ax1.plot(ip[12:],'b',linewidth=2.)
line2, = ax3.plot(dip,'b',linewidth=2.)
line3, = ax2.plot(imp_price[12:],'k',linewidth=2.)
line4, = ax4.plot(dimp_price,'k',linewidth=2.)
ax1.set_title('Industrial Production(SA)')
ax2.set_title('Petroleum Products Index')
ax3.set_title('Growth Rate, Industrial Production')
ax4.set_title('The Growth Rate of Import Price Index for Petroleum Products')
ax1.set_ylabel('2010=100')
ax2.set_ylabel('%',rotation=-45)
ax1.set_xticks(list(range(9,165,12)),labels)
ax1.set_xticklabels(labels,rotation=-45)
ax1.grid(True)
ax2.set_xticks(list(range(9,165,12)),labels)
ax2.set_xticklabels(labels,rotation=-45)
ax2.grid(True)
ax3.set_xticks(list(range(9,165,12)),labels)
ax3.set_xticklabels(labels,rotation=-45)
ax3.grid(True)
ax4.set_xticks(list(range(9,165,12)),labels)
ax4.set_xticklabels(labels,rotation=-45)
ax4.grid(True)

plt.show()