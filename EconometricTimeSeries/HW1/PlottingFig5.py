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


# labels = ['01.Jan','02.Jan','03.Jan','04.Jan','05.Jan','06.Jan','07.Jan','08.Jan','09.Jan','10.Jan','11.Jan','12.Jan','13.Jan']
labels = ['04.Jan','05.Jan','06.Jan','07.Jan','08.Jan','09.Jan','10.Jan','11.Jan','12.Jan','13.Jan']
fig,(ax1,ax2) = plt.subplots(2,1,sharex=False)
# line, = ax1.plot(cpi[12:],linewidth=2.)
line, = ax1.plot(cpi[36:],linewidth=2.)
line2, = ax2.plot(inf[24:],linewidth=2.)
ax1.xaxis.set_ticks(list(range(9,129,12)))
ax1.set_xticklabels(labels,rotation=-45)
ax2.xaxis.set_ticks(list(range(9,129,12)))
ax2.set_xticklabels(labels,rotation=-45)
ax1.set_ylabel('2010=100',fontsize=10)
ax2.set_ylabel('%',fontsize=10,rotation=0)
ax1.set_title('CPI',fontsize=15)
ax2.set_title('Inflation Rate',fontsize=15)
plt.show()



