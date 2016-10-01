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
fig,ax = plt.subplots()
fig.suptitle('Bank of Korea Base Rate')
line, = ax.plot(base[12:],label='Base Rate',linestyle='--')
plt.ylim([1.5,np.max(base[12:])+.5])
#ax = plt.axes()
plt.xticks(list(range(9,165,12)),labels)
# start,end = ax.get_xlim()
# ax.xaxis.set_ticks(np.arange(start,end,(end-start)/13.))
ax.set_xticklabels(labels,rotation=-45)
ax.annotate('Lehman Brothers files for bankruptcy protection',xy=(100.441,5.21404),xytext=(105,5.5),
            arrowprops=dict(facecolor='black',arrowstyle='->'))
ax.annotate('BOK cuts base rate',xy=(100.688,4.99432),xytext=(105,4.8),
            arrowprops=dict(facecolor='black',arrowstyle='->'))
legend = plt.legend(handles=[line],loc=4)

plt.ylabel('%')
for lab in ax.xaxis.get_majorticklabels():
  lab.customShiftValue = 5.
  lab.set_x = types.MethodType(lambda self,x:  txt.Text.set_x(self, x+self.customShiftValue), 
                                    lab)
plt.show()