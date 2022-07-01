#!/usr/bin/python
# -*- coding:utf-8 -*-
import numpy as np
import matplotlib as mpl
mpl.use('Agg') #silent mode
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker

title_set='$GaAs$ band structure'
de=0  # The difference between the Fermi energies of self- and non self-consistant calculation, respectively. nonself - self
ybottom=-6.0  # Energy limit
yupper=10.0   # Energy limit
label_set="DFT"
group_labels=[]
x=[]

with open("KLABELS", 'r') as reader:
    lines=reader.readlines()[1:]
for i in lines:
    s=i.encode('utf-8')#.decode('latin-1')
    if len(s.split())==2 and not s.decode('utf-8','ignore').startswith("*"):
        group_labels.append(s.decode('utf-8','ignore').split()[0])
        x.append(float(s.split()[1]))

#hsp=np.loadtxt(open("KLABELS", encoding='utf8'), dtype=np.string_,skiprows=1,usecols = (0,1))
df=np.loadtxt("BAND.dat",dtype=np.float64)
#wf=np.loadtxt("wannier90_band.dat",dtype=np.float64)
for index in range(len(group_labels)):
    if group_labels[index]=="GAMMA" or group_labels[index]=="G":
        group_labels[index]="Γ"

#except:
    #print("failed to open KLABELS containing High symmetry point!")



axe = plt.subplot(111)
axe.plot(df[:,0],df[:,1:]+de,linewidth=1.2,color='blue',label=label_set)
#axe.scatter(wf[:,0],wf[:,1:]-2.78554,color='red',marker='.',s=2.0,label="Wannier")
axe.legend()
font = {'family' : 'DejaVu Sans',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 20,  
        }  

plt.legend(loc='center right')
axe.set_title(title_set,fontsize=20,fontname='DejaVu Sans')
#axe.set_xlabel(r'$\mathbf{k}$-points',fontsize=13)
axe.set_ylabel(r'${E}$-$E_{f}$ (eV)',fontdict=font)
axe.set_xticks(x)
plt.yticks(fontsize=15,fontname='DejaVu Sans')
#plt.yticks(fontdict=ticksfont)
axe.set_xticklabels(group_labels, rotation=0,fontsize=12,fontname='DejaVu Sans')
axe.axhline(y=0, xmin=0, xmax=1,linestyle= '--',linewidth=1.0,color='red')
for i in x[1:-1]:
	axe.axvline(x=i, ymin=0, ymax=1,linestyle= '-',linewidth=0.6,color='0.5')
axe.set_xlim((x[0], x[-1]))
plt.ylim((ybottom, yupper)) # set y limits manually
fig = plt.gcf()
fig.set_size_inches(8,6)
#plt.show()
plt.savefig('band.png',dpi=1000)
