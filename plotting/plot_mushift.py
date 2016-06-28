#! /usr/bin/python

import numpy as np
import math
import matplotlib.pyplot as plt


from pylab import rcParams
#rcParams['figure.figsize'] = 10, 8
plt.rc('font',family='Times New Roman')
plt.rc('font', size=18)
    



vecs=[np.loadtxt("../output/mu2.5_fields.txt"),np.loadtxt("../output/mu5_fields.txt"),np.loadtxt("../output/mu7.5_fields.txt")]

radii=[]
sample=[]

for field in vecs:
	radii.append(field[:,0])
	sample.append(field[:,1])

def normgaus(x,mu):
	ret = (1/(math.sqrt(2*3.141592)))*math.exp(-0.5*(x-mu)**2)
	return ret
vec_normgaus = np.vectorize(normgaus)

kspace = np.linspace(0,10,100)
#fig,ax=plt.subplots()       



axes=[]

axes.append(plt.subplot2grid((3,2),(0,0),rowspan=1,colspan=1))
axes[0].plot(radii[0],sample[0],"o-",color="red",label=r"$\mu=2.5$",lw=2)
#axes[0].set_xlabel("Radius (A.U.)")
axes[0].set_ylabel("Field Value (A.U.)")
axes[0].set_ylim([-20,40])

axes.append(plt.subplot2grid((3,2),(0,1),rowspan=1,colspan=1))
axes[1].plot(kspace,vec_normgaus(kspace,2.5),"-",color="red",label=r"$\mu=2.5$",lw=2)
#axes[1].set_xlabel("K (A.U.)")
axes[1].set_ylabel("Power (A.U.)")
axes[1].set_ylim([0,0.75])


axes.append(plt.subplot2grid((3,2),(1,0),rowspan=1,colspan=1))
axes[2].plot(radii[1],sample[1],"o-",color="green",label=r"$\mu=2.5$",lw=2)
#axes[2].set_xlabel("Radius (A.U.)")
axes[2].set_ylabel("Field Value (A.U.)")
axes[2].set_ylim([-20,40])

axes.append(plt.subplot2grid((3,2),(1,1),rowspan=1,colspan=1))
axes[3].plot(kspace,vec_normgaus(kspace,5.0),"-",color="green",label=r"$\mu=5.0$",lw=2)
#axes[3].set_xlabel("K (A.U.)")
axes[3].set_ylabel("Power (A.U.)")
axes[3].set_ylim([0,0.75])

axes.append(plt.subplot2grid((3,2),(2,0),rowspan=1,colspan=1))
axes[4].plot(radii[2],sample[2],"o-",color="blue",label=r"$\mu=2.5$",lw=2)
axes[4].set_xlabel("Radius (A.U.)")
axes[4].set_ylabel("Field Value (A.U.)")
axes[4].set_ylim([-20,40])

axes.append(plt.subplot2grid((3,2),(2,1),rowspan=1,colspan=1))
axes[5].plot(kspace,vec_normgaus(kspace,7.5),"-",color="blue",label=r"$\mu=7.5$",lw=2)
axes[5].set_xlabel("K (A.U.)")
axes[5].set_ylabel("Power (A.U.)")
axes[5].set_ylim([0,0.75])

plt.suptitle("Radial Field Component: "+r"$\phi_{0,0}$",fontsize=24)
legend1 = axes[1].legend(loc='upper right',fancybox=True, shadow=False)
legend2 = axes[3].legend(loc='upper right',fancybox=True, shadow=False)
legend3 = axes[5].legend(loc='upper right',fancybox=True, shadow=False)
"""
axes.append(plt.subplot2grid((4,4),(0,0),colspan=4,sharex=axes[0]))
axes[1].plot(np.multiply(valvec,param.GeV/param.TeV),np.divide(np.subtract(pert_vec0,num_vec),num_vec),color="red",lw=2,ls='--')
axes[1].plot(np.multiply(valvec,param.GeV/param.TeV),np.divide(np.subtract(pert_vec1,num_vec),num_vec),color="green",lw=2,ls='--')
axes[1].plot(np.multiply(valvec,param.GeV/param.TeV),np.divide(np.subtract(pert_vec2,num_vec),num_vec),color="blue",lw=2,ls='--')
axes[1].plot(np.multiply(valvec,param.GeV/param.TeV),np.divide(np.subtract(num_vec,num_vec),num_vec),color="black",lw=2,ls='-')
#axes[1].set_yscale('log')
axes[1].set_xlim([np.multiply(valvec,param.GeV/param.TeV)[0],np.multiply(valvec,param.GeV/param.TeV)[-1]])
axes[1].set_ylabel("Relative Residual")
axes[1].yaxis.set_major_formatter(mpl.ticker.ScalarFormatter(useMathText=True, useOffset=False))
plt.setp(axes[1].get_xticklabels(), visible=False)



ax.plot(radii[0],sample[0],"o-",color="red",label=r"$\mu=2.5$")
ax.plot(radii[1],sample[1],"o-",color="green",label=r"$\mu=5.0$")
ax.plot(radii[2],sample[2],"o-",color="blue",label=r"$\mu=7.5$")

#ax.set_yscale('log')
#ax.set_ylim([1e-1,1])
#ax.set_ylim([0,1.0])
#ax.set_xlim([0,1.65e4])
ax.set_xlabel("Radius (A.U.)")
ax.set_ylabel("Field Value (A.U.)")
#ax.set_title("Antineutrino With Decay, With Matter")


# Now add the legend with some customizations.
legend = ax.legend(loc='center right', shadow=False)
legend = ax.legend(bbox_to_anchor=(0.9, 0.6),fancybox=True, bbox_transform=plt.gcf().transFigure, shadow=False)
"""
"""
legend = ax.legend(loc='center right', shadow=False)

# The frame is matplotlib.patches.Rectangle instance surrounding the legend.
frame = legend.get_frame()
frame.set_facecolor('1.0')

# Set the fontsize
for label in legend.get_texts():
    label.set_fontsize(24)

for label in legend.get_lines():
    label.set_linewidth(1.5)  # the legend line width
"""
"""

ax.tick_params(axis='x', which='major', labelsize=16)
ax.tick_params(axis='y', which='major', labelsize=16)
"""
plt.show()
#print sp.AtmosphericNeutrinoOscillationProbability(1,1,100*param.GeV,param.PI,param)

