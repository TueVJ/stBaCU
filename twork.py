from init import *
import numpy as np
import os
import random
from mpl_toolkits.mplot3d import Axes3D


dict={}
peaklist=[]
path='slts/'
dirList=os.listdir(path)
count=0.0
tot=len(dirList)*1.0
for fname in dirList:
	ID=fname.rstrip('.npy')
	ts=np.load(path+fname)
	maxim=np.max(ts)
	parts=ID.split('_')
	try:
		ID=parts[0]+' '+str(round(float(parts[1]),2))+\
	' '+str(round(float(parts[2]),2))
	except ValueError:
		ID=parts[0]+' '+parts[1]+\
	' '+parts[2]
		
	gamma=parts[3]
	try:
		dict[ID][0].append(gamma)
		dict[ID][1].append(maxim)
	except KeyError:
		dict[ID]=[[],[]]
		dict[ID][0].append(gamma)
		dict[ID][1].append(maxim)
	if parts[0]=='wavelet' and gamma=='1.0':
		peaklist.append([float(parts[1]),float(parts[2]),maxim])
	count+=1
	print('Loaded '+ str(round(count/tot*100,2)) + '%')

interesting=random.sample(dict.keys(),10)
#interesting=[]
#for ID in dict.keys():
#	if ID.count('us')>0:
#		interesting.append(ID)
linestyles=['-','--','-.']
i=0
fig1=figure()
for ID in interesting:
	plot(dict[ID][0],dict[ID][1],label=ID,\
	linestyle=linestyles[i%len(linestyles)])
	i+=1
yscale('log')
#legend()
interactive(1)
show()

peaklist=np.array(peaklist)
peaklist=peaklist.transpose()
xx=peaklist[0]
yy=peaklist[1]
zz=peaklist[2]
fig=figure()
ax=Axes3D(fig)
ax.scatter3D(np.log10(xx),np.log10(yy),np.log10(zz))
ax.set_xlabel('log(low)')
ax.set_ylabel('log(high)')
ax.set_zlabel('log(Peak height)')
show()
