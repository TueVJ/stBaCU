import os
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

peaklist=[]
inpath='slts/'
outpath='hist/'
dirList=os.listdir(inpath)
count=0.0
tot=len(dirList)*1.0
tslist=[]
for fname in dirList:
	ID=fname.rstrip('.npy')
	ts=np.load(inpath+fname)
	parts=ID.split('_')
	if parts[0]=='us':
#		if float(parts[3])>1.0:
#			ts -= max(ts)
		tslist.append((float(parts[3]),ts[:]))
#		np.savez(outpath+ID,bins=bins,vals=vals)
	count+=1
	print('Handled '+ str(round(count/tot*100,2)) + '%')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

i=0
colors=['r','g','b','c','m','y','k']
for gamma,ts in tslist:
	vals,bins=np.histogram(ts,bins=50)
	vals=np.log(vals+1)-np.log(2)
	wide=bins[1]-bins[0]
	ax.bar(bins[:-1],vals,zs=gamma,zdir='y',width=0.8*wide,color=colors[i%len(colors)],alpha=0.8)
	i+=1
