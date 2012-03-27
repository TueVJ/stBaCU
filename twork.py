import init.py
import numpy as np
import pp
import sys

def get_scale_limits(scales):
	output={}
	for low in scales:
		for high in scales:
			if low < high:
				output.add((low,high))
	return output

N=101
gammasub=np.linspace(0.5,1.0,N)
gammasup=2*gammasub
gammasub=gammasub.tolist()
dummy=gammasub.pop()
gammasub=np.array(gammasub)
gammas=np.concatenate(gammasub,gammasup)

#Interesting scales are from 1 hour to 100000 hours (Covers entire dataset)
scales=np.logspace(0,5,25)
scalelimits=get_scale_limits(scales)

'''
#pp setup
ppservers = ()
if len(sys.argv) > 1:
    ncpus = int(sys.argv[1])
    # Creates jobserver with ncpus workers
    job_server = pp.Server(ncpus, ppservers=ppservers)
else:
    # Creates jobserver with automatically detected number of workers
    job_server = pp.Server(ppservers=ppservers)

print "Starting pp with", job_server.get_ncpus(), "workers"
'''
tsdict={}
for gamma in gammas:
	tsdict[gamma]=get_mismatch(gamma)

for gamma in gammas:
	save('slts/us_low=None_high=None_gamma='+...
	str(gamma)+'.dat',tsdict[gamma])
