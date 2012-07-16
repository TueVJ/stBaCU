import numpy
from mismatch import *
from mpolicies import *
from wavelet_smoother import wavelet_smoother
###
#	Returns a list of tuples (begin,end) with the beginning
#	and end indices of intervals of length > 1 where the
#	array Ch is non-zero.
#	Returns an empty list if none such intervals were found.
###
def get_nonzero_periods(Ch,tolerance=10.**-8):
	n=len(Ch)
	nonzero=[0,]*n
	out=[]
	start=0
	toggle=False
	for i in range(n):
		nonzero[i]=int(abs(Ch[i])>tolerance)
	for i in range(n-1):
		if nonzero[i]==1:
			if toggle==False and nonzero[i+1]==1:
				#We are beginning a streak!
				toggle=True
				start=i
			elif toggle==True and nonzero[i+1]==0:
				#We are ending a streak.
				toggle=False
				out.append((start,i))
	if toggle == True:
		#We ended on a streak.
		out.append((start,n-1))
	return out

##
#	Returns an array of the maximum values over
#	each nonzero interval in Ch
##
def get_max_of_nonzero(Ch,return_indices=False,mintolerance=10.**-8):
	indices=get_nonzero_periods(Ch,tolerance=mintolerance)
	absCh=abs(Ch)
	out=[]
	max_indices=[]
	for indextuple in indices:
		out.append(max(absCh[indextuple[0]:indextuple[1]]))
		max_indices.append(indextuple[0]+numpy.argmax(Ch[indextuple[0]:indextuple[1]]))
	if return_indices:
		return out,max_indices
	else:
		return out

##
#	Returns the intermittent maxima of the storage level
#	at the chosen value of renewable penetration gamma,
#	and the wind component of the mix, a.
#	Note, that choosing gamma close to 1 will give very
#	few maxima.
##
def get_maximal_storages(gamma=0.7,a=0.6):
	ts=get_mismatch(gamma,a)
	dummy,storage_capacity,storage_level=get_policy_2_storage(ts,return_storage_filling_time_series=True)
	if gamma>1:
		storage_level=storage_level-storage_capacity
	maxes,ind_max=get_max_of_nonzero(storage_level,return_indices=True)
	return maxes
##
#	Same as above, but with an intermittent storage of size
#	storage_size. Policy is to always use the intermittent
#	storage first.
##
def get_maximal_storages_with_extra_storage(gamma=0.7,storage_size=6,tolerance=10.**-11):
	ts=get_mismatch(gamma)
	ts6h,dummy=get_policy_2_storage(ts,storage_capacity=storage_size)
	dummy,storage_capacity,storage_level=get_policy_2_storage(ts6h,return_storage_filling_time_series=True)
	if gamma>1:
		storage_level=storage_level-storage_capacity
	plot(storage_level)
	maxes,ind_max=get_max_of_nonzero(storage_level,return_indices=True,mintolerance=tolerance)
	return maxes

def get_scale_limits(scales):
	output=[]
	for low in scales:
		for high in scales:
			if low < high:
				output.append((low,high))
	return output

# Saves the smoothed slts defined through gamma and scales.
def save_smoothed_slts():
	gammas=np.arange(0.5,2.0,0.005)
	
	#Interesting scales are from 1 hour to 100000 hours (Covers entire dataset)
#	scales=np.logspace(0.0,5.0,16)
#	scalelimits=get_scale_limits(scales)
#	print(scalelimits)
	for gamma in gammas:
		ts=get_mismatch(gamma)
		dummy,dummytoo,slts=get_policy_2_storage(ts,return_storage_filling_time_series=True)
		np.save('slts/us_None_None_'+ \
		str(gamma),slts)
		'''
		for limitvec in scalelimits:
			smooth=wavelet_smoother(ts,\
			long_cut=limitvec[1],short_cut=limitvec[0])
			try:
				dummy,dummytoo,slts=get_policy_2_storage(smooth,return_storage_filling_time_series=True)
			except ValueError:
				print('Empty storage for '+str(gamma)+' '+str(limitvec[0])+' '+str(limitvec[1]))
			except TypeError:
				print('Crossing problem for '+str(gamma)+' '+str(limitvec[0])+' '+str(limitvec[1]))
			else:
				np.save('slts/wavelet_'+str(limitvec[0])+\
				'_'+str(limitvec[1])+'_'+\
				str(gamma),slts)
		'''
	return 0


