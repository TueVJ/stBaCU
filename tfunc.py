import numpy
from mfunc import *
from mpolicies import *
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
		storage_level=storage_capacity-storage_level
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
		storage_level=storage_capacity-storage_level
	maxes,ind_max=get_max_of_nonzero(storage_level,return_indices=True,mintolerance=tolerance)
	return maxes

##
#	Estimates the Extreme Value Theory index of the
#	timeseries X using the Pickands estimator.
##
def gamma_estimator_Pickands(X,return_k=False,excluder=0):
	X.sort() #Sorts list by smallest first
	X=X[:len(X)-excluder]
	N=len(X)/4
	ks=range(1+excluder,N)
	est=[]
	log2=np.log(2)
	for k in ks:
		est.append(np.log((X[-k]-X[-2*k])/(X[-2*k]-X[-4*k]))/log2)
	if return_k==True:
		return est,ks
	else:
		return est

##
#	Estimates the Extreme Value Theory index of the
#	timeseries X using the Hill estimator.
#	Note that this only works for gamma > 0.
##
def gamma_estimator_Hill(X,return_k=False,excluder=0):
	X.sort() #Sorts list by smallest first
	X=X[:len(X)-excluder]
	N=len(X)/2
	ks=range(1,N)
	est=[]
	for k in ks:
		est.append(sum(log(X[-k:]))/k-log(X[-k]))
	if return_k==True:
		return est,ks
	else:
		return est

##
#	Estimates the Extreme Value Theory index of the
#	timeseries X using the adapted Hill estimator.
#	This works for any gamma.
##
def gamma_estimator_adapted_Hill(X,return_k=False,excluder=0):
	X.sort()
	X=X[::-1] #Sorts list by largest first
	X=X[excluder:]
	N=len(X)-1
	Us=np.zeros(N)
	est=np.zeros(N-1)
	for i in range(N):
		Us[i]=X[i+1]*(sum(np.log(X[:i+1]))/(i+1)-log(X[i+1]))
	for k in range(N-1):
		est[k]=sum(np.log(Us[:k+1]))/(k+1)-np.log(Us[k+1])
	if return_k==True:
		return est,np.arange(N-1)+1
	else:
		return est

