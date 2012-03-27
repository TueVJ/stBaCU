import numpy
from mfunc import *
from mpolicies import *
import pywt
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

##
#	Estimates the Extreme Value Theory index of the
#	timeseries X using the Pickands estimator.
##
def gamma_estimator_Pickands(X,return_k=False,excluder=0):
	X.sort() #Sorts list by smallest first
	X=X[:len(X)-excluder]
	N=len(X)/4
	if N<2+excluder: #Need at least two k values to produce proper output.
		print "Insufficient values in Pickands estimator"
		if return_k:
			return [0,0],[1,2] 
		else:
			return [0,0]
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
	if len(X)<=2+excluder:
		print "Insufficient values in adapted Hill estimator"
		if return_k:
			return [0,0],[1,2] 
		else:
			return [0,0]
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


def wavelet_smoother(ts,long_cut=np.inf,short_cut=0,wavelet='db4',normalize=True):
	''' Smooth ts using the wavelet specified.

		@long_cut is the longest wavelength to keep.
		@short_cut is the shortest wavelength to keep.
		@normalize will normalize positive and negative
		contributions after smoothing.
	'''
	smooth_long=smooth_short=True
	if long_cut >= len(ts):
		smooth_long=False
	if short_cut < 1:
		smooth_short=False
	
	# No smoothing needed
	if smooth_long == False and smooth_short == False:
		return ts
	
	# Oversmoothed; Output 0.
	if long_cut <= short_cut or short_cut >= len(ts):
		return ts*0
	if normalize==True:
		avg=np.average(ts)
		ts-=avg
	coeffs=pywt.wavedec(ts,wavelet)
	
	#indices of short and long wavelengths.
	if smooth_short==True:
		si=np.log2((1.0*len(ts))/short_cut)
		if si >= len(coeffs):
			smooth_short=False
	if smooth_long == True:
		li=np.log2((1.0*len(ts))/long_cut)
	#No information on short wavelength

	if smooth_long==True:
		intpart=int(li)
		fracpart=li-intpart
		for i in range(intpart):
			coeffs[i]*=0
		coeffs[intpart]*=fracpart

	if smooth_short==True:
		intpart=int(si)
		fracpart=si-intpart
		for i in range(intpart+1,len(coeffs)):
			coeffs[i]*=0
		coeffs[intpart]*=fracpart
	
	smoothed=pywt.waverec(coeffs,wavelet)

	# Normalization step, ensures the positive and negative
	# parts of smoothed match those of ts. "Preserves gamma"
	if normalize==True:
		smoothed+=avg
		ts+=avg
	return smoothed
