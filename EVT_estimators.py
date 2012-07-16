import numpy as np

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
#	Note that this only works for EV index > 0.
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
#	This works for any EV index.
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

