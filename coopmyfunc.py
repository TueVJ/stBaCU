import numpy
###
#	Returns a list of tuples (begin,end) with the beginning
#	and end indices of intervals of length > 1 where the
#	array Ch is non-negative.
#	Returns an empty list if none such intervals were found.
###
def get_nonzero_periods(Ch):
	n=len(Ch)
	nonzero=[0,]*n
	out=[]
	start=0
	toggle=False
	for i in range(n):
		nonzero[i]=int(abs(Ch[i])>10.**-10)
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

def get_max_of_nonzero(Ch,return_indices=False):
	indices=get_nonzero_periods(Ch)
	out=[]
	max_indices=[]
	for indextuple in indices:
		out.append(max(Ch[indextuple[0]:indextuple[1]]))
		max_indices.append(indextuple[0]+numpy.argmax(Ch[indextuple[0]:indextuple[1]]))
	if return_indices:
		return out,max_indices
	else:
		return out

print get_nonzero_periods([0,1,0,1,1,0])
listofrand=2*numpy.random.rand(500)*(numpy.random.rand(500)<0.9).astype(int)
boollistofrand=2*(listofrand<0.9).astype(int)
print get_nonzero_periods(listofrand)
print get_max_of_nonzero(listofrand)
