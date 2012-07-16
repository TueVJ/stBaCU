import pywt
import numpy as np

def wavelet_smoother(ts,long_cut=np.inf,short_cut=0,wavelet='db4',preserve_gamma=True):
	''' Smooth ts using the wavelet specified.

		@long_cut is the longest wavelength to keep.
		@short_cut is the shortest wavelength to keep.
		@normalize will normalize positive and negative
		contributions after smoothing.
	'''
	smooth_long=True
	smooth_short=True
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
	if preserve_gamma==True:
		avg=np.average(ts)
		rawts=ts-avg
	coeffs=pywt.wavedec(rawts,wavelet)
	
	#indices of short and long wavelengths.
	if smooth_short==True:
		si=np.log2((1.0*len(ts))/short_cut)
		if si >= len(coeffs):
			smooth_short=False
	if smooth_long == True:
		li=np.log2((1.0*len(ts))/long_cut)
		if li >= len(coeffs):
			return ts*0
	if smooth_long==True:
		intpart=int(li)
		fracpart=li-intpart
		for i in range(intpart):
			coeffs[i]*=0
		coeffs[intpart]*=1-fracpart

	if smooth_short==True:
		intpart=int(si)
		fracpart=si-intpart
		for i in range(intpart+1,len(coeffs)):
			coeffs[i]*=0
		coeffs[intpart]*=fracpart
	
	smoothts=pywt.waverec(coeffs,wavelet)

	# Normalization step, ensures the positive and negative
	# parts of smoothed match those of ts. "Preserves gamma"
	if preserve_gamma==True:
		smoothts+=avg
	return np.array(smoothts)


