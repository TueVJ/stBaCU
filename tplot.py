from mpolicies import *
from tfunc import *
from mismatch import *
import matplotlib.pyplot as plt
import copy
import numpy as np
from wavelet_smoother import wavelet_smoother
from EVT_estimators import *

def plot_gamma_1_storage():
	plt.figure()
	ts=get_mismatch(1.0)
	dummy,storage_capacity,storage_level=get_policy_2_storage(ts,return_storage_filling_time_series=True)
	plt.plot(storage_level)
	plt.interactive(1)
	plt.show()

def plot_storage_with_smoothing():
	N=101
	storage_capacities_raw=empty(N)
	storage_capacities_6h=empty(N)
	storage_capacities_24h=empty(N)
	storage_capacities_1w=empty(N)
	storage_capacities_2w=empty(N)
	storage_capacities_1m=empty(N)
	storage_capacities_6m=empty(N)
	storage_capacities_1y=empty(N)
	storage_capacities_2y=empty(N)
	gammas=linspace(0.5,2.5,N)
	for i in arange(N):
		ts=get_mismatch(gammas[i])
		dummy,storage_capacities_raw[i]=get_policy_2_storage(ts)
		dummy,storage_capacities_6h[i]=get_policy_2_storage(get_smoothed(ts,6))
		dummy,storage_capacities_24h[i]=get_policy_2_storage(get_smoothed(ts,24))
		dummy,storage_capacities_1w[i]=get_policy_2_storage(get_smoothed(ts,7*24))
		dummy,storage_capacities_2w[i]=get_policy_2_storage(get_smoothed(ts,14*24))
		dummy,storage_capacities_1m[i]=get_policy_2_storage(get_smoothed(ts,30*24))
		dummy,storage_capacities_6m[i]=get_policy_2_storage(get_smoothed(ts,6*30*24))
		dummy,storage_capacities_1y[i]=get_policy_2_storage(get_smoothed(ts,365*24))
		dummy,storage_capacities_2y[i]=get_policy_2_storage(get_smoothed(ts,2*365*24))
	
	clf()
	plt.plot(gammas,storage_capacities_raw,label="EU raw")
	plt.plot(gammas,storage_capacities_6h,label="EU 6h")
	plt.plot(gammas,storage_capacities_24h,label="EU 24h")
	plt.plot(gammas,storage_capacities_1w,label="EU 1w")
	plt.plot(gammas,storage_capacities_2w,label="EU 2w")
	plt.plot(gammas,storage_capacities_1m,label="EU 1m")
	plt.plot(gammas,storage_capacities_6m,label="EU 6m")
	plt.plot(gammas,storage_capacities_1y,label="EU 1y")
	plt.plot(gammas,storage_capacities_2y,label="EU 2y")
	plt.setp(plt.gca(),yscale="log")
	legend(loc="upper right")
	interactive(1)
	show()

def plot_EVT_weather_Pickands():
	plt.figure()
	gammas=np.linspace(0.5,0.9,5)
	dict={}
	for gamma in gammas:
		maxes=get_maximal_storages(gamma)
		EVI,ks=gamma_estimator_Pickands(maxes,return_k=True)
		dict[gamma]=(ks[:],EVI[:])
	for gamma in gammas:
		maxk=max(dict[gamma][0])*1.0
		relks=[k/maxk for k in dict[gamma][0]]
		plt.plot(relks,dict[gamma][1],label=str(gamma)+", "+str(maxk))
	plt.gca().axhline(0,linestyle=":")
	plt.gca().set_xlabel("k/kmax")
	plt.gca().set_ylabel("$\hat{\gamma}$")
	plt.gca().set_title("Pickands")
	plt.gca().set_ylim((-1,2))
	plt.gca().set_xlim((0,1))
	plt.legend(loc="upper right")
	plt.savefig("weather_pickands_sub.png")
	plt.interactive(1)
	plt.show()

def plot_EVT_weather_Pickands_w_6h():
	plt.figure()
	gammas=np.linspace(0.9,0.94,5)
	dict={}
	for gamma in gammas:
		maxes=get_maximal_storages_with_extra_storage(gamma)
		EVI,ks=gamma_estimator_Pickands(maxes,return_k=True)
		dict[gamma]=(ks[:],EVI[:])
	for gamma in gammas:
		maxk=max(dict[gamma][0])*1.0
		relks=[k/maxk for k in dict[gamma][0]]
		plt.plot(relks,dict[gamma][1],label=str(gamma)+", "+str(maxk))
	plt.gca().axhline(0,linestyle=":")
	plt.gca().set_xlabel("k/kmax")
	plt.gca().set_ylabel("$\hat{\gamma}$")
	plt.gca().set_title("Pickands")
	plt.gca().set_ylim((-1,2))
	plt.gca().set_xlim((0,1))
	plt.legend(loc="upper right")
	plt.savefig("weather_pickands_sub_6h.png")
	plt.interactive(1)
	plt.show()

def plot_EVT_weather_adHill(exclude=10):
	plt.figure()
	gammas=np.linspace(0.5,0.9,5)
	dict={}
	for gamma in gammas:
		maxes=get_maximal_storages(gamma)
		EVI,ks=gamma_estimator_adapted_Hill(maxes,return_k=True,excluder=exclude)
		dict[gamma]=(ks[:],EVI[:])
	for gamma in gammas:
		maxk=max(dict[gamma][0])*1.0
		relks=[k/maxk for k in dict[gamma][0]]
		plt.plot(relks,dict[gamma][1],label=str(gamma)+", "+str(maxk))
	plt.gca().axhline(0,linestyle=":")
	plt.gca().set_xlabel("k/kmax")
	plt.gca().set_ylabel("$\hat{\gamma}$")
	plt.gca().set_title("Adapted Hill")
	plt.gca().set_ylim((-1,2))
	plt.gca().set_xlim((0.001,1))
	plt.xscale("log")
	plt.legend(loc="upper left")
	plt.savefig("weather_adHill_sub.png")
	plt.interactive(1)
	plt.show()

def plot_EVT_weather_adHill_w_6h(exclude=10):
	plt.figure()
	gammas=np.linspace(0.9,0.94,5)
	dict={}
	for gamma in gammas:
		maxes=get_maximal_storages_with_extra_storage(gamma)
		EVI,ks=gamma_estimator_adapted_Hill(maxes,return_k=True,excluder=exclude)
		dict[gamma]=(ks[:],EVI[:])
	for gamma in gammas:
		maxk=max(dict[gamma][0])*1.0
		relks=[k/maxk for k in dict[gamma][0]]
		plt.plot(relks,dict[gamma][1],label=str(gamma)+", "+str(maxk))
	plt.gca().axhline(0,linestyle=":")
	plt.gca().set_xlabel("k/kmax")
	plt.gca().set_ylabel("$\hat{\gamma}$")
	plt.gca().set_title("Adapted Hill")
	plt.gca().set_ylim((-1,2))
	plt.gca().set_xlim((0.001,1))
	plt.xscale("log")
	plt.legend(loc="upper left")
	plt.savefig("weather_adHill_sub_6h.png")
	plt.interactive(1)
	plt.show()

def plot_EVT_weather():
	close("all")
	maxes=get_maximal_storages(0.7)
	gammaP07,Pks07=gamma_estimator_Pickands(maxes,return_k=True)
	gammaH07,Hks07=gamma_estimator_adapted_Hill(maxes,return_k=True,excluder=10)
	maxes=get_maximal_storages(0.8)
	gammaP08,Pks08=gamma_estimator_Pickands(maxes,return_k=True)
	gammaH08,Hks08=gamma_estimator_adapted_Hill(maxes,return_k=True,excluder=10)
	plt.plot(Pks07,gammaP07,label="Pickands, gamma=0.7",)
	plt.plot(Hks07,gammaH07,label="Adapted Hill, gamma=0.7")
	plt.plot(Pks08,gammaP08,label="Pickands, gamma=0.8")
	plt.plot(Hks08,gammaH08,label="Adapted Hill, gamma=0.8")
	plt.title("Storage level maxima")
	plt.legend(loc="lower right")
	plt.xscale("log")
	plt.gca().axhline(0,linestyle=":")
	plt.gca().set_ylim((-2,2))
	plt.gca().set_xlabel("k")
	plt.gca().set_ylabel("$\hat{\gamma}$")
	plt.savefig("weather.png")

def plot_EVT_uniform():
	somerands=np.random.rand(1000)
	unP,unPks=gamma_estimator_Pickands(somerands, return_k=True)
	unH,unHks=gamma_estimator_adapted_Hill(somerands,return_k=True)
	somerands=np.random.rand(1000)
	unP2,unPks2=gamma_estimator_Pickands(somerands, return_k=True)
	unH2,unHks2=gamma_estimator_adapted_Hill(somerands,return_k=True)
	plt.plot(unPks,unP,label="Pickands")
	plt.plot(unHks,unH,label="Adapted Hill")
	plt.plot(unPks2,unP2,label="Pickands #2")
	plt.plot(unHks2,unH2,label="Adapted Hill #2")
	plt.title("Uniform Distribution")
	plt.xscale("log")
	plt.gca().axhline(0,linestyle=":")
	plt.legend(loc="upper left")
	plt.gca().set_xlabel("k")
	plt.gca().set_ylabel("$\hat{\gamma}$")
	plt.savefig("uniform.png")

def plot_smoothed_singularity():
	N=11
	small_storage_size=6
	gammas=linspace(0.4,2.5,N)
	smoothers=[1,6,24,7*24,14*24]
	smoothernames=["no smooth","6h","24h","1w","2w"]
	num_smoothers=5#len(smoothers)
	
	storage_capacities_raw=[[0 for i in range(N)] for j in range(num_smoothers)]
	storage_capacities_small=[[0 for i in range(N)] for j in range(num_smoothers)]
	for j in arange(num_smoothers):
		for i in arange(N):
			ts=get_mismatch(gammas[i])
			smoothed_ts=get_smoothed(ts ,smoothers[j])
			dummy,storage_capacities_raw[j][i]=get_policy_2_storage(smoothed_ts)
			after_small,dummy=get_policy_2_storage(ts,storage_capacity=small_storage_size)
			dummy,storage_capacities_small[j][i]=get_policy_2_storage(get_smoothed(after_small,smoothers[j]))
	
	fig=figure()
	ax=subplot(111)
	for j in arange(num_smoothers):
		plt.plot(gammas,storage_capacities_raw[j],label="No storage, "+smoothernames[j])
		plt.plot(gammas,storage_capacities_small[j],label=str(small_storage_size)+"h storage, "+smoothernames[j])
	plt.gca().set_ylim(bottom=1)
	yscale("log")
	box = ax.get_position()
	legend(bbox_to_anchor=(0.0,-0.2,1.,0.1),ncol=3, borderaxespad=0.,prop={'size':10},fancybox=True)
	ax.set_position([box.x0, box.y0 + box.height * 0.1,\
	                 box.width, box.height * 0.9])
	interactive(1)
	show()

def plot_wavelets_smoothed_capacities():
	wavelet='db6'
	gammas=np.linspace(1.0,2.0,11)
	scales=np.logspace(1,np.log10(24*365*9),100)
	results={}
	for gamma in gammas:
		ts = get_mismatch(gamma)
		maxes=[]
		for scale in scales:
			print str(scale) + ', '+str(gamma)
			sts=wavelet_smoother(ts,long_cut=scale,wavelet=wavelet)
			try:
				dummy,storage_capacity=get_policy_2_storage(sts)
			except TypeError:
				storage_capacity=0
			maxes.append([scale,storage_capacity])
		results[gamma]=np.array(maxes)
	fig=figure()
	ax=subplot(111)
	for gamma in gammas:
		plot(results[gamma][:,0],results[gamma][:,1],label=str(round(gamma,2)))
	yscale('log')
	box = ax.get_position()
	legend(bbox_to_anchor=(-0.1,-0.25,1.,0.1),ncol=6, borderaxespad=0.,prop={'size':10},fancybox=True)
	ax.set_position([box.x0, box.y0 + box.height * 0.1,\
	                 box.width, box.height * 0.9])
	xlabel('Highpass Scale [h]')
	ylabel('Min Capacity')
	title('Storage Capacities With Highpass db6 Filtering')
	ylim(10**-1,10**3)
	axvline(365*24,ls='--',c='r')
	annotate('Yearly',xy=(365*24*1.1,1.1*10**-1))
	axvline(120*24,ls='--',c='r')
	annotate('Seasonal',xy=(120*24*1.1,1.5*10**-1))
	axvline(30*24,ls='--',c='r')
	annotate('Monthly',xy=(30*24*1.1,1.1*10**-1))
	axvline(14*24,ls='--',c='r')
	annotate('Synoptic',xy=(14*24*1.1,1.5*10**-1))
	axvline(7*24,ls='--',c='r')
	annotate('Weekly',xy=(7*24*1.1,1.1*10**-1))
	axvline(24,ls='--',c='r')
	annotate('Daily',xy=(24*1.1,1.5*10**-1))
	xscale('log')
	interactive(1)
	show()

def plot_mismatch_moments():
   def cbrt(x):
      return np.power(np.abs(x),1./3)*np.sign(x)
   
   gammas=np.linspace(0.0,2.0,400)
   mean=[]
   var=[]
   skewness=[]
   kurtosis=[]
   for gamma in gammas:
      ts=get_mismatch(gamma)
      n=len(ts)
      curmean=np.sum(ts)/n
      mean.append(curmean)
      var.append(np.sum((ts-curmean)**2)/n)
      skewness.append(np.sum((ts-curmean)**3)/n)
      kurtosis.append(np.sum((ts-curmean)**4)/n)
   
   mean=np.array(mean)
   var=np.array(var)
   skewness=np.array(skewness)
   kurtosis=np.array(kurtosis)
   
   fig=plt.figure()
   ax=plt.axes()
   plt.plot(gammas,mean,label='Mean')
   plt.plot(gammas,var,label='Variance')
   plt.plot(gammas,skewness,label='Skewness')
   plt.plot(gammas,kurtosis,label='Kurtosis')
   plt.xlabel("Gamma")
   plt.ylabel("Moment")
   plt.title("N'th central moments")
   plt.legend(loc='lower right')
   plt.axhline(0,linestyle='--',c='k')
   
   fig=plt.figure()
   ax=plt.axes()
   plt.plot(gammas,mean,label='Mean')
   plt.plot(gammas,var**(1./2),label='Variance')
   plt.plot(gammas,cbrt(skewness),label='Skewness')
   plt.plot(gammas,kurtosis**(1./4),label='Kurtosis')
   plt.xlabel("Gamma")
   plt.ylabel("Moment^(1/n)")
   plt.title("N'th roots of N'th central moments")
   plt.legend(loc='lower right')
   plt.axhline(0,linestyle='--',c='k')
   
   plt.interactive(1)
   plt.show()


def plot_slts_selection():
   fig=plt.figure()
   gammas=[0.8,0.95,1.0,1.05]
   ts=get_mismatch(gammas[0])
   dummy,dummy,slts1=get_policy_2_storage(ts,return_storage_filling_time_series=True)
   ts=get_mismatch(gammas[1])
   dummy,dummy,slts2=get_policy_2_storage(ts,return_storage_filling_time_series=True)
   ts=get_mismatch(gammas[2])
   dummy,dummy,slts3=get_policy_2_storage(ts,return_storage_filling_time_series=True)
   ts=get_mismatch(gammas[3])
   dummy,dummy,slts4=get_policy_2_storage(ts,return_storage_filling_time_series=True)
   
   plt.subplot(221)
   plt.plot(np.arange(len(slts1)*1.0)/(364*24),slts1)
   #plt.xlabel("Year")
   plt.ylabel("Storage Level")
   plt.title("Storage Level, gamma="+str(gammas[0]))
   plt.xlim((0,8))
   
   plt.subplot(222)
   plt.plot(np.arange(len(slts2)*1.0)/(364*24),slts2)
   #plt.xlabel("Year")
   #plt.ylabel("Storage Level")
   plt.title("Storage Level, gamma="+str(gammas[1]))
   plt.xlim((0,8))
   
   plt.subplot(223)
   plt.plot(np.arange(len(slts3)*1.0)/(364*24),slts3)
   plt.xlabel("Year")
   plt.ylabel("Storage Level")
   plt.title("Storage Level, gamma="+str(gammas[2]))
   plt.xlim((0,8))
   
   plt.subplot(224)
   plt.plot(np.arange(len(slts4)*1.0)/(364*24),slts4)
   plt.xlabel("Year")
   #plt.ylabel("Storage Level")
   plt.title("Storage Level, gamma="+str(gammas[3]))
   plt.xlim((0,8))
   
   plt.interactive(1)
   plt.show()


