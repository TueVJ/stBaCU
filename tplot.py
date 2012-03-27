from mfunc import *
from mpolicies import *
#from coopplot import *
from tfunc import *
#import countries
import matplotlib.pyplot as plt
import copy
import pywt
import numpy as np

#ts=get_mismatch();
#reduced_mismatch,used_storage=get_policy_2_storage(ts,storage_capacity=6);
#some_zeros,total_storage=get_policy_2_storage(ts);
#plt.plot(ts)
#plt.plot(reduced_mismatch)
#print used_storage
#print total_storage
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
	num_smoothers=2#len(smoothers)
	
	storage_capacities_raw=[[0 for i in range(N)] for j in range(num_smoothers)]
	storage_capacities_small=[[0 for i in range(N)] for j in range(num_smoothers)]
	for j in arange(num_smoothers):
		for i in arange(N):
			ts=get_mismatch(gammas[i])
			smoothed_ts=get_smoothed(ts ,smoothers[j])
			dummy,storage_capacities_raw[j][i]=get_policy_2_storage(smoothed_ts)
			after_small,dummy=get_policy_2_storage(ts,storage_capacity=small_storage_size)
			dummy,storage_capacities_small[j][i]=get_policy_2_storage(get_smoothed(after_small,smoothers[j]))
	
	for j in arange(num_smoothers):
		plt.plot(gammas,storage_capacities_raw[j],label="No storage, "+smoothernames[j])
		plt.plot(gammas,storage_capacities_small[j],label=str(small_storage_size)+"h storage, "+smoothernames[j])
	plt.gca().set_ylim(bottom=1)
	yscale("log")
	legend(loc="upper right")
	interactive(1)
	show()

def plot_wavelets_smoothed_capacities():
	gammas=np.linspace(0.5,1,11)
	results={}
	for gamma in gammas:
		ts = get_mismatch(gamma)
		coeffs=pywt.wavedec(ts,'db2')
		steps=len(coeffs)-np.arange(len(coeffs)-3)-3
		maxes={}
		for step in steps:
			print str(step) + ', '+str(gamma)
			coeffs[step]*=0
			ts=pywt.waverec(coeffs,'db2')
			dummy,storage_capacity=get_policy_2_storage(ts)
			maxes[step]=storage_capacity
		results[gamma]=copy.deepcopy(maxes)
	
	for gamma in gammas:
		plot(results[gamma].keys(),results[gamma].values(),label=str(round(gamma,2)))
	yscale('log')
	#legend(loc='lower right')
	interactive(1)
	show()
