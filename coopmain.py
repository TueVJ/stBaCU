from coopgrid import *
from cooppolicies import *
#from coopplot import *
import countries
import matplotlib.pyplot as plt


#ts=get_mismatch();
#reduced_mismatch,used_storage=get_policy_2_storage(ts,storage_capacity=6);
#some_zeros,total_storage=get_policy_2_storage(ts);
#plt.plot(ts)
#plt.plot(reduced_mismatch)
#print used_storage
#print total_storage

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

#save raw data as gnuplot data file
#file=open("gnuplot.data","w")
#string="gamma raw\n"
#for i in arange(N):
#    string+=str(gammas[i])+" "+str(storage_capacities_raw[i])+"\n"
#print string
#file.write(string)
#file.close
