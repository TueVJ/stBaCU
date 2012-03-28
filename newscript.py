from init import *

ts=get_mismatch(1.01)
ss24h=wavelet_smoother(ts,short_cut=24*2)
sl3m=wavelet_smoother(ts,long_cut=24*31*3)
dummy,scus,slus=get_policy_2_storage(ts,return_storage_filling_time_series=True)
dummy,scss24h,slss24h=get_policy_2_storage(ss24h,return_storage_filling_time_series=True)
dummy,scsl3m,slsl3m=get_policy_2_storage(sl3m,return_storage_filling_time_series=True)
plot(slus,label='unsmoothed')
plot(slss24h,label='short24h')
plot(slsl3m,label='long3m')
legend()
interactive(1)
show()
