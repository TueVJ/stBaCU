from pylab import NaN
from numpy import array
import numpy
from scipy.optimize import *
from scipy.integrate import *
from scipy.stats.mstats import mquantiles
from policy_helpers import *
from utils import pos,neg


## Policy: Minimimze balancing power on connected components.
def get_storage_min_balancing_power(mismatch, offset = 0, fast_eta_in = 1., fast_eta_out = 1., slow_eta_in = .6, slow_eta_out = .6, full_output = True):

    slow_excess = pos(get_smoothed(mismatch - offset, 24))
    slow_deficiency = neg(get_smoothed(mismatch - offset, 24))
    #print "sum(slow_storing): ", sum(slow_storing)
    h = find_height(-mismatch, (slow_eta_in * slow_eta_out) * sum(slow_storing))
    if h < 0:
        print "WARNING!"
    slow_extracting = pos(-mismatch - h)
    #print "sum(slow_extracting): ", sum(slow_extracting)
 
    excess = pos(mismatch)
 
    fast_excess = pos(excess - slow_storing)

    deficiency = neg(mismatch) + neg(excess - slow_storing)

    residual_deficiency = pos(deficiency - slow_extracting)

    fast_storing = fast_excess + neg(deficiency - slow_extracting)
    #print "sum(fast_storing): ", sum(fast_storing)

    k = find_height(residual_deficiency, sum(fast_storing)*(fast_eta_in * fast_eta_out))
    if k < 0:
        print "WARNING!"
    fast_extracting = pos(residual_deficiency - k)
    #print "sum(fast_extracting): ", sum(fast_extracting)


    new_mismatch = mismatch + (slow_extracting + fast_extracting) - (slow_storing + fast_storing)

    if full_output:
        return new_mismatch, fast_storing, fast_extracting, slow_storing, slow_extracting
    else: return new_mismatch


## Assign storage_fraction of the mismatch to the storage
## the rest goes to balancing.
def get_storage_balancing_first(ts, eta_in, eta_out, storage_capacity = NaN, storage_fraction=0.5):
    fractional_balancing = (1. - storage_fraction) * ts
    fractional_storage = storage_fraction * ts
    storage_ts, used_storage = get_storage_max_power(fractional_storage, eta_in, eta_out, storage_capacity)
    return storage_ts + fractional_balancing, used_storage


###
# Policy 2: Dumb storage:
# If mismatch is positive, throw that much into storage.
# If mismatch is negative, throw as much as we can 
# out of storage.
#	@param: storage_capacity: Maximum storage capacity
#			in units of average consumption per hour.
###
def get_storage_max_power(ts, eta_in = 1., eta_out = 1., storage_capacity = NaN,return_storage_filling_time_series=False):
    """Policy 2"""
    if ts.min() >= -1e-10 or ts.max() <= 1e-10:
        return ts, 0
    indices, integrals = get_indices_and_integrals(ts)
    storing_potential = eta_in * eta_out * sum(integrals * (integrals > 0))
    extracting_potential = -sum(integrals * (integrals < 0))

    storage_start, used_storage = get_storage_start(integrals, eta_in, eta_out, 
                                      storage_capacity, 
                                      storing_potential > extracting_potential)
    #print storage_start, used_storage
    
    storage_usage_ts = get_storage_usage_ts(integrals, storage_start, eta_in, eta_out, storage_capacity, storing_potential > extracting_potential)
    
    #print (pos(storage_usage_ts) * eta_in - neg(storage_usage_ts) / eta_in).sum()

    extracting = zeros_like(ts)
    storing = zeros_like(ts)

	## If start and end have the same sign, we need to ensure distribution
	## across the cyclic boundaries.
    end_to_start_slice = concatenate((ts[indices[-1] + 1:],ts[:indices[0] + 1]))
    if storage_usage_ts[0] < 0:
        end_to_start_extracting = -burn_off(-end_to_start_slice, -storage_usage_ts[0])
        extracting[indices[-1] + 1:] = end_to_start_extracting[:len(extracting[indices[-1] + 1:])]
        extracting[:indices[0] + 1] = end_to_start_extracting[len(extracting[indices[-1] + 1:]):]
    else:
        end_to_start_storing = burn_off(pos(end_to_start_slice), pos(storage_usage_ts[0]))
        storing[indices[-1] + 1:] = end_to_start_storing[:len(extracting[indices[-1] + 1:])]
        storing[:indices[0] + 1] = end_to_start_storing[len(extracting[indices[-1] + 1:]):]
    ## Distribute across all other intervals.
    for x in arange(len(indices)-1) + 1:
        current_slice = ts[indices[x - 1] + 1: indices[x] + 1]
        #print x
        if storage_usage_ts[x] < 0:
            #print indices[x]
            extracting[indices[x - 1] + 1: indices[x] + 1] = -burn_off(-current_slice, -storage_usage_ts[x])
        else:
            storing[indices[x - 1] + 1: indices[x] + 1] = burn_off(pos(current_slice), pos(storage_usage_ts[x]))

#    figure()
    #plot(storing, lw = 5)
#    #plot(extracting, lw = 10)
    #plot(ts)
    #plot(ts - storing - extracting+ 1)
#    figure()
    storage_filling_with_offset = eta_in * storing.cumsum() + extracting.cumsum() / eta_out
    storage_filling_without_offset = storage_filling_with_offset - storage_filling_with_offset.min()
    used_storage = storage_filling_with_offset.max() - storage_filling_with_offset.min()
    #plot(storage_filling_with_offset)
    if return_storage_filling_time_series:
	    return ts - storing - extracting, used_storage, storage_filling_without_offset
    else:
        return ts - storing - extracting, used_storage
		

## Policy: minimize storage extraction/insertion power
## on connected components.
def get_storage_min_storage_power(ts, eta_in = 1., eta_out = 1., storage_capacity = NaN):
    indices, integrals = get_indices_and_integrals(ts)
    storing_potential = eta_in * eta_out * sum(integrals * (integrals > 0))
    extracting_potential = -sum(integrals * (integrals < 0))

    storage_start, used_storage = get_storage_start(integrals, eta_in, eta_out, 
                                      storage_capacity, 
                                      storing_potential > extracting_potential)
    #print storage_start, used_storage
    
    storage_usage_ts = get_storage_usage_ts(integrals, storage_start, eta_in, eta_out, storage_capacity, storing_potential > extracting_potential)
    
    extracting = zeros_like(ts)
    storing = zeros_like(ts)

    end_to_start_slice = concatenate((ts[indices[-1] + 1:],ts[:indices[0] + 1]))
    if storage_usage_ts[0] < 0:
        end_to_start_extracting = -end_to_start_slice - pos(-end_to_start_slice - find_height(-end_to_start_slice, -end_to_start_slice.sum() + storage_usage_ts[0]))
        extracting[indices[-1] + 1:] = -end_to_start_extracting[:len(extracting[indices[-1] + 1:])]
        extracting[:indices[0] + 1] = -end_to_start_extracting[len(extracting[indices[-1] + 1:]):]
    else:
        end_to_start_storing = pos(end_to_start_slice - pos(end_to_start_slice - find_height(end_to_start_slice, end_to_start_slice.sum() - storage_usage_ts[0])))
        storing[indices[-1] + 1:] = end_to_start_storing[:len(extracting[indices[-1] + 1:])]
        storing[:indices[0] + 1] = end_to_start_storing[len(extracting[indices[-1] + 1:]):]
    
    for x in arange(len(indices)-1) + 1:
        current_slice = ts[indices[x - 1] + 1: indices[x] + 1]
        #print x
        if storage_usage_ts[x] < 0:
            #print indices[x]
            extracting[indices[x - 1] + 1: indices[x] + 1] = current_slice + pos(-current_slice - find_height(-current_slice, -current_slice.sum() + storage_usage_ts[x]))
        else:
            storing[indices[x - 1] + 1: indices[x] + 1] = pos(current_slice - pos(current_slice - find_height(current_slice, current_slice.sum() - storage_usage_ts[x])))

#    figure()
    #plot(storing, lw = 5)
#    #plot(extracting, lw = 10)
    #plot(ts)
    #plot(ts - storing - extracting+ 1)
#    figure()
    storage_filling_with_offset = eta_in * storing.cumsum() + extracting.cumsum() / eta_out
    used_storage = storage_filling_with_offset.max() - storage_filling_with_offset.min()
    return ts - storing - extracting, used_storage

