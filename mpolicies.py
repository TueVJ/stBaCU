from pylab import NaN
from numpy import array
import numpy
from scipy.optimize import *
from scipy.integrate import *
from scipy.stats.mstats import mquantiles
import constant
from copy import deepcopy
from policies import *
from policy_helpers import *
from mutils import pos,neg

## Wrapper
## Policy: Minimimze balancing power on connected components.
def get_storage_and_new_mismatch(mismatch, offset = 0, fast_eta_in = 1., fast_eta_out = 1., slow_eta_in = .6, slow_eta_out = .6, full_output = True):
	return get_storage_min_balancing_power(mismatch, offset = offset, fast_eta_in = fast_eta_in, fast_eta_out = fast_eta_out, slow_eta_in = slow_beta_in, slow_eta_out = slow_eta_out, full_output = full_output)

## Wrapper
def get_policy_1_storage(ts, eta_in, eta_out, storage_capacity = NaN):
	return get_storage_balancing_first(ts,constant.policy_1_storage_fraction*eta_in/constant.policy_1_fraction,eta_out,storage_capacity=storage_capacity,storage_fraction=constant.policy_1_fraction)

def get_policy_5_storage(ts, eta_in, eta_out, storage_capacity):
#    plot(ts, 'k', lw = 2)
    storage_part = ts * (ts <= constant.policy_5_in_capacity) * (ts >= -constant.policy_5_out_capacity) + constant.policy_5_in_capacity * (ts > constant.policy_5_in_capacity) - constant.policy_5_out_capacity * (ts < -constant.policy_5_out_capacity)
    remainder = ts - storage_part
#    plot(remainder + 2)
    storage_ts, used_storage = get_policy_2_storage(storage_part, eta_in, eta_out, storage_capacity)
#    plot(storage_ts + remainder, 'r')

#    plot(storage_ts + remainder - ts - 2)
    return storage_ts + remainder, used_storage

## Wrapper
###
# Policy 2: Dumb storage:
# If mismatch is positive, throw that much into storage.
# If mismatch is negative, throw as much as we can 
# out of storage.
#	@param: storage_capacity: Maximum storage capacity
#			in units of average consumption per hour.
###
def get_policy_2_storage(ts, eta_in = 1., eta_out = 1., storage_capacity = NaN,return_storage_filling_time_series=False):
	return get_storage_max_power(ts=ts,eta_in=eta_in,eta_out=eta_out,storage_capacity=storage_capacity,return_storage_filling_time_series=return_storage_filling_time_series)

## Wrapper
## Policy: minimize storage extraction/insertion power
## on connected components.
def get_policy_4_storage(ts, eta_in = 1., eta_out = 1., storage_capacity = NaN):
    """Policy 4"""
    return get_storage_min_storage_power(ts=ts,eta_in=eta_in,eta_out=eta_out,storage_capacity = storage_capacity)


def get_policy_7_storage(ts, eta_in, eta_out, storage_capacity):
    storage_part = pos(ts) - pos(neg(ts) - constant.policy_7_capacity)
    if storage_part.min() >= 0:
        return ts, 0
    remainder = ts - storage_part
    storage_ts, used_storage = get_policy_2_storage(storage_part, eta_in, eta_out, storage_capacity)
    return storage_ts + remainder, used_storage

def get_policy_8_storage(ts, eta_in, eta_out, storage_capacity):
    storage_part = pos(ts) - (neg(ts) - constant.policy_8_capacity)
    if storage_part.min() >= 0:
        return ts, 0
    remainder = ts - storage_part
    storage_ts, used_storage = get_policy_2_storage(storage_part, eta_in, eta_out, storage_capacity)
    return storage_ts + remainder, used_storage

def get_policy_9_storage(ts, eta_in, eta_out, storage_capacity):
    if pos(ts).sum() * eta_in * eta_out < neg(ts).sum():
        beta = (pos(ts).sum() * eta_in * eta_out) / neg(ts).sum()
    else:
        beta = 1.
    storage_part = pos(ts) - beta * neg(ts)
    remainder = ts - storage_part
    storage_ts, used_storage = get_policy_2_storage(storage_part, eta_in, eta_out, storage_capacity)
    print used_storage
    return storage_ts + remainder, used_storage



def old_storage(ts, eta_in = 1., eta_out = 1., storage_capacity = NaN):
    indices, integrals = get_indices_and_integrals(ts)
    storing_potential = eta_in * eta_out * sum(integrals * (integrals > 0))
    extracting_potential = -sum(integrals * (integrals < 0))

    storage_start, used_storage = get_storage_start(integrals, eta_in, eta_out, 
                                      storage_capacity, 
                                      storing_potential > extracting_potential)
    
    storage_usage_ts = get_storage_usage_ts(integrals, storage_start, eta_in, eta_out, storage_capacity, storing_potential > extracting_potential)
    
    extracting = zeros_like(ts)
    storing = zeros_like(ts)

    end_to_start_slice = concatenate((ts[indices[-1] + 1:],ts[:indices[0] + 1]))
    if storage_usage_ts[0] < 0:
        end_to_start_extracting = pos(-end_to_start_slice - find_height(-end_to_start_slice, -storage_usage_ts[0]))
        extracting[indices[-1] + 1:] = -end_to_start_extracting[:len(extracting[indices[-1] + 1:])]
        extracting[:indices[0] + 1] = -end_to_start_extracting[len(extracting[indices[-1] + 1:]):]
    else:
        end_to_start_storing = pos(end_to_start_slice - find_height(end_to_start_slice, storage_usage_ts[0]))
        storing[indices[-1] + 1:] = end_to_start_storing[:len(extracting[indices[-1] + 1:])]
        storing[:indices[0] + 1] = end_to_start_storing[len(extracting[indices[-1] + 1:]):]
    
    for x in arange(len(indices)-1) + 1:
        #print x
        if storage_usage_ts[x] < 0:
            #print storage_usage_ts[x]
            #print indices[x]
            extracting[indices[x - 1] + 1: indices[x] + 1] = -pos(-ts[indices[x - 1] + 1: indices[x] + 1] - find_height(-ts[indices[x - 1] + 1: indices[x] + 1], -storage_usage_ts[x]))
        else:
            storing[indices[x - 1] + 1: indices[x] + 1] = pos(ts[indices[x - 1] + 1: indices[x] + 1] - find_height(ts[indices[x - 1] + 1: indices[x] + 1], storage_usage_ts[x]))

    #plot(storing, lw = 2)
    #plot(extracting, lw = 2)
    #plot(ts)
    #plot(ts - storing - extracting+ 1)

    return(ts - storing - extracting)

def get_mismatch_after_double_storage(mismatch, fast_eta_in = .9, fast_eta_out = .9, slow_eta_in = .6, slow_eta_out = .6, fast_storage_capacity = NaN, slow_storage_capacity = NaN, cut_time = 24., returnmismatchonly = False, slow_storage_speed = NaN):
    smoothed = get_smoothed(mismatch, cut_time)
    #print mismatch
    #print "cut_time: ", cut_time
    #print smoothed
    hf = mismatch - smoothed
    if isnan(slow_storage_speed):
        slow_storage_speed = cut_time

    offset = find_offset(hf, fast_eta_in, fast_eta_out)
    
    if fast_storage_capacity <= 0.:
        fast_mismatch = hf - offset
        fast_storage_size = 0.
    else:
        fast_mismatch, fast_storage_size = get_storage(hf - offset, fast_eta_in, fast_eta_out, fast_storage_capacity)

    #print "slow_storage_speed: ", slow_storage_speed
    lf = smoothed + get_smoothed(fast_mismatch * (hf - offset - fast_mismatch), slow_storage_speed)
    #print lf - smoothed
    hf_mismatch = fast_mismatch + offset + smoothed - lf
    if (lf.max() + offset <= 0) or (lf.min() +offset >= 0):
        slow_mismatch = lf + offset
        slow_storage_size = 0.
    else:
        #print "lf, offset, slow_eta's, slow_storage_capacity: ", lf, offset, slow_eta_in, slow_eta_out, slow_storage_capacity
        slow_mismatch, slow_storage_size = get_storage(lf + offset, slow_eta_in, slow_eta_out, slow_storage_capacity)
    lf_mismatch = slow_mismatch - offset
    if returnmismatchonly:
        return hf_mismatch + lf_mismatch
    return hf_mismatch + lf_mismatch, fast_storage_size, slow_storage_size

def get_mismatch_after_policy_3_storage(mismatch, eta_in = .6, eta_out = .6, storage_capacity = NaN):
    new_mismatch, storage_size = get_storage(mismatch, eta_in, eta_out, storage_capacity)
    return new_mismatch
    
def get_mismatch_after_policy_4_storage(mismatch, eta_in = 1., eta_out = 1., storage_capacity = NaN):
    new_mismatch, storage_size = get_policy_4_storage(mismatch, eta_in, eta_out, storage_capacity)
    return new_mismatch

def get_mismatch_after_policy_2_storage(mismatch, eta_in = 1., eta_out = 1., storage_capacity = NaN):
    new_mismatch, storage_size = get_policy_2_storage(mismatch, eta_in, eta_out, storage_capacity)
    return new_mismatch

def get_mismatch_after_policy_5_storage(mismatch, eta_in = 1., eta_out = 1., storage_capacity = NaN):
    new_mismatch, storage_size = get_policy_5_storage(mismatch, eta_in, eta_out, storage_capacity)
    return new_mismatch

def get_mismatch_after_policy_1_storage(mismatch, eta_in = 1., eta_out = 1., storage_capacity = NaN):
    new_mismatch, storage_size = get_policy_1_storage(mismatch, eta_in, eta_out, storage_capacity)
    return new_mismatch

def get_mismatch_after_policy_7_storage(mismatch, eta_in = 1., eta_out = 1., storage_capacity = NaN):
    new_mismatch, storage_size = get_policy_7_storage(mismatch, eta_in, eta_out, storage_capacity)
    return new_mismatch

def get_mismatch_after_policy_8_storage(mismatch, eta_in = 1., eta_out = 1., storage_capacity = NaN):
    new_mismatch, storage_size = get_policy_8_storage(mismatch, eta_in, eta_out, storage_capacity)
    return new_mismatch

def get_mismatch_after_policy_9_storage(mismatch, eta_in = 1., eta_out = 1., storage_capacity = NaN):
    new_mismatch, storage_size = get_policy_9_storage(mismatch, eta_in, eta_out, storage_capacity)
    return new_mismatch

