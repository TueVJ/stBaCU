from pylab import *
import numpy
from scipy.optimize import *
from scipy.integrate import *
from utils import pos


###
#   DESC:
#   Given an array of data, returns an array of the same
#   length, containing the original data up until the
#   index where cumsum(array)>number. Array(index) is
#   replaced by the remainder af subtraction from number,
#   and trailing zeros are inserted.
#   USAGE:
#   Given a storage level of number, and an array of 
#   power consumption data, return an array of power
#   usage when running completely on the storage.
###

def burn_off(array, number):
    """array and number must be positive and array.sum() >= number"""
    if array.sum() > number:
        index = argmax(array.cumsum() > number)
    else:
        return array
    return_array = zeros_like(array)
    return_array[:index] = array[:index]
    remainder = number - return_array.sum()
    return_array[index] = remainder
    return return_array


###
#	@param: negative_numbers: BOOL indicating whether ts bottoms out.
###
def get_storage_usage_ts(integrals, storage_start, eta_in, eta_out, storage_capacity, negative_numbers):
#    print negative_numbers
#    print storage_start
    in_out_data = eta_in * integrals * (integrals > 0) + integrals * (integrals < 0) / eta_out
    if negative_numbers:
        full_storage = 0.
        empty_storage = -storage_capacity
    else:
        full_storage = storage_capacity
        empty_storage = 0.
    suts = empty_like(integrals) #Storage usage time series
    old_storage_level = storage_start
    if size(in_out_data) == 1:
        in_out_data = array([in_out_data])
#    print in_out_data[0]
    if in_out_data[0] >= 0:
        for x in 2 * arange(len(in_out_data)/2):
            new_storage_level = old_storage_level + in_out_data[x]
            if new_storage_level > full_storage:
                new_storage_level = full_storage
            suts[x] = (new_storage_level - old_storage_level) / eta_in
            old_storage_level = new_storage_level
            new_storage_level = old_storage_level + in_out_data[x + 1]
            if new_storage_level < empty_storage:
                new_storage_level = empty_storage
            suts[x + 1] = (new_storage_level - old_storage_level) * eta_out
            old_storage_level = new_storage_level
    else:
        for x in 2 * arange(len(in_out_data)/2):
            new_storage_level = old_storage_level + in_out_data[x]
            if new_storage_level < empty_storage:
                new_storage_level = empty_storage
            suts[x] = (new_storage_level - old_storage_level) * eta_out
            old_storage_level = new_storage_level
            new_storage_level = old_storage_level + in_out_data[x + 1]
            if new_storage_level > full_storage:
                new_storage_level = full_storage
            suts[x + 1] = (new_storage_level - old_storage_level) / eta_in
            old_storage_level = new_storage_level
#    plot((pos(suts) * eta_in - neg(suts) / eta_out).cumsum())
    return suts
    

###
#	@param:	storage_capacity: indicates maximum storage capacity (possibly NaN)
#			start_at_top: BOOL indicating whether to start with a full storage or not.
#	@return: storage_level: The amount which the storage should start from
#			in order to be neutral, or cyclic.
#			used_storage:	LEGACY: A first guess on how much storage capacity is used.
###

def get_storage_start(integrals, eta_in, eta_out, storage_capacity, start_at_top):
    in_out_data = eta_in * integrals * (integrals > 0) + integrals * (integrals < 0) / eta_out

    used_storage = 0.

#    print start_at_top
    
    if start_at_top:
        full_storage = 0.
        empty_storage = -storage_capacity
    else:
        full_storage = storage_capacity
        empty_storage = 0.
    if size(in_out_data) == 1:
        in_out_data = array([in_out_data,0])
    if (in_out_data[0] < 0 ^ start_at_top):
        storage_level = 0.
    else:
        storage_level = min(in_out_data[0], full_storage)
#        used_storage = max(abs(used_storage), abs(storage_level))
    if in_out_data[0] > 0:
        storage_level = storage_level + in_out_data[1]
        if storage_level < empty_storage:
            storage_level = empty_storage
#        used_storage = max(abs(used_storage), abs(storage_level))
        start_index = 2
    else:
#        used_storage = max(abs(used_storage), abs(storage_level))
        start_index = 1
#    print start_index, storage_level, full_storage, empty_storage
    for x in 2 * arange(len(in_out_data)/2 - 1) + start_index:
        storage_level = storage_level + in_out_data[x]
        if storage_level > full_storage:
            storage_level = full_storage
        used_storage = max(abs(used_storage), abs(storage_level))
        storage_level = storage_level + in_out_data[x + 1]
        if storage_level < empty_storage:
            storage_level = empty_storage
        used_storage = max(abs(used_storage), abs(storage_level))

    if start_index == 1:
        storage_level = storage_level + in_out_data[len(in_out_data) - 1]
        if storage_level > full_storage:
            storage_level = full_storage
        used_storage = max(abs(used_storage), abs(storage_level))

#    print storage_level, used_storage
    return storage_level, used_storage

def get_indices_and_integrals(ts):
    if min(ts) >= 0:
        print "INTET NEGATIVT!"
    if sum(ts == 0) != 0:
        #print "zeros found!"
        new_ts = deepcopy(ts)
        for i in find(ts == 0):
            if i + 1 < len(ts):
                new_ts[i] = (ts[i-1]+ts[i+1])/2 + 1e-10
            else:
                new_ts[i] = (ts[i-1]+ts[0])/2 + 1e-10
        if sum(new_ts == 0) != 0:
            #print "Crap!!"
            return
        else:
            indices = find(diff(sign(new_ts))!=0)
    else:
        indices = find(diff(sign(ts))!=0)
    length = len(indices)
    if length < 2:
        print "Too few crossings:", length, min(ts), indices
        return 0, sum(ts)
    integral = ts.cumsum()
    component_integral = integral[indices]
    if length/2 == length/2.:
        #print "odd"
        component_integrals = empty(length)
        component_integrals[0] = component_integral[0] + (ts.sum() - 
                                                          component_integral[length - 1])
        component_integrals[1:] = diff(component_integral)
        return indices, component_integrals
    else:
        #print "even"
        component_integrals = empty(length + 1)
        component_integrals[0] = component_integral[0]
        component_integrals[1:-1] = diff(component_integral)
        component_integrals[length] = ts.sum() - component_integral[length - 1]
        return concatenate([indices,[len(ts) - 1]]), component_integrals


def find_height(array, constant):
    if constant < 0:
        return array.max()
    def funct(k):
        return sum(pos(array - k)) - constant
#    #print funct(-sum(array)*constant), funct(-1e10), -sum(array)*constant
    integral = sum(array)
    if not constant == 0:
        if not integral == 0:
            return brenth(funct, -(integral * constant + 1000), integral * constant + 1000)
        else:
            return brenth(funct, -constant**2, constant**2)
    else:
        return brenth(funct, -max(array)*1000, max(array)*1000)

def find_offset(ts, eta_in, eta_out):
    #print eta_in, eta_out, ts
    def funct(k):
        return eta_in * eta_out * pos(ts - k).sum() - neg(ts - k).sum()
    return brenth(funct, ts.min(), ts.max(), xtol = 0.)

