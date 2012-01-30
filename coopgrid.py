from pylab import *
import numpy
from scipy.optimize import *
from scipy.integrate import *
from scipy.stats.mstats import mquantiles
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import countries
import constant
from copy import deepcopy

def get_time_sum(D, time = 24):
    return D.reshape((-1, time)).sum(1)

def integrate_arrays(X,Y):
    if len(X) != len(Y):
        return NaN
    length = len(X)
    sum = 0.
    for d in arange(length - 1):
        sum = sum + (X[d + 1] - X[d]) * (Y[d] + Y[d + 1]) / 2
    return sum

def get_average_slopes_change_rates_and_projected_endpoints(X,Y):
    length = len(X)
    A = ones_like(X)
    B = ones_like(X)
    C = empty_like(X)
    for i in arange(length - 1):
        A[i] = 2 * (integrate_arrays(X[i:], Y[i:]) - (X[-1] - X[i]) * Y[-1]) / (Y[i] - Y[-1])
        B[i] = Y[i] - Y[-1]
    for i in arange(length - 1):
        C[i] = (A[i] / B[i] - A[i + 1] / B[i + 1]) / (A[i] / B[i])
    return A / B, C, X + A

def get_final_slope_endpoint_and_end_balancing(X, Y):
    t1, t2, t3 = get_average_slopes_change_rates_and_projected_endpoints(X, Y)
    i = len(X) - 2
    while not abs((t1*t2)[i]) < 100:
        i = i - 1
    return t1[i], t3[i], interp(t3[i], X, Y)

def get_characteristic_data_of_convex_curve(X, Y):
    balancing_only = Y[0]
    final_slope, storage_only, base_balancing = get_final_slope_endpoint_and_end_balancing(X, Y)
    integral = integrate_arrays(X, Y) - base_balancing * storage_only

    balancing_at_half_storage = interp(storage_only / 2, X, Y)
    storage_at_half_balancing = interp(balancing_only / 2, flipud(Y), flipud(X))

    initial_slope = integrate_arrays(concatenate([X[X < storage_at_half_balancing], [storage_at_half_balancing]]), concatenate([Y[X < storage_at_half_balancing],[balancing_only / 2]])) - storage_at_half_balancing * balancing_only / 2

    final_slope = integrate_arrays(concatenate([[storage_only / 2], X[X > storage_only / 2]]), concatenate([[balancing_at_half_storage], Y[X > storage_only / 2]])) - storage_only / 2 * base_balancing

    normalized_final_slope = 8 * final_slope / (storage_only * balancing_only)
    normalized_initial_slope = 8 * initial_slope / (storage_only * balancing_only)
    normalized_integral = 2 * integral / (storage_only * balancing_only)
    normalized_balancing_only = balancing_only / 70128
    normalized_storage_only = storage_only / 8766
    normalized_base_balancing = base_balancing / 70128
    return normalized_balancing_only, normalized_storage_only, normalized_initial_slope, normalized_final_slope, normalized_integral, normalized_base_balancing

def get_Wind(D):
    return D[:,0] * D[:,4] / D[:,0].mean()

def get_PV(D):
    if D[:,1].sum()==0:
        return D[:,1]
    return D[:,1] * D[:,4] / D[:,1].mean()

def get_Load(D):
    return D[:,2] * D[:,4] / D[:,2].mean()

def get_mismatch(gamma = 1., a = .6):
    D = countries.Region
    Wind = get_Wind(D = D)
    PV = get_PV(D = D)
    Load = get_Load(D = D)
    return gamma * (a * Wind + (1. - a) * PV) - Load

def get_smoothed(ts,cut_time=24):

    if cut_time == 0:
        return ts
   
    time = arange(len(ts));
    freq = fftfreq(len(time))[:len(time)/2+1]
    ts_fft = rfft(ts);

    sigma = cut_time*2;
    kernel = exp(-2*(pi*freq/(2*pi)*sigma)**2);
   
    return irfft(ts_fft*kernel);

def pos(array):
    return array * (array > 0) 

def neg(array):
    return -array * (array < 0)

def get_storage_and_new_mismatch(mismatch, offset = 0, fast_eta_in = 1., fast_eta_out = 1., slow_eta_in = .6, slow_eta_out = .6, full_output = True):

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
    suts = empty_like(integrals)
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
#def get_storage_usage_ts(integrals, storage_start, eta_in, eta_out, storage_capacity, negative_numbers):
#    in_out_data = eta_in * integrals * (integrals > 0) + integrals * (integrals < 0) / eta_out
#    if negative_numbers:
#        full_storage = 0.
#        empty_storage = -storage_capacity
#    else:
#        full_storage = storage_capacity
#        empty_storage = 0.
#    suts = empty_like(integrals)
#    old_storage_level = storage_start
#    if in_out_data[0] > 0:
#        for x in 2 * arange(len(in_out_data)/2):
#            new_storage_level = old_storage_level + in_out_data[x]
#            if new_storage_level > full_storage:
#                new_storage_level = full_storage
#            suts[x] = new_storage_level / eta_in
#            old_storage_level = new_storage_level
#            new_storage_level = old_storage_level + in_out_data[x + 1]
#            if new_storage_level < empty_storage:
#                new_storage_level = empty_storage
#            suts[x + 1] = new_storage_level * eta_out
#            old_storage_level = new_storage_level
#    else:
#        for x in 2 * arange(len(in_out_data)/2):
#            new_storage_level = old_storage_level + in_out_data[x]
#            if new_storage_level < empty_storage:
#                new_storage_level = empty_storage
#            suts[x] = new_storage_level * eta_out
#            old_storage_level = new_storage_level
#            new_storage_level = old_storage_level + in_out_data[x + 1]
#            if new_storage_level > full_storage:
#                new_storage_level = full_storage
#            suts[x + 1] = new_storage_level / eta_in
#            old_storage_level = new_storage_level
#    #plot(suts)
#    return suts

def get_storage(ts, eta_in = 1., eta_out = 1., storage_capacity = NaN):
    """Policy 3"""
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

#    figure()
    #plot(storing, lw = 5)
#    #plot(extracting, lw = 10)
    #plot(ts)
    #plot(ts - storing - extracting+ 1)
#    figure()
    storage_filling_with_offset = eta_in * storing.cumsum() + extracting.cumsum() / eta_out
    used_storage = storage_filling_with_offset.max() - storage_filling_with_offset.min()
    return ts - storing - extracting, used_storage

def get_policy_4_storage(ts, eta_in = 1., eta_out = 1., storage_capacity = NaN):
    """Policy 4"""
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

def get_policy_2_storage(ts, eta_in = 1., eta_out = 1., storage_capacity = NaN):
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

    end_to_start_slice = concatenate((ts[indices[-1] + 1:],ts[:indices[0] + 1]))
    if storage_usage_ts[0] < 0:
        end_to_start_extracting = -burn_off(-end_to_start_slice, -storage_usage_ts[0])
        extracting[indices[-1] + 1:] = end_to_start_extracting[:len(extracting[indices[-1] + 1:])]
        extracting[:indices[0] + 1] = end_to_start_extracting[len(extracting[indices[-1] + 1:]):]
    else:
        end_to_start_storing = burn_off(pos(end_to_start_slice), pos(storage_usage_ts[0]))
        storing[indices[-1] + 1:] = end_to_start_storing[:len(extracting[indices[-1] + 1:])]
        storing[:indices[0] + 1] = end_to_start_storing[len(extracting[indices[-1] + 1:]):]
    
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
    used_storage = storage_filling_with_offset.max() - storage_filling_with_offset.min()
    #plot(storage_filling_with_offset)
    return ts - storing - extracting, used_storage

def get_policy_1_storage(ts, eta_in, eta_out, storage_capacity = NaN):
    balancing_fraction = (1. - constant.policy_1_fraction) * ts
    storage_fraction = constant.policy_1_fraction * ts
    storage_ts, used_storage = get_policy_2_storage(storage_fraction, constant.policy_1_storage_fraction * eta_in / constant.policy_1_fraction, eta_out, storage_capacity)
    return storage_ts + balancing_fraction, used_storage

def get_policy_5_storage(ts, eta_in, eta_out, storage_capacity):
#    plot(ts, 'k', lw = 2)
    storage_part = ts * (ts <= constant.policy_5_in_capacity) * (ts >= -constant.policy_5_out_capacity) + constant.policy_5_in_capacity * (ts > constant.policy_5_in_capacity) - constant.policy_5_out_capacity * (ts < -constant.policy_5_out_capacity)
    remainder = ts - storage_part
#    plot(remainder + 2)
    storage_ts, used_storage = get_policy_2_storage(storage_part, eta_in, eta_out, storage_capacity)
#    plot(storage_ts + remainder, 'r')

#    plot(storage_ts + remainder - ts - 2)
    return storage_ts + remainder, used_storage

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
    #print storage_start
    
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


def sliding_min(ts, windowsize = 12):
    length = len(ts)
    temp_ts = empty(length + windowsize)
    temp_ts[:length] = ts
    temp_ts[length:] = ts[:windowsize]
    sliding_min = empty(length)
    for i in range(length):
        sliding_min[i] = temp_ts[i:i + windowsize].min()
    return sliding_min

def smooth_min(ts, windowsize = 12):
    length = len(ts)
    bool_indices = empty(length, bool)
    bool_indices[1:-1] = (diff(ts[:-1]) < 0) * (diff(ts[1:]) > 0)
    bool_indices[0] = (ts[-1] - ts[0] > 0) * (ts[0] - ts[1] < 0)
    bool_indices[-1] = (ts[-2] - ts[-1] > 0) * (ts[-1] - ts[0] < 0)
    indices = arange(length)[bool_indices]
    for i in arange(len(indices[1:])):
        if indices[i] - indices[i - 1] < windowsize:
            if ts[indices[i]] > ts[indices[i - 1]]:
                bool_indices[indices[i]] = False
            else:
                bool_indices[indices[i - 1]] = False
        
    return interp(arange(length), arange(length)[bool_indices],ts[bool_indices])

def get_base_mismatch(ts, cut_time = 24, windowsize = 9):
    hf_offset = get_smoothed(smooth_min(ts, windowsize) - get_smoothed(ts, cut_time), cut_time)
    return get_smoothed(ts, cut_time) + hf_offset


def find_offset(ts, eta_in, eta_out):
    #print eta_in, eta_out, ts
    def funct(k):
        return eta_in * eta_out * pos(ts - k).sum() - neg(ts - k).sum()
    return brenth(funct, ts.min(), ts.max(), xtol = 0.)

def arch_function(x, unit, begin, end):
    arch = sqrt(1 - (x / unit - 1)**2)
    return begin * (1 - arch) + end * arch
 
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

def get_balancing_from_mismatch(ts):
    return neg(ts).sum()


def find_balancing_at_optimal_alpha(gamma, fast_eta_in = .9, fast_eta_out = .9, slow_eta_in = .6, slow_eta_out = .6, fast_storage_capacity = NaN, slow_storage_capacity = NaN, cut_time = 24, slow_storage_speed = NaN, returnalpha = False):
    def funct(a):
        #print "here: ", a, get_mismatch(gamma, a)
        return get_balancing_from_mismatch(get_mismatch_after_double_storage(get_mismatch(gamma, a), fast_eta_in, fast_eta_out, slow_eta_in, slow_eta_out, fast_storage_capacity, slow_storage_capacity, cut_time, slow_storage_speed = slow_storage_speed, returnmismatchonly = True))
    alpha, balancing, ierr, iterations = fminbound(funct,0.,1. , xtol = 1e-2, full_output = True)
    print "alpha: ", alpha, "iterations: ", iterations
    sys.stdout.flush()
    if returnalpha:
        return balancing, alpha 
    return balancing
    
def plot_slow_storage_vs_balance_graph(gamma, fast_eta_in = .9, fast_eta_out = .9, slow_eta_in = .6, slow_eta_out = .6, fast_storage_capacity = 6, slow_storage_capacity_array = array([0., 5., 10., 50., 100., 200., 500., 1000.]), cut_time = 24):
    used_storages = empty_like(slow_storage_capacity_array)
    optimal_alphas = empty_like(slow_storage_capacity_array)
    balancing_energies = empty_like(slow_storage_capacity_array)
    for ssc in arange(len(slow_storage_capacity_array)):
        balancing_energies[ssc], optimal_alphas[ssc] = find_balancing_at_optimal_alpha(gamma, fast_eta_in, fast_eta_out, slow_eta_in, slow_eta_out, fast_storage_capacity, slow_storage_capacity_array[ssc], cut_time, returnalpha = True)
        mismatch, smallstorage_usage, used_storages[ssc] = get_mismatch_after_double_storage(get_mismatch(gamma, optimal_alphas[ssc]), fast_eta_in, fast_eta_out, slow_eta_in, slow_eta_out, fast_storage_capacity, slow_storage_capacity_array[ssc])

    used_storages = 370. * used_storages
    slow_storage_capacity_array = 370. * slow_storage_capacity_array 
    balancing_energies = balancing_energies / 70128.
        
    #print optimal_alphas

    plot(used_storages, balancing_energies)
    for i in arange(len(slow_storage_capacity_array)):
        text(used_storages[i], balancing_energies[i], "$alpha = %.2f$"%optimal_alphas[i])
    
def plot_fast_storage_vs_balance_graph(gamma, fast_eta_in = .9, fast_eta_out = .9, slow_eta_in = .6, slow_eta_out = .6, fast_storage_capacity_array = array([0.,3.,6.,9.,12.,24.]), slow_storage_capacity = 200., cut_time = 24.):
    used_storages = empty_like(fast_storage_capacity_array)
    optimal_alphas = empty_like(fast_storage_capacity_array)
    balancing_energies = empty_like(fast_storage_capacity_array)
    for fsc in arange(len(fast_storage_capacity_array)):
        balancing_energies[fsc], optimal_alphas[fsc] = find_balancing_at_optimal_alpha(gamma, fast_eta_in, fast_eta_out, slow_eta_in, slow_eta_out, fast_storage_capacity_array[fsc], slow_storage_capacity, cut_time, returnalpha = True)
        mismatch, used_storages[fsc], largestorage_usage = get_mismatch_after_double_storage(get_mismatch(gamma, optimal_alphas[fsc]), fast_eta_in, fast_eta_out, slow_eta_in, slow_eta_out, fast_storage_capacity_array[fsc], slow_storage_capacity)

    used_storages = 370. * used_storages
    fast_storage_capacity_array = 370. * fast_storage_capacity_array 
    balancing_energies = balancing_energies / 70128.
        
    #print optimal_alphas

    plot(used_storages, balancing_energies)
    for i in arange(len(fast_storage_capacity_array)):
        text(used_storages[i], balancing_energies[i], "$alpha = %.2f$"%optimal_alphas[i])
    
def get_logged_scale(N = 25, M = 1.):
    #partition: 1 + N/5 + R
    logpart = N / 5 + 1
    rest = N - logpart
    return M * concatenate([[0], logspace(-2 - log(N / 25.), log10(1), logpart)[:-1] / 10, linspace(.1, 1, rest)])

def get_balance_array_from_ssc(gamma, ssc_array, fsc = 6, fast_eta_in = .9, fast_eta_out = .9, slow_eta_in = .6, slow_eta_out = .6, cut_time = 24, slow_storage_speed = NaN, returnalphas = True):
    alphas = empty_like(ssc_array)
    balancings = empty_like(ssc_array)
    for x in arange(len(ssc_array)):
        balancings[x], alphas[x] = find_balancing_at_optimal_alpha(gamma, fast_eta_in, fast_eta_out, slow_eta_in, slow_eta_out, fsc, ssc_array[x], cut_time, returnalpha = True, slow_storage_speed = slow_storage_speed)
    if returnalphas:
        return balancings, alphas
    else:
        return balancings

def get_storage_balance_and_alpha_arrays(gamma, N = 25, fast_storage_capacity = 6., fast_eta_in = .9, fast_eta_out = .9, slow_eta_in = .6, slow_eta_out = .6, cut_time = 24, slow_storage_speed = NaN):
    #find optimal alpha for unlimited slow storage:
    t1, optimal_alpha = find_balancing_at_optimal_alpha(gamma, fast_eta_in, fast_eta_out, slow_eta_in, slow_eta_out, fast_storage_capacity, NaN, cut_time, returnalpha = True, slow_storage_speed = slow_storage_speed)
    print "BALANCING: ", t1
    sys.stdout.flush()
    t1, t2, M = get_mismatch_after_double_storage(get_mismatch(gamma, optimal_alpha), fast_eta_in, fast_eta_out, slow_eta_in, slow_eta_out, fast_storage_capacity, slow_storage_capacity = NaN, cut_time = cut_time, slow_storage_speed = slow_storage_speed)
    print "NEEDED STORAGE: ", M
    sys.stdout.flush()
    X = get_logged_scale(N, M)
    print "X: ", X
    sys.stdout.flush()
    Y, alphas = get_balance_array_from_ssc(gamma, X, fast_storage_capacity, fast_eta_in, fast_eta_out, slow_eta_in, slow_eta_out, cut_time, slow_storage_speed = slow_storage_speed, returnalphas = True)
    return X, Y, alphas

def get_countrydata(gamma = 1.2, N = 50, fast_storage_capacity = 6., fast_eta_in = .9, fast_eta_out = .9, slow_eta_in = .6, slow_eta_out = .6, cut_time = 24., slow_storage_speed = NaN):
    countrydata = []
    for i in arange(len(countries.Countries)):
        countries.Region = countries.Countries[i][1]
        print countries.Countries[i][0]
        X, Y, ALPHAS = get_storage_balance_and_alpha_arrays(gamma, N, fast_storage_capacity, fast_eta_in, fast_eta_out, slow_eta_in, slow_eta_out, cut_time, slow_storage_speed)
        countrydata.append([countries.Countries[i][0], get_characteristic_data_of_convex_curve(X,Y), X, Y, ALPHAS])
    return countrydata

def get_curve_for_fixed_alpha(gamma = 1., alpha = .6, N = 100, fast_eta_in = .9, fast_eta_out = .9, slow_eta_in = .6, slow_eta_out = .6, fast_storage_capacity = 0, cut_time = 0, slow_storage_speed = NaN):
    mismatch = get_mismatch(gamma, alpha)
    t1, t2, M = get_mismatch_after_double_storage(mismatch, fast_eta_in, fast_eta_out, slow_eta_in, slow_eta_out, fast_storage_capacity, slow_storage_capacity = NaN, cut_time = cut_time, slow_storage_speed = slow_storage_speed)
    print "NEEDED STORAGE: ", M
    X = get_logged_scale(N, M)
    Y = empty(N)
    for i in arange(N):
        Y[i] = get_balancing_from_mismatch(get_mismatch_after_double_storage(mismatch, fast_eta_in, fast_eta_out, slow_eta_in, slow_eta_out, fast_storage_capacity, X[i], cut_time, True, slow_storage_speed))
    return X, Y

def get_initial_slopes(countrydata):
    inislopes = []
    for i in arange(len(countrydata)):
        inislopes.append((countrydata[i][1])[2])
    return inislopes

def get_final_slopes(countrydata):
    finalslopes = []
    for i in arange(len(countrydata)):
        finalslopes.append((countrydata[i][1])[3])
    return finalslopes

def get_integrals(countrydata):
    integrals = []
    for i in arange(len(countrydata)):
        integrals.append((countrydata[i][1])[4])
    return integrals

def get_normalized_balancings(countrydata):
    balancings = []
    for i in arange(len(countrydata)):
        balancings.append((countrydata[i][1])[0])
    return balancings

def get_normalized_storages(countrydata):
    storages = []
    for i in arange(len(countrydata)):
        storages.append((countrydata[i][1])[1])
    return storages

def get_base_balancings(countrydata):
    balancings = []
    for i in arange(len(countrydata)):
        balancings.append((countrydata[i][1])[5])
    return balancings

def get_country_names(countrydata):
    names = []
    for i in arange(len(countrydata)):
        names.append(countrydata[i][0])
    return names

def plot_slopes(countrydata):
    N = len(countrydata)
    inislopes = get_initial_slopes(countrydata)
    finalslopes = get_final_slopes(countrydata)
    integrals = get_integrals(countrydata)
    names = get_country_names(countrydata)

    f = open("plot_slopes.txt", 'w')
    print >>f, inislopes, finalslopes, integrals, names
    f.close()
    
    ind = arange(N)
    width = 0.27

    subplot(111)
    rects1 = bar(ind, inislopes, width, color = 'r')
    rects2 = bar(ind + width, integrals, width, color = 'g')
    rects3 = bar(ind + 2 * width, finalslopes, width, color = 'b')

    xticks(ind + 1.5*width, names, rotation = 90)

    legend((rects1[0], rects2[0], rects3[0]), ("Normalized initial slopes", "Normalized integrals", "Normalized final slopes"), loc = 9)

def plot_balances_storages_bases(countrydata):
    N = len(countrydata)
    balancings = get_normalized_balancings(countrydata)
    storages = get_normalized_storages(countrydata)
    bases = get_base_balancings(countrydata)
    names = get_country_names(countrydata)

    f = open("plot_balances.txt", 'w')
    print >>f, balancings, storages, bases, names
    f.close()
    
    ind = arange(N)
    width = 0.27

    subplot(111)
    rects1 = bar(ind, balancings, width, color = 'r')
    rects2 = bar(ind + width, storages, width, color = 'g')
    rects3 = bar(ind + 2 * width, bases, width, color = 'b')

    xticks(ind + 1.5*width, names, rotation = 90)

    legend((rects1[0], rects2[0], rects3[0]), ("Balancings", "Storages", "Base balancings"), loc = 2)

def get_alphas_from_countrydata(countrydata):
    alphas = []
    for i in arange(len(countrydata)):
        alphas.append(countrydata[i][4])
    return alphas

def get_X_from_countrydata(countrydata):
    X = []
    for i in arange(len(countrydata)):
        X.append(countrydata[i][2])
    return X

def get_Y_from_countrydata(countrydata):
    Y = []
    for i in arange(len(countrydata)):
        Y.append(countrydata[i][3])
    return Y

def plot_alphas(cd):
    countrydata = deepcopy(cd)
    N = len(countrydata)
    storages = 8766*array(get_normalized_storages(countrydata))
    Xes = get_X_from_countrydata(countrydata)
    alphas = get_alphas_from_countrydata(countrydata)

    for i in arange(N):
        new_alpha = interp(storages[i], Xes[i], alphas[i])
        alphas[i][Xes[i] > storages[i]] = new_alpha
        Xes[i][Xes[i] > storages[i]] = storages[i]

    f = open("plot_alphas.txt", 'w')
    print >>f, N, storages, Xes, alphas

    for i in arange(N - 1):
        plot(Xes[i] / storages[i], alphas[i])
    plot(Xes[N - 1] / storages[N - 1], alphas[N - 1], lw = 4)

def plot_alphas_unnormalized(cd):
    countrydata = deepcopy(cd)
    N = len(countrydata)
    storages = 8766*array(get_normalized_storages(countrydata))
    Xes = get_X_from_countrydata(countrydata)
    alphas = get_alphas_from_countrydata(countrydata)

    for i in arange(N):
        new_alpha = interp(storages[i], Xes[i], alphas[i])
        alphas[i][Xes[i] > storages[i]] = new_alpha
        Xes[i][Xes[i] > storages[i]] = storages[i]

    f = open("plot_alphas.txt", 'w')
    print >>f, N, storages, Xes, alphas

    for i in arange(N - 1):
        plot(Xes[i], alphas[i])
    plot(Xes[N - 1], alphas[N - 1], lw = 4)

def plot_balances(cd):
    countrydata = deepcopy(cd)
    N = len(countrydata)
    storages = 8766*array(get_normalized_storages(countrydata))
    Xes = get_X_from_countrydata(countrydata)
    Ys = get_Y_from_countrydata(countrydata)
    print storages
    print Xes[1]
    print Ys[1]

    for i in arange(N):
        new_balance = interp(storages[i], Xes[i], Ys[i])
        print new_balance
        Ys[i][Xes[i] > storages[i]] = new_balance
        Xes[i][Xes[i] > storages[i]] = storages[i]

    f = open("plot_balances.txt", 'w')
    print >>f, N, storages, Xes, Ys

    for i in arange(N - 1):
        plot(Xes[i] / storages[i], Ys[i] / Ys[i][0])
    plot(Xes[N - 1] / storages[N - 1], Ys[N - 1] / Ys[N - 1][0], lw = 4)

def plot_balances_unnormalized(cd):
    countrydata = deepcopy(cd)
    N = len(countrydata)
    storages = 8766*array(get_normalized_storages(countrydata))
    Xes = get_X_from_countrydata(countrydata)
    Ys = get_Y_from_countrydata(countrydata)
    print storages
    print Xes[1]
    print Ys[1]

    for i in arange(N):
        new_balance = interp(storages[i], Xes[i], Ys[i])
        print new_balance
        Ys[i][Xes[i] > storages[i]] = new_balance
        Xes[i][Xes[i] > storages[i]] = storages[i]

    f = open("plot_balances.txt", 'w')
    print >>f, N, storages, Xes, Ys

    for i in arange(N - 1):
        plot(Xes[i], Ys[i])
    plot(Xes[N - 1], Ys[N - 1], lw = 4)

def get_storage_balance_power_arrays(gamma, X, alphas, fast_storage_capacity = 6., fast_eta_in = .9, fast_eta_out = .9, slow_eta_in = .6, slow_eta_out = .6, cut_time = 24, slow_storage_speed = NaN):
    N = len(X)
    storage_in_power = empty(N)
    storage_out_power = empty(N)
    balance_power = empty(N)
    for i in arange(N):
        mismatch = get_mismatch(gamma, alphas[i])
        new_mismatch = get_mismatch_after_double_storage(mismatch, fast_eta_in, fast_eta_out, slow_eta_in, slow_eta_out, fast_storage_capacity, slow_storage_capacity = X[i], cut_time = cut_time, returnmismatchonly = True, slow_storage_speed = slow_storage_speed)
        balance_power[i] = (-new_mismatch).max()
        storage_in_power[i] = (mismatch - new_mismatch).max()
        storage_out_power[i] = (new_mismatch - mismatch).max()
    return storage_in_power, storage_out_power, balance_power

def get_storage_balance_quantile_power_arrays(gamma, X, alphas, fast_storage_capacity = 6., fast_eta_in = .9, fast_eta_out = .9, slow_eta_in = .6, slow_eta_out = .6, cut_time = 24, slow_storage_speed = NaN):
    N = len(X)
    storage_in_power = empty(N)
    storage_out_power = empty(N)
    balance_power = empty(N)
    for i in arange(N):
        mismatch = get_mismatch(gamma, alphas[i])
        new_mismatch = get_mismatch_after_double_storage(mismatch, fast_eta_in, fast_eta_out, slow_eta_in, slow_eta_out, fast_storage_capacity, slow_storage_capacity = X[i], cut_time = cut_time, returnmismatchonly = True, slow_storage_speed = slow_storage_speed)
        balance_power[i] = double(mquantiles(-new_mismatch, .99))
        storage_in_power[i] = double(mquantiles(mismatch - new_mismatch, .99))
        storage_out_power[i] = double(mquantiles(new_mismatch - mismatch, .99))
    return storage_in_power, storage_out_power, balance_power

def plot_many_storage_fillings(gamma, alpha, N, eta_in = 1., eta_out = 1., policy = 3):
    storage_capacity = plot_storage_filling_etc(gamma, alpha, NaN, eta_in, eta_out, NaN, policy)
    a = linspace(0, storage_capacity, N)[:-1]
    print a
    for x in a:
        plot_storage_filling_etc(gamma, alpha, x, eta_in, eta_out, storage_capacity, policy)

def plot_smallest_ten_percent(gamma, alpha, N, eta_in = 1., eta_out = 1., policy = 3):
    storage_capacity = plot_storage_filling_etc(gamma, alpha, NaN, eta_in, eta_out, NaN, policy)
    a = linspace(0, storage_capacity/10, N + 1)
    print a
    for x in a:
        plot_storage_filling_etc(gamma, alpha, x, eta_in, eta_out, storage_capacity/10, policy)

def plot_storage_filling_etc(gamma, alpha, storage_capacity, eta_in = 1., eta_out = 1., maxvalue = NaN, policy = 3):
    mismatch = get_mismatch(gamma, alpha)
    if policy == 4:
        new_mismatch = get_mismatch_after_policy_4_storage(mismatch, eta_in, eta_out, storage_capacity)
    elif policy == 2:
        new_mismatch = get_mismatch_after_policy_2_storage(mismatch, eta_in, eta_out, storage_capacity)
    elif policy == 1:
        new_mismatch = get_mismatch_after_policy_1_storage(mismatch, eta_in, eta_out, storage_capacity)
    elif policy == 5:
        new_mismatch = get_mismatch_after_policy_5_storage(mismatch, eta_in, eta_out, storage_capacity)
    else:
        new_mismatch = get_mismatch_after_policy_3_storage(mismatch, eta_in, eta_out, storage_capacity)
    storage_transactions = mismatch - new_mismatch

    
    storage_filling = eta_in * pos(storage_transactions).cumsum() - neg(storage_transactions).cumsum() / eta_out
    storage_filling = storage_filling - storage_filling.min()
    storage_capacity = storage_filling.max()

    if isnan(maxvalue):
        maxvalue = storage_capacity

    balancing = neg(new_mismatch)
    
    excess = pos(new_mismatch)

    if storage_capacity != 0:
        normfactor = maxvalue / 2
    else:
        normfactor = 1.
    
    clf()
    subplot(311)
    title(r"$\gamma = %.2f, \alpha_w = %.2f, C_S = %.2f, \eta_{\mathrm{in}}=%.2f, \eta_{\mathrm{out}}=%.2f, \langle E_{\mathrm{xs}}\rangle=%.3f, \langle E_{\mathrm{xs}}\rangle/Q_{\mathrm{xs}}=%.3f$"%(gamma, alpha, storage_capacity, eta_in, eta_out, excess.mean(), excess.mean() / float(mquantiles(excess, .99))) + "\n" + r"${\langle E_{S,\mathrm{in}}\rangle}/{Q_{S,\mathrm{in}}}=%.4f, Q_B = %.3f, {\langle E_B\rangle}/{Q_B}=%.3f, {\langle S_L\rangle}/{C_S}=%.3f, {\langle \Delta {S_{L,+}}\rangle}/{C_S}=%.3g$"%(pos(storage_transactions).mean() / float(mquantiles(pos(storage_transactions), .99)), float(mquantiles(balancing, .99)), balancing.mean() / float(mquantiles(balancing, .99)), storage_filling.mean() / storage_capacity, (pos(storage_transactions) * eta_in).mean() * 8766 / storage_capacity))
    plot(storage_filling, lw = 1, color = "b", label = "storage filling")
    axis(ymax = maxvalue, ymin = 0, xmax = 70128)
    subplot(312)
    plot(storage_transactions, color = 'g', label = "storage transactions", alpha = .3)
    axis(ymax = max(mismatch), ymin = min(mismatch), xmax = 70128)
#    legend(loc = 4)
    subplot(313)
    plot(balancing, lw = 1, color = 'r', label = "balancing")
    plot(excess, lw = 1, label = "excess", color = 'y', alpha = .7)

    axis(xmax = 70128)
    legend(loc = 1)
    print r"Policy_%d_gamma%.2falpha_w%.2fC_S%04d.%02deta_in%.2feta_out%.2f.pdf"%(policy, gamma, alpha, storage_capacity, 100 * (storage_capacity - int(storage_capacity)), eta_in, eta_out)
    savefig(r"Policy_%d_gamma%.2falpha_w%.2fC_S%04d.%02deta_in%.2feta_out%.2f.pdf"%(policy, gamma, alpha, storage_capacity, 100 * (storage_capacity - int(storage_capacity)), eta_in, eta_out))

    return storage_capacity

def get_storage_data(gamma, alpha, storage_capacity, eta_in = 1., eta_out = 1., policy = 3, quantile = .99, mismatch = None):
    if not isnan(gamma):
        mismatch = get_mismatch(gamma, alpha)
    if mismatch.max() < 0 or mismatch.min() > 0:
        print "Warning!"
        new_mismatch = mismatch
    elif policy == 4:
        new_mismatch = get_mismatch_after_policy_4_storage(mismatch, eta_in, eta_out, storage_capacity)
    elif policy == 2:
        new_mismatch = get_mismatch_after_policy_2_storage(mismatch, eta_in, eta_out, storage_capacity)
    elif policy == 1:
        new_mismatch = get_mismatch_after_policy_1_storage(mismatch, eta_in, eta_out, storage_capacity)
    elif policy == 5:
        new_mismatch = get_mismatch_after_policy_5_storage(mismatch, eta_in, eta_out, storage_capacity)
    elif policy == 7:
        new_mismatch = get_mismatch_after_policy_7_storage(mismatch, eta_in, eta_out, storage_capacity)
    elif policy == 8:
        new_mismatch = get_mismatch_after_policy_8_storage(mismatch, eta_in, eta_out, storage_capacity)
    elif policy == 9:
        new_mismatch = get_mismatch_after_policy_9_storage(mismatch, eta_in, eta_out, storage_capacity)
    else:
        new_mismatch = get_mismatch_after_policy_3_storage(mismatch, eta_in, eta_out, storage_capacity)
    storage_transactions = mismatch - new_mismatch
    
    storage_filling = eta_in * pos(storage_transactions).cumsum() - neg(storage_transactions).cumsum() / eta_out
    storage_filling = storage_filling - storage_filling.min()

    storage_capacity = storage_filling.max()

    balancing = neg(new_mismatch)
    if policy == 8:
        balancing = balancing * (balancing > constant.policy_8_capacity) + constant.policy_8_capacity * (balancing <= constant.policy_8_capacity)
    
    #plot(storage_filling)

    excess = pos(new_mismatch)
    return storage_capacity, excess.mean(), mquantiles(excess, quantile), pos(storage_transactions).mean(), mquantiles(pos(storage_transactions), quantile), neg(storage_transactions).mean(), mquantiles(neg(storage_transactions), quantile), balancing.mean(), mquantiles(balancing, quantile), storage_filling.mean()


                                                                                
def plot_storages_data(gamma, alpha, N = 30, eta_in = 1., eta_out = 1., policy = 1):
    storages, mean_excesses, Q_excesses, mean_store_ins, Q_store_ins, mean_store_outs, Q_store_outs, mean_Bs, Q_Bs, mean_Ses = get_storages_data(gamma, alpha, N, eta_in, eta_out, policy)
    figure(1)
    clf()
    plot(storages, mean_excesses, color = 'b', label = "Mean excess")
    plot(storages, Q_excesses, '--', color = 'b', label = "99% quantile excess")
#    plot(storages, mean_excesses / Q_excesses, '-.', color = 'b')
    plot(storages, mean_store_ins, color = 'y', label = "Mean $E_{S,\mathrm{in}}$")
    plot(storages, Q_store_ins, "--", color = 'y', label = "99% quantile $E_{S,\mathrm{in}}$")
#    plot(storages, mean_store_ins / Q_store_ins, '-.', color = 'y')
    plot(storages, mean_store_outs, color = 'g', label = "Mean $E_{S,\mathrm{out}}$")
    plot(storages, Q_store_outs, "--", color = 'g', label = "99% quantile $E_{S,\mathrm{out}}$")
#    plot(storages, mean_store_outs / Q_store_outs, '-.', color = 'g')
    plot(storages, mean_Bs, color = 'r', label = "Mean $E_B$")
    plot(storages, Q_Bs, "--", color = 'r', label = "99% quantile $E_B$")
#    plot(storages, mean_Bs / Q_Bs, '-.', color = 'r')
    legend()
    savefig("storage_data_policy%d_gamma%.2f_alpha_w%.2f.pdf"%(policy, gamma, alpha))
    return storages, mean_Ses

def plot_policies_data(gamma, alpha, N = 30, eta_in = 1., eta_out = 1.):
    p1storages, p1means = plot_storages_data(gamma, alpha, N, eta_in, eta_out, policy = 1)
    p2storages, p2means = plot_storages_data(gamma, alpha, N, eta_in, eta_out, policy = 2)
    p3storages, p3means = plot_storages_data(gamma, alpha, N, eta_in, eta_out, policy = 3)
    p4storages, p4means = plot_storages_data(gamma, alpha, N, eta_in, eta_out, policy = 4)
    p5storages, p5means = plot_storages_data(gamma, alpha, N, eta_in, eta_out, policy = 5)
    figure(2)
    clf()
    title("Mean storage level ratio as function of storage capacity")
    plot(p1storages, p1means / p1storages, label = "Policy 1")
    plot(p2storages, p2means / p2storages, label = "Policy 2")
    plot(p3storages, p3means / p3storages, label = "Policy 3")
    plot(p4storages, p4means / p4storages, label = "Policy 4")
    plot(p5storages, p5means / p5storages, label = "Policy 5")
    legend()
    savefig("storage_policies_gamma%.2f_alpha_w%.2f.pdf"%(gamma, alpha))

def get_storages_data(gamma, alpha, N = 30, eta_in = 1., eta_out = 1., policy = 3, quantile = .99, Narray = None, mismatch = None):
    quantiles = array(quantile)
    L = quantiles.size
    mean_excesses = empty(N)
    Q_excesses = empty([N, L])
    mean_store_ins = empty(N)
    Q_store_ins = empty([N, L])
    mean_store_outs = empty(N)
    Q_store_outs = empty([N, L])
    mean_Bs = empty(N)
    Q_Bs = empty([N, L])
    mean_Ses = empty(N)
    storage_capacity, mean_excesses[N - 1], Q_excesses[N - 1], mean_store_ins[N - 1], Q_store_ins[N - 1], mean_store_outs[N - 1], Q_store_outs[N - 1], mean_Bs[N - 1], Q_Bs[N - 1], mean_Ses[N - 1] = get_storage_data(gamma, alpha, NaN, eta_in, eta_out, policy, quantile, mismatch)
    if Narray is None:
        storages = storage_capacity * exp(linspace(-N/10.,N/10.,N))/exp(N/10.)
    else:
        storages = storage_capacity * Narray
    print storages
    for i in arange(N - 1):
        storages[i], mean_excesses[i], Q_excesses[i], mean_store_ins[i], Q_store_ins[i], mean_store_outs[i], Q_store_outs[i], mean_Bs[i], Q_Bs[i], mean_Ses[i] = get_storage_data(gamma, alpha, storages[i], eta_in, eta_out, policy, quantile, mismatch)
    return storages, mean_excesses, Q_excesses, mean_store_ins, Q_store_ins, mean_store_outs, Q_store_outs, mean_Bs, Q_Bs, mean_Ses

def get_no_storage_data(alpha, eta_in = 1., eta_out = 1., N = 31, policy = 3, quantile = .99):
    mean_excess = empty(N)
    Q_excess = empty(N)
    mean_balancing = empty(N)
    Q_balancing = empty(N)
    gamma = linspace(.5, 1.5, N)
    for i in arange(N):
        a, mean_excess[i], Q_excess[i], b, c, d, e, mean_balancing[i], Q_balancing[i], f = get_storage_data(gamma[i], alpha, 0., eta_in, eta_out, policy, quantile)
    return mean_excess, Q_excess, mean_balancing, Q_balancing

def plot_no_storage_data(alpha, eta_in = 1., eta_out = 1., N = 31, policy = 2):
    mean_excess, Q9_excess, mean_balancing, Q9_balancing = get_no_storage_data(alpha, eta_in, eta_out, N, policy, .9)
    mean_excess, Q99_excess, mean_balancing, Q99_balancing = get_no_storage_data(alpha, eta_in, eta_out, N, policy, .99)
    mean_excess, Q999_excess, mean_balancing, Q999_balancing = get_no_storage_data(alpha, eta_in, eta_out, N, policy, .999)
    mean_excess, Q9999_excess, mean_balancing, Q9999_balancing = get_no_storage_data(alpha, eta_in, eta_out, N, policy, .9999)
    mean_excess, max_excess, mean_balancing, max_balancing = get_no_storage_data(alpha, eta_in, eta_out, N, policy, 1.)
    gamma = linspace(.5, 1.5, N)
    clf()
    title(r"Functions of $\gamma$ in a no-storage scenario" + "\n" + r"$\alpha = %.2f$"%(alpha))
    plot(gamma, mean_excess, "b", lw = 2, label = "Mean excess")
    plot(gamma, Q9_excess, "b-", label = ".9 quantile")
    plot(gamma, Q99_excess, "b--", label = ".99 quantile")
    plot(gamma, Q999_excess, "b-.", label = ".999 quantile")
    plot(gamma, Q9999_excess, "b:", label = ".9999 quantile")
    plot(gamma, max_excess, "b.", label = "max excess")
    plot(gamma, mean_balancing, "r", lw = 2, label = "Mean balancing")
    plot(gamma, Q9_balancing, "r-", label = ".9 quantile")
    plot(gamma, Q99_balancing, "r--", label = ".99 quantile")
    plot(gamma, Q999_balancing, "r-.", label = ".999 quantile")
    plot(gamma, Q9999_balancing, "r:", label = ".9999 quantile")
    plot(gamma, max_balancing, "r.", label = "max balancing")
    legend(loc = 2)
    axis("tight")
    xlabel(r"$\gamma$")
    ylabel("units of average hourly load")


def get_no_balancing_data(alpha, eta_in = 1., eta_out = 1., N = 31, policy = 2, quantile = .99):
    quantiles = array(quantile)
    L = quantiles.size
    storage_capacity = empty(N)
    mean_excess = empty(N)
    Q_excess = empty([N, L])
    mean_balancing = empty(N)
    Q_balancing = empty([N, L])
    gamma = linspace(0.5, 1.5, N)
    for i in arange(N):
        storage_capacity[i], mean_excess[i], Q_excess[i], b, c, d, e, mean_balancing[i], Q_balancing[i], f = get_storage_data(gamma[i], alpha, NaN, eta_in, eta_out, policy, quantile)
    return storage_capacity, mean_excess, Q_excess, mean_balancing, Q_balancing
    
def get_daily_no_balancing_data(alpha, eta_in = 1., eta_out = 1., N = 31, policy = 2, quantile = .99):
    quantiles = array(quantile)
    L = quantiles.size
    storage_capacity = empty(N)
    mean_excess = empty(N)
    Q_excess = empty([N, L])
    mean_balancing = empty(N)
    Q_balancing = empty([N, L])
    gamma = linspace(0., 2., N)
    for i in arange(N):
        daily_mismatch = get_time_sum(get_mismatch(gamma[i], alpha))
        storage_capacity[i], mean_excess[i], Q_excess[i], b, c, d, e, mean_balancing[i], Q_balancing[i], f = get_storage_data(NaN, alpha, NaN, eta_in, eta_out, policy, quantile, daily_mismatch)
    return storage_capacity, mean_excess, Q_excess, mean_balancing, Q_balancing

def get_daily_bounded_storage_data(alpha, storage_bound, eta_in = 1., eta_out = 1., N = 31, policy = 2, quantile = .99):
    quantiles = array(quantile)
    L = quantiles.size
    storage_capacity = empty(N)
    mean_excess = empty(N)
    Q_excess = empty([N, L])
    mean_balancing = empty(N)
    Q_balancing = empty([N, L])
    gamma = linspace(0., 2., N)
    for i in arange(N):
        daily_mismatch = get_time_sum(get_mismatch(gamma[i], alpha))
        storage_capacity[i], mean_excess[i], Q_excess[i], b, c, d, e, mean_balancing[i], Q_balancing[i], f = get_storage_data(NaN, alpha, storage_bound, eta_in, eta_out, policy, quantile, daily_mismatch)
    return storage_capacity, mean_excess, Q_excess, mean_balancing, Q_balancing


def plot_no_balancing_data(alpha, eta_in = 1., eta_out = 1., N = 31, policy = 3):
    storage_capacity, mean_excess, Q99_excess, mean_balancing, Q99_balancing = get_no_balancing_data(alpha, eta_in, eta_out, N, policy, .99)
    gamma = linspace(.5, 1.5, N)
    clf()
    title(r"Functions of $\gamma$ in a no-constraints-on-storage scenario" + "\n" + r"$\alpha = %.2f,\, \eta_{\mathrm{in}} = %.2f,\, \eta_{\mathrm{out}} = %.2f$"%(alpha, eta_in, eta_out))
    plot(gamma, storage_capacity * max(Q99_balancing.max(), Q99_excess.max()) / storage_capacity.max(), "g", lw = 3, label = "Storage capacity\n(scaled $%4e$)"%(max(Q99_balancing.max(), Q99_excess.max()) / storage_capacity.max()))
    plot(gamma, mean_excess, "b", lw = 2, label = "Mean excess")
    plot(gamma, Q99_excess, "b--", label = ".99 quantile")
    plot(gamma, mean_balancing, "r", lw = 2, label = "Mean balancing")
    plot(gamma, Q99_balancing, "r--", label = ".99 quantile")
    legend(loc = 2)
    axis("tight")
    xlabel(r"$\gamma$")
    ylabel("units of average hourly load")

def plot_E_B_vs_C_S(gamma, alpha, N = 30, eta_in = 1., eta_out = 1., policy = 3):
    storages, mean_excesses, Q_excesses, mean_store_ins, Q_store_ins, mean_store_outs, Q_store_outs, mean_Bs, Q_Bs, mean_Ses = get_storages_data(gamma, alpha, N, eta_in, eta_out, policy, [.9, .99, .999, .9999, 1.])
    plot(storages, mean_Bs, "r", lw = 2, label = "Average balancing")

def plot_Q_Bs_vs_C_S(gamma, alpha, N = 30, eta_in = 1., eta_out = 1., policy = 3):
    storages, mean_excesses, Q_excesses, mean_store_ins, Q_store_ins, mean_store_outs, Q_store_outs, mean_Bs, Q_Bs, mean_Ses = get_storages_data(gamma, alpha, N, eta_in, eta_out, policy, [.9, .99, .999, .9999, 1.], linspace(0, 1, N))
    plot(storages, Q_Bs[:,0], "r-", label = ".9 quantile")
    plot(storages, Q_Bs[:,1], "r--", label = ".99 quantile")
    plot(storages, Q_Bs[:,2], "r-.", label = ".999 quantile")
    plot(storages, Q_Bs[:,3], "r:", label = ".9999 quantile")
    plot(storages, Q_Bs[:,4], "r.", label = "max balancing")
    legend()

def plot_20110423_figure1(policy = 3):
    clf()
    plot_E_B_vs_C_S(1,.6, 50, policy = policy)
    plot_Q_Bs_vs_C_S(1,.6, 301, policy = policy)
    title("Balancing as a function of storage energy capacity\n" + r"$\gamma = 1.00,\,\alpha_W = 0.60,\,\eta_{\mathrm{in}}=1.00,\,\eta_{\mathrm{out}}=1.00$, policy %d"%(policy))
    savefig("Martin20110423Figure1policy%d.pdf"%(policy))

def plot_mean_store_in_vs_C_S(gamma, alpha, N = 30, eta_in = 1., eta_out = 1., policy = 3):
    storages, mean_excesses, Q_excesses, mean_store_ins, Q_store_ins, mean_store_outs, Q_store_outs, mean_Bs, Q_Bs, mean_Ses = get_storages_data(gamma, alpha, N, eta_in, eta_out, policy, [.9, .99, .999, .9999, 1.])
    plot(storages, mean_store_ins, "b", lw = 2, label = "Average $S_+$")

def plot_Q_store_ins_vs_C_S(gamma, alpha, N = 30, eta_in = 1., eta_out = 1., policy = 3):
    storages, mean_excesses, Q_excesses, mean_store_ins, Q_store_ins, mean_store_outs, Q_store_outs, mean_Bs, Q_Bs, mean_Ses = get_storages_data(gamma, alpha, N, eta_in, eta_out, policy, [.9, .99, .999, .9999, 1.], linspace(0, 1, N))
    plot(storages, Q_store_ins[:,0], "b-", label = ".9 quantile")
    plot(storages, Q_store_ins[:,1], "b--", label = ".99 quantile")
    plot(storages, Q_store_ins[:,2], "b-.", label = ".999 quantile")
    plot(storages, Q_store_ins[:,3], "b:", label = ".9999 quantile")
    plot(storages, Q_store_ins[:,4], "b.", label = r"$\mathrm{max}S_+$")
    legend()

def plot_20110423_figure2(policy = 3):
    clf()
    plot_mean_store_in_vs_C_S(1,.6, 50, policy = policy)
    plot_Q_store_ins_vs_C_S(1,.6, 301, policy = policy)
    title("Storage input as a function of storage energy capacity\n" + r"$\gamma = 1.00,\,\alpha_W = 0.60,\,\eta_{\mathrm{in}}=1.00,\,\eta_{\mathrm{out}}=1.00$, policy %d"%(policy))
    savefig("Martin20110423Figure2policy%d.pdf"%(policy))

def plot_mean_store_out_vs_C_S(gamma, alpha, N = 30, eta_in = 1., eta_out = 1., policy = 3):
    storages, mean_excesses, Q_excesses, mean_store_ins, Q_store_ins, mean_store_outs, Q_store_outs, mean_Bs, Q_Bs, mean_Ses = get_storages_data(gamma, alpha, N, eta_in, eta_out, policy, [.9, .99, .999, .9999, 1.])
    plot(storages, mean_store_outs, "b", lw = 2, label = "Average $S_-$")

def plot_Q_store_outs_vs_C_S(gamma, alpha, N = 30, eta_in = 1., eta_out = 1., policy = 3):
    storages, mean_excesses, Q_excesses, mean_store_ins, Q_store_ins, mean_store_outs, Q_store_outs, mean_Bs, Q_Bs, mean_Ses = get_storages_data(gamma, alpha, N, eta_in, eta_out, policy, [.9, .99, .999, .9999, 1.], linspace(0, 1, N))
    plot(storages, Q_store_outs[:,0], "b-", label = ".9 quantile")
    plot(storages, Q_store_outs[:,1], "b--", label = ".99 quantile")
    plot(storages, Q_store_outs[:,2], "b-.", label = ".999 quantile")
    plot(storages, Q_store_outs[:,3], "b:", label = ".9999 quantile")
    plot(storages, Q_store_outs[:,4], "b.", label = r"$\mathrm{max}S_-$")
    legend()

def plot_20110423_figure3(policy = 3):
    clf()
    plot_mean_store_out_vs_C_S(1,.6, 50, policy = policy)
    plot_Q_store_outs_vs_C_S(1,.6, 301, policy = policy)
    title("Storage output as a function of storage energy capacity\n" + r"$\gamma = 1.00,\,\alpha_W = 0.60,\,\eta_{\mathrm{in}}=1.00,\,\eta_{\mathrm{out}}=1.00$, policy %d"%(policy))
    savefig("Martin20110423Figure3policy%d.pdf"%(policy))

def plot_mean_excess_vs_C_S(gamma, alpha, N = 30, eta_in = 1., eta_out = 1., policy = 3):
    storages, mean_excesses, Q_excesses, mean_store_ins, Q_store_ins, mean_store_outs, Q_store_outs, mean_Bs, Q_Bs, mean_Ses = get_storages_data(gamma, alpha, N, eta_in, eta_out, policy, [.9, .99, .999, .9999, 1.])
    plot(storages, mean_excesses, "k", lw = 2, label = "Average excess")

def plot_Q_excesses_vs_C_S(gamma, alpha, N = 30, eta_in = 1., eta_out = 1., policy = 3):
    storages, mean_excesses, Q_excesses, mean_store_ins, Q_store_ins, mean_store_outs, Q_store_outs, mean_Bs, Q_Bs, mean_Ses = get_storages_data(gamma, alpha, N, eta_in, eta_out, policy, [.9, .99, .999, .9999, 1.], linspace(0, 1, N))
    plot(storages, Q_excesses[:,0], "k-", label = ".9 quantile")
    plot(storages, Q_excesses[:,1], "k--", label = ".99 quantile")
    plot(storages, Q_excesses[:,2], "k-.", label = ".999 quantile")
    plot(storages, Q_excesses[:,3], "k:", label = ".9999 quantile")
    plot(storages, Q_excesses[:,4], "k.", label = r"$max excess")
    legend()

def plot_20110423_figure4(policy = 3):
    clf()
    plot_mean_excess_vs_C_S(1,.6, 50, policy = policy)
    plot_Q_excesses_vs_C_S(1,.6, 301, policy = policy)
    title("Excess as a function of storage energy capacity\n" + r"$\gamma = 1.00,\,\alpha_W = 0.60,\,\eta_{\mathrm{in}}=1.00,\,\eta_{\mathrm{out}}=1.00$, policy %d"%(policy))
    savefig("Martin20110423Figure4policy%d.pdf"%(policy))

def find_crossing_point(a1,b1,a2,b2,x1,y1,x2,y2):
    print "crossing",
    if not (a1 <= a2 and a2 <= x1 and x1 <= x2):
        print "Weird: ", a1, a2, x1, x2
    if x1 == x2:
        print "x1 = x2!"
        a = 0
        b = y1
    else:
        a = (y1 - y2) / (x1 - x2)
        b = (y2 * x1 - y1 * x2) / (x1 - x2)
    if a1 == a2:
        crossing = a2
    else:
        c = (b1 - b2) / (a1 - a2)
        d = (b2 * a1 - b1 * a2) / (a1 - a2)
    if a - c == 0:
        print "a = c! ", a1, b1, a2, b2, x1, y1, x2, y2
        crossing = (a2 + x1) / 2
    else:
        crossing = (d - b) / (a - c)
    if crossing > x1:
        print a1, b1, a2, b2, x1, y1, x2, y2
        crossing = x1
    sys.stdout.flush()
    print crossing
    return crossing
#    else:
#        print a1, b1, a2, b2, x1, y1, x2, y2
#        return

def get_gamma_storage_balance_data(alpha, N = 31, eta_in = 1., eta_out = 1., policy = 2, M = NaN):
    if isnan(M) or M < 5:
        M = min(max(5,N),31)
    gamma = linspace(.5, 1.5, N)
    gamma_range = -ones([M,N])
    storage = -ones([M,N])
    balancing = -ones([M,N])
    for i in arange(N):
        storage[0, i], a, b, c, d, e, f, balancing[0, i], g, h = get_storage_data(gamma[i], alpha, 0, eta_in, eta_out, policy, 1.)
        storage[M-2, i], a, b, c, d, e, f, balancing[M-2, i], g, h = get_storage_data(gamma[i], alpha, NaN, eta_in, eta_out, policy, 1.)
        storage[1, i], a, b, c, d, e, f, balancing[1, i], g, h = get_storage_data(gamma[i], alpha, max(1e-13, max(2,storage[M - 2, i])**-2 * 2**-M), eta_in, eta_out, policy, 1.)
        crosspoint = find_crossing_point(storage[0, i], balancing[0, i], storage[1, i], balancing[1, i], storage[M - 2, i], balancing[M - 2, i], storage[M - 2, i], balancing[M - 2, i])
        storage[M-3, i], a, b, c, d, e, f, balancing[M-3, i], g, h = get_storage_data(gamma[i], alpha, a * crosspoint + (1 - a) * storage[M - 2, i], eta_in, eta_out, policy, 1.)
    for i in arange(M):
        gamma_range[i,:] = gamma
    storage[M-1,:] = storage[M-2,:].max() * ones(N)
    balancing[M-1,:] = balancing[M-2,:]
    for m in arange(M-5):
        print m
        sys.stdout.flush()
        if m/2 == m/2.:
            for i in arange(N):
                print i
                crosspoint = find_crossing_point(storage[m/2, i], balancing[m/2, i], storage[m/2 + 1, i], balancing[m/2 + 1, i], storage[M - 3 - m/2, i], balancing[M - 3 - m/2, i], storage[M - 2 - m/2, i], balancing[M - 2 - m/2, i])
                a = 1. / (1 + (M - 4 - m) / 2)
                storage[2 + m/2, i], a, b, c, d, e, f, balancing[2 + m/2, i], g, h = get_storage_data(gamma[i], alpha, (1 - a) * crosspoint + a * storage[1 + m/2,i], eta_in, eta_out, policy, 1.)
        if m/2 != m/2.:
            for i in arange(N):
                print i
                crosspoint = find_crossing_point(storage[1 + m/2, i], balancing[1 + m/2, i], storage[2 + m/2, i], balancing[2 + m/2, i], storage[M - 3 - m/2, i], balancing[M - 3 - m/2, i], storage[M - 2 - m/2, i], balancing[M - 2 - m/2, i])
                a = 1. / (1 + (M - 4 - m) / 2)
                storage[M - 4 - m/2, i], a, b, c, d, e, f, balancing[M - 4 - m/2, i], g, h = get_storage_data(gamma[i], alpha, a * crosspoint + (1 - a) * storage[M - 3 - m/2, i], eta_in, eta_out, policy, 1.)

    return gamma_range, storage, balancing


def old_get_gamma_storage_balance_data(alpha, N = 31, eta_in = 1., eta_out = 1., policy = 2):
    max_storage, me, qx, a, b, c, d, mb, qb, e = get_storage_data(1., alpha, NaN, eta_in, eta_out, policy, 1.)
    storage_limit = concatenate([array([0]), max_storage * logspace(- log(N - 1 / 10.), 1., N - 1) / 10.])
    gamma = linspace(.5, 1.5, N)
    storage_range = empty([N,N])
    gamma_range = empty([N,N])
    storage = empty([N,N])
    balancing = empty([N,N])
    for i in arange(N):
        print "storage_limit: ", storage_limit[i]
        for j in N - 1 - arange(N):
            print "gamma: ", gamma[j]
            storage_range[i,j] = storage_limit[i]
            gamma_range[i,j] = gamma[j]
#            if i + 1 < N:
#                if storage_range[i, j] > storage[i + 1, j]:
#                    storage[i, j] = storage[i + 1, j]
#                    balancing[i, j] = balancing[i + 1, j]
#                else:
#                    storage[i, j], a, b, c, d, e, f, balancing[i, j], g, h = get_storage_data(gamma[j], alpha, storage_limit[i], eta_in, eta_out, policy, 1.)
#            else:
            storage[i, j], a, b, c, d, e, f, balancing[i, j], g, h = get_storage_data(gamma[j], alpha, storage_limit[i], eta_in, eta_out, policy, 1.)
    return gamma_range, storage_range, balancing, storage

def plot_3d_test(x, y, z):
    fig = figure()
    ax = Axes3D(fig)
    ax.scatter3D(z, y, z)

def plot_contour_test(N = 4):
    X, Y, Z1, Z2 = get_gamma_storage_balance_data(.6, N)
    print X, Y, Z1, Z2
    figure()
    contourf(X, Y, Z1)
#    ax.plot_wireframe(X, Y, Z1, alpha = .5, color = 'r')
#    ax.plot_surface(X, Y, Z2/Z2.max(), alpha = .8, color = 'b')
#    ax.plot_wireframe(X, Y, Z2/Z2.max())

def plot_x_y_z(x, y, z):
    fig = figure()
    ax = Axes3D(fig)
    ax.plot_surface(x, y, z, alpha = .3, cmap = cm.jet, rstride=1, cstride=1)
    ax.plot_surface(x[:-1], y[:-1], z[:-1], alpha = 1, cmap = cm.jet, rstride=1, cstride=1)
    ax.plot(x[-2], y[-2], z[-2], color = 'k', lw = .5)

def get_policy_3_C_S_E_B_Q_S_and_Q_B(gamma, alpha, N = 31, eta_in = 1., eta_out = 1., quantile = .99, mismatch = None):
    if N < 5:
        N = 5
    E_B = empty(N)
    C_S = empty(N)
    Q_B = empty(N)
    Q_S = empty(N)
    C_S[0], a, b, c, d, e, Q_S[0], E_B[0], Q_B[0], f = get_storage_data(gamma, alpha, 0, eta_in, eta_out, 3, quantile, mismatch)
    print E_B[0]
    C_S[N - 1], a, b, c, d, e, Q_S[N - 1], E_B[N - 1], Q_B[N - 1], f = get_storage_data(gamma, alpha, NaN, eta_in, eta_out, 3, quantile, mismatch)
    C_S[1], a, b, c, d, e, Q_S[1], E_B[1], Q_B[1], f = get_storage_data(gamma, alpha, 1e-6, eta_in, eta_out, 3, quantile, mismatch)
    print C_S[0], E_B[0], C_S[1], E_B[1], C_S[N - 1], E_B[N - 1]
    crosspoint = find_crossing_point(C_S[0], E_B[0], C_S[1], E_B[1], C_S[N - 1], E_B[N - 1], C_S[N - 1], E_B[N - 1])
    a = 1. / (1 + (N - 3) / 2)
    C_S[N - 2], a, b, c, d, e, Q_S[N - 2], E_B[N - 2], Q_B[N - 2], f = get_storage_data(gamma, alpha, a * crosspoint + (1 - a) * C_S[N - 1], eta_in, eta_out, 3, quantile, mismatch)
    print "C_S[N - 2]: ", C_S[0], C_S[1], C_S[N - 2], C_S[N - 1]
    for i in arange(N - 4):
        print i
        sys.stdout.flush()
        if i/2 != i/2.:
            print C_S[N - 2]
            crosspoint = find_crossing_point(C_S[i/2], E_B[i/2], C_S[i/2 + 1], E_B[i/2 + 1], C_S[N - 2 - i/2], E_B[N - 2 - i/2], C_S[N - 1 - i/2], E_B[N - 1 - i/2])
            a = 1. / (1 + (N - 3 - i) / 2)
            C_S[2 + i/2], a, b, c, d, e, Q_S[2 + i/2], E_B[2 + i/2], Q_B[2 + i/2], f = get_storage_data(gamma, alpha, a * crosspoint + (1 - a) * C_S[1 + i/2], eta_in, eta_out, 3, quantile, mismatch)
        else:
            crosspoint = find_crossing_point(C_S[i/2], E_B[i/2], C_S[1 + i/2], E_B[1 + i/2], C_S[N - 2 - i/2], E_B[N - 2 - i/2], C_S[N - 1 - i/2], E_B[N - 1 - i/2])
            a = 1. / (1 + (N - 3 - i) / 2)
            C_S[N - 3 - i/2], a, b, c, d, e, Q_S[N - 3 - i/2], E_B[N - 3 - i/2], Q_B[N - 3 - i/2], f = get_storage_data(gamma, alpha, a * crosspoint + (1 - a) * C_S[N - 2 - i/2], eta_in, eta_out, 3, quantile, mismatch)

    return C_S, E_B, Q_S, Q_B

def get_policy_5_C_S_E_B_Q_S_and_Q_B(gamma, alpha, N = 31, eta_in = 1., eta_out = 1., quantile = .99, mismatch = None):
    quantiles = array(quantile)
    L = quantiles.size
    if N < 5:
        N = 5
    E_B = empty(N)
    C_S = empty(N)
    Q_B = empty([N, L + 1])
#    QN = empty(N)
    Q_S = empty([N, L + 1])
    constant.policy_5_out_capacity = 1e10
    C_S[0], a, b, c, d, e, Q_S[0], E_B[0], Q_B[0], f = get_storage_data(gamma, alpha, 0., eta_in, eta_out, 2, concatenate([quantiles, [1.]]), mismatch)
    print "HEY!"
    constant.policy_5_out_capacity = 1e10
    C_S[N - 1], a, b, c, d, e, Q_S[N - 1], E_B[N - 1], Q_B[N - 1], f = get_storage_data(gamma, alpha, NaN, eta_in, eta_out, 5, concatenate([quantiles, [1. - (1./N)**2]]), mismatch)
    print "HEP!"
    for i in arange(N - 2):
        print i, N - 2
        constant.policy_5_out_capacity = max(1e-5,Q_S[N - 1 - i, L])
        print Q_S[N - 1 - i, L]
        C_S[N - 2 - i], a, b, c, d, e, Q_S[N - 2 - i], E_B[N - 2 - i], Q_B[N - 2 - i], f = get_storage_data(gamma, alpha, NaN, eta_in, eta_out, 5, concatenate([quantile, [1. - ((1. + i)/N)**2]]), mismatch)
        print "Made it!"
        sys.stdout.flush()
    return C_S, E_B, Q_S[:,:-1], Q_B[:,:-1]

def get_policy_5_data(gamma, alpha, N = 31, eta_in = 1., eta_out = 1., quantile = .99, mismatch = None):
    quantiles = array(quantile)
    L = quantiles.size
    if L == 1:
        quantiles = [quantiles]
    if N < 5:
        N = 5
    E_B = empty(N)
    E_X = empty(N)
    Q_X = empty([N, L + 1])
    E_S_in = empty(N)
    Q_S_in = empty([N, L + 1])
    E_S_out = empty(N)
    C_S = empty(N)
    Q_B = empty([N, L + 1])
#    QN = empty(N)
    Q_S_out = empty([N, L + 1])
    S_L = empty(N)
    constant.policy_5_out_capacity = 1e10
    constant.policy_5_in_capacity = 1e10
    C_S[0], E_X[0], Q_X[0], E_S_in[0], Q_S_in[0], E_S_out[0], Q_S_out[0], E_B[0], Q_B[0], S_L[0] = get_storage_data(gamma, alpha, 0., eta_in, eta_out, 2, concatenate([quantiles, [1.]]), mismatch)
    constant.policy_5_out_capacity = 1e10
    C_S[N - 1], E_X[N - 1], Q_X[N - 1], E_S_in[N - 1], Q_S_in[N - 1], E_S_in[N - 1], Q_S_out[N - 1], E_B[N - 1], Q_B[N - 1], S_L[N - 1] = get_storage_data(gamma, alpha, NaN, eta_in, eta_out, 5, concatenate([quantiles, [1. - (1./N)**2]]), mismatch)
    for i in arange(N - 2):
        constant.policy_5_out_capacity = max(1e-5,Q_S_out[N - 1 - i, L])
        C_S[N - 2 - i], E_X[N - 2 - i], Q_X[N - 2 -i], E_S_in[N - 2 - i], Q_S_in[N - 2 - i], E_S_out[N - 2 - i], Q_S_out[N - 2 - i], E_B[N - 2 - i], Q_B[N - 2 - i], S_L[N - 2 - i] = get_storage_data(gamma, alpha, NaN, eta_in, eta_out, 5, concatenate([quantiles, [1. - ((1. + i)/N)**2]]), mismatch)
    return C_S, E_X, Q_X[:,:-1], E_S_in, Q_S_in[:,:-1], E_S_out, Q_S_out[:,:-1], E_B, Q_B[:,:-1], S_L

def get_policy_55_data(gamma, alpha, N = 31, eta_in = 1., eta_out = 1., quantile = .99, mismatch = None):
    quantiles = array(quantile)
    L = quantiles.size
    if L == 1:
        quantiles = [quantiles]
    if N < 5:
        N = 5
    E_B = empty(N)
    E_X = empty(N)
    Q_X = empty([N, L + 1])
    E_S_in = empty(N)
    Q_S_in = empty([N, L + 1])
    E_S_out = empty(N)
    C_S = empty(N)
    Q_B = empty([N, L + 1])
#    QN = empty(N)
    Q_S_out = empty([N, L + 1])
    S_L = empty(N)
    constant.policy_5_out_capacity = 1e10
    constant.policy_5_in_capacity = constant.policy_5_ratio * constant.policy_5_out_capacity
    C_S[0], E_X[0], Q_X[0], E_S_in[0], Q_S_in[0], E_S_out[0], Q_S_out[0], E_B[0], Q_B[0], S_L[0] = get_storage_data(gamma, alpha, 0., eta_in, eta_out, 2, concatenate([quantiles, [1.]]), mismatch)
    constant.policy_5_out_capacity = 1e10
    constant.policy_5_in_capacity = constant.policy_5_ratio * constant.policy_5_out_capacity
    C_S[N - 1], E_X[N - 1], Q_X[N - 1], E_S_in[N - 1], Q_S_in[N - 1], E_S_in[N - 1], Q_S_out[N - 1], E_B[N - 1], Q_B[N - 1], S_L[N - 1] = get_storage_data(gamma, alpha, NaN, eta_in, eta_out, 5, concatenate([quantiles, [1. - (1./N)**2]]), mismatch)
    for i in arange(N - 2):
        constant.policy_5_out_capacity = max(1e-5,Q_S_out[N - 1 - i, L])
        constant.policy_5_in_capacity = constant.policy_5_ratio * constant.policy_5_out_capacity
        print constant.policy_5_in_capacity, constant.policy_5_out_capacity
        sys.stdout.flush()
        C_S[N - 2 - i], E_X[N - 2 - i], Q_X[N - 2 -i], E_S_in[N - 2 - i], Q_S_in[N - 2 - i], E_S_out[N - 2 - i], Q_S_out[N - 2 - i], E_B[N - 2 - i], Q_B[N - 2 - i], S_L[N - 2 - i] = get_storage_data(gamma, alpha, NaN, eta_in, eta_out, 5, concatenate([quantiles, [1. - ((1. + i)/N)**2]]), mismatch)
    return C_S, E_X, Q_X[:,:-1], E_S_in, Q_S_in[:,:-1], E_S_out, Q_S_out[:,:-1], E_B, Q_B[:,:-1], S_L


def get_policy_7_data(mismatch, N, eta_in, eta_out, quantile = .99):
    K_B = linspace(0,(-min(mismatch))**.2, N)**5
    quantiles = array(quantile)
    L = quantiles.size
    E_B = empty(N)
    E_X = empty(N)
    Q_X = empty([N, L])
    E_S_in = empty(N)
    Q_S_in = empty([N, L])
    E_S_out = empty(N)
    C_S = empty(N)
    Q_B = empty([N, L])
    Q_S_out = empty([N, L])
    S_L = empty(N)
    for i in arange(N):
        constant.policy_7_capacity = K_B[i]
        C_S[i], E_X[i], Q_X[i], E_S_in[i], Q_S_in[i], E_S_out[i], Q_S_out[i], E_B[i], Q_B[i], S_L[i] = get_storage_data(NaN, 1, NaN, eta_in, eta_out, 7, quantile, mismatch)
    return C_S, E_X, Q_X, E_S_in, Q_S_in, E_S_out, Q_S_out, E_B, Q_B, S_L

def get_policy_8_data(mismatch, N, eta_in, eta_out, quantile = .99):
    K_B = linspace(0,(-min(mismatch))**.25, N)**4
    quantiles = array(quantile)
    L = quantiles.size
    E_B = empty(N)
    E_X = empty(N)
    Q_X = empty([N, L])
    E_S_in = empty(N)
    Q_S_in = empty([N, L])
    E_S_out = empty(N)
    C_S = empty(N)
    Q_B = empty([N, L])
    Q_S_out = empty([N, L])
    S_L = empty(N)
    for i in arange(N):
        constant.policy_8_capacity = K_B[i]
        print K_B[i]
        C_S[i], E_X[i], Q_X[i], E_S_in[i], Q_S_in[i], E_S_out[i], Q_S_out[i], E_B[i], Q_B[i], S_L[i] = get_storage_data(NaN, 1, NaN, eta_in, eta_out, 8, quantile, mismatch)
    return C_S, E_X, Q_X, E_S_in, Q_S_in, E_S_out, Q_S_out, E_B, Q_B, S_L

def get_policy_6_data(gamma, alpha, N = 31, eta_in = 1., eta_out = 1., quantile = .99, mismatch = None):
    quantiles = array(quantile)
    L = quantiles.size
    if N < 5:
        N = 5
    E_B = empty(N)
    E_X = empty(N)
    Q_X = empty([N, L])
    E_S_in = empty(N)
    Q_S_in = empty([N, L])
    E_S_out = empty(N)
    C_S = empty(N)
    Q_B = empty([N, L])
    QN = max(mismatch) * (exp(linspace(-N/10.,N/10.,N))**.25/exp(N/10.)**.25)[1:]
    Q_S_out = empty([N, L])
    S_L = empty(N)
    constant.policy_5_in_capacity = 1e10
    constant_policy_5_out_capacity = 1e10
    C_S[0], E_X[0], Q_X[0], E_S_in[0], Q_S_in[0], E_S_out[0], Q_S_out[0], E_B[0], Q_B[0], S_L[0] = get_storage_data(gamma, alpha, 0., eta_in, eta_out, 2, quantiles, mismatch)
    for i in arange(N - 1):
        constant.policy_5_in_capacity = QN[i]
        print constant.policy_5_in_capacity, constant.policy_5_out_capacity
        sys.stdout.flush()
        C_S[1 + i], E_X[1 + i], Q_X[1 + i], E_S_in[1 + i], Q_S_in[1 + i], E_S_out[1 + i], Q_S_out[1 + i], E_B[1 + i], Q_B[1 + i], S_L[1 + i] = get_storage_data(gamma, alpha, NaN, eta_in, eta_out, 5, quantiles, mismatch)
    return C_S, E_X, Q_X, E_S_in, Q_S_in, E_S_out, Q_S_out, E_B, Q_B, S_L

def plot_20110505_daily_figure(alpha):
    constant.policy_5_capacity = 1e10
    sc1, ex1, qx1, eb1, qb1 = get_daily_no_balancing_data(alpha, .6, .6, 301, policy = 5)
    constant.policy_5_capacity = 24.
    sc2, ex2, qx2, eb2, qb2 = get_daily_no_balancing_data(alpha, .6, .6, 301, policy = 5)
    constant.policy_5_capacity = 18.
    sc3, ex3, qx3, eb3, qb3 = get_daily_no_balancing_data(alpha, .6, .6, 301, policy = 5)
    constant.policy_5_capacity = 12.
    sc4, ex4, qx4, eb4, qb4 = get_daily_no_balancing_data(alpha, .6, .6, 301, policy = 5)
    constant.policy_5_capacity = 6.
    sc5, ex5, qx5, eb5, qb5 = get_daily_no_balancing_data(alpha, .6, .6, 301, policy = 5)
    constant.policy_5_capacity = 3.
    sc6, ex6, qx6, eb6, qb6 = get_daily_no_balancing_data(alpha, .6, .6, 301, policy = 5)
    constant.policy_5_capacity = 2.
    sc7, ex7, qx7, eb7, qb7 = get_daily_no_balancing_data(alpha, .6, .6, 301, policy = 5)
    constant.policy_5_capacity = 1.
    sc8, ex8, qx8, eb8, qb8 = get_daily_no_balancing_data(alpha, .6, .6, 301, policy = 5)
    constant.policy_5_capacity = .5
    sc9, ex9, qx9, eb9, qb9 = get_daily_no_balancing_data(alpha, .6, .6, 301, policy = 5)
    figure()
    subplot(211)
    plot(linspace(.5,1.5,301), sc2/8766, 'g-', lw = 2, label = r"$K_S=1.000$")
    plot(linspace(.5,1.5,301), sc3/8766, 'g--', lw = 2, label = r"$K_S=0.750$")
    plot(linspace(.5,1.5,301), sc4/8766, 'g-.', lw = 2, label = r"$K_S=0.500$")
    plot(linspace(.5,1.5,301), sc5/8766, 'g:', lw = 2, label = r"$K_S=0.250$")
    plot(linspace(.5,1.5,301), sc6/8766, 'g-', lw = 1, label = r"$K_S=0.125$")
    plot(linspace(.5,1.5,301), sc7/8766, 'g--', lw = 1, label = r"$K_S=0.083$")
    plot(linspace(.5,1.5,301), sc8/8766, 'g-.', lw = 1, label = r"$K_S=0.042$")
    plot(linspace(.5,1.5,301), sc9/8766, 'g:', lw = 1, label = r"$K_S=0.021$")
    plot(linspace(.5,1.5,301), sc1/8766, 'k:', lw = 1, label = r"$K_S=\infty$")
    leg1 = legend(loc = 2, labelspacing = .1)
    setp(leg1.get_texts(), fontsize = "small")
    xticks(linspace(.5,1.5,11))
    axis("tight")
    ylabel("Storage energy capacity\n[average yearly consumption]")
    subplot(212)
    plot(linspace(.5,1.5,301), eb2/8766, 'r-', lw = 2, label = r"balancing")
    plot(linspace(.5,1.5,301), eb3/8766, 'r--', lw = 2)
    plot(linspace(.5,1.5,301), eb4/8766, 'r-.', lw = 2)
    plot(linspace(.5,1.5,301), eb5/8766, 'r:', lw = 2)
    plot(linspace(.5,1.5,301), eb6/8766, 'r-', lw = 1)
    plot(linspace(.5,1.5,301), eb7/8766, 'r--', lw = 1)
    plot(linspace(.5,1.5,301), eb8/8766, 'r-.', lw = 1)
    plot(linspace(.5,1.5,301), eb9/8766, 'r:', lw = 1)
    plot(linspace(.5,1.5,301), ex2/8766, 'c-', lw = 2, label = r"excess")
    plot(linspace(.5,1.5,301), ex3/8766, 'c--', lw = 2)
    plot(linspace(.5,1.5,301), ex4/8766, 'c-.', lw = 2)
    plot(linspace(.5,1.5,301), ex5/8766, 'c:', lw = 2)
    plot(linspace(.5,1.5,301), ex6/8766, 'c-', lw = 1)
    plot(linspace(.5,1.5,301), ex7/8766, 'c--', lw = 1)
    plot(linspace(.5,1.5,301), ex8/8766, 'c-.', lw = 1)
    plot(linspace(.5,1.5,301), ex9/8766, 'c:', lw = 1)
    plot(linspace(.5,1.5,301), eb1/8766, 'k:', lw = 1)
    plot(linspace(.5,1.5,301), ex1/8766, 'k:', lw = 1)
    leg2 = legend(loc = "upper center", labelspacing = .1)
    setp(leg2.get_texts(), fontsize = "small")
    xlabel("Average renewable generation in average load hours")
    ylabel("[average load hours]")
    axis("tight")
    xticks(linspace(.5,1.5,11))
    subplot(211)

def plot_all_policies(mismatch, N, eta_in, eta_out, quantile = .99):
    #policy 2:
    p2sc, p2me, p2qe, p2msi, p2qsi, p2mso, p2qso, p2mb, p2qb, p2msl = get_storages_data(NaN, 1, N, eta_in, eta_out, 2, quantile, concatenate([[0],exp(linspace(-31/10.,(N-1)/10.,(N-1)))/exp((N-1)/10.)]), mismatch)
    #policy 3:
    p3sc, p3me, p3qe, p3msi, p3qsi, p3mso, p3qso, p3mb, p3qb, p3msl = get_storages_data(NaN, 1, N, eta_in, eta_out, 3, quantile, concatenate([[0],exp(linspace(-31/10.,(N-1)/10.,(N-1)))/exp((N-1)/10.)]), mismatch)
    #policy 4:
    p4sc, p4me, p4qe, p4msi, p4qsi, p4mso, p4qso, p4mb, p4qb, p4msl = get_storages_data(NaN, 1, N, eta_in, eta_out, 4, quantile, concatenate([[0],exp(linspace(-31/10.,(N-1)/10.,(N-1)))/exp((N-1)/10.)]), mismatch)
    #policy 5:
    constant.policy_5_in_capacity = 1e10
    p5sc, p5me, p5qe, p5msi, p5qsi, p5mso, p5qso, p5mb, p5qb, p5msl = get_policy_5_data(NaN, 1, N, eta_in, eta_out, quantile, mismatch)
    #policy 5.5:
    p55sc, p55me, p55qe, p55msi, p55qsi, p55mso, p55qso, p55mb, p55qb, p55msl = get_policy_55_data(NaN, 1, N, eta_in, eta_out, quantile, mismatch)
    #policy 6:
    constant.policy_5_out_capacity = 1e10
    p6sc, p6me, p6qe, p6msi, p6qsi, p6mso, p6qso, p6mb, p6qb, p6msl = get_policy_6_data(NaN, 1, N, eta_in, eta_out, quantile, mismatch)
    #policy 7:
    p7sc, p7me, p7qe, p7msi, p7qsi, p7mso, p7qso, p7mb, p7qb, p7msl = get_policy_7_data(mismatch, N, eta_in, eta_out, quantile)
    #policy 8:
    p8sc, p8me, p8qe, p8msi, p8qsi, p8mso, p8qso, p8mb, p8qb, p8msl = get_policy_8_data(mismatch, N, eta_in, eta_out, quantile)
#    clf()
    print len(p2sc), len(p2mb)
    plot(p8sc,p8mb, lw = 8)
    plot(p7sc,p7mb, lw = 7)
    plot(p6sc,p6mb, lw = 6)
    plot(p5sc,p5mb, lw = 5)
    plot(p55sc, p55mb, lw = 5.5)
    plot(p4sc,p4mb, lw = 4)
    plot(p3sc,p3mb, lw = 3)
    plot(p2sc,p2mb, lw = 2)
    plot(p2sc,p2mb, 'x')
    plot(p3sc,p3mb, 'x')
    plot(p4sc,p4mb, 'x')
    plot(p5sc,p5mb, 'x')
    plot(p55sc, p55mb, 'x')
    plot(p6sc,p6mb, 'x')
    plot(p7sc,p7mb, 'x')
    plot(p8sc,p8mb, 'x')
    p2 = (p2sc, p2me, p2qe, p2msi, p2qsi, p2mso, p2qso, p2mb, p2qb, p2msl)
    p3 = (p3sc, p3me, p3qe, p3msi, p3qsi, p3mso, p3qso, p3mb, p3qb, p3msl)
    p4 = (p4sc, p4me, p4qe, p4msi, p4qsi, p4mso, p4qso, p4mb, p4qb, p4msl)
    p5 = (p5sc, p5me, p5qe, p5msi, p5qsi, p5mso, p5qso, p5mb, p5qb, p5msl)
    p55 = (p55sc, p55me, p55qe, p55msi, p55qsi, p55mso, p55qso, p55mb, p55qb, p55msl)
    p6 = (p6sc, p6me, p6qe, p6msi, p6qsi, p6mso, p6qso, p6mb, p6qb, p6msl)
    p7 = (p7sc, p7me, p7qe, p7msi, p7qsi, p7mso, p7qso, p7mb, p7qb, p7msl)
    p8 = (p8sc, p8me, p8qe, p8msi, p8qsi, p8mso, p8qso, p8mb, p8qb, p8msl)
    return p2, p3, p4, p5, p55, p6, p7, p8

def plot_E_B_vs_C_S(p2,p3,p4,p5,p55,p6,p7,p8):
    p2sc, p2me, p2qe, p2msi, p2qsi, p2mso, p2qso, p2mb, p2qb, p2msl = p2
    p3sc, p3me, p3qe, p3msi, p3qsi, p3mso, p3qso, p3mb, p3qb, p3msl = p3
    p4sc, p4me, p4qe, p4msi, p4qsi, p4mso, p4qso, p4mb, p4qb, p4msl = p4
    p5sc, p5me, p5qe, p5msi, p5qsi, p5mso, p5qso, p5mb, p5qb, p5msl = p5
    p55sc, p55me, p55qe, p55msi, p55qsi, p55mso, p55qso, p55mb, p55qb, p55msl = p55
    p6sc, p6me, p6qe, p6msi, p6qsi, p6mso, p6qso, p6mb, p6qb, p6msl = p6
    p7sc, p7me, p7qe, p7msi, p7qsi, p7mso, p7qso, p7mb, p7qb, p7msl = p7
    p8sc, p8me, p8qe, p8msi, p8qsi, p8mso, p8qso, p8mb, p8qb, p8msl = p8

    plot(p2sc, p2mb, lw = 3, label = "Policy 2", color = '#FF1493')
    plot(p3sc, p3mb, '--', lw = 3, label = "Policy 3", color = "#FF0000") 
    plot(p4sc, p4mb, ':', lw = 3, label = "Policy 4", color = "#FFA500") 
    plot(p5sc, p5mb, lw = 3, label = "Policy 5", color = "#FFFF00")
    plot(p55sc, p55mb, lw = 3, label = "Policy 5.5", color = "#008000")
    plot(p6sc, p6mb, lw = 3, label = "Policy 6", color = "#0000FF")
    plot(p7sc, p7mb, lw = 3, label = "Policy 7", color = "#000000") 
    plot(p8sc, p8mb, lw = 3, label = "Policy 8", color = "#800080") 

def plot_E_B_vs_K_S_out(p2,p3,p4,p5,p55,p6,p7,p8):
    p2sc, p2me, p2qe, p2msi, p2qsi, p2mso, p2qso, p2mb, p2qb, p2msl = p2
    p3sc, p3me, p3qe, p3msi, p3qsi, p3mso, p3qso, p3mb, p3qb, p3msl = p3
    p4sc, p4me, p4qe, p4msi, p4qsi, p4mso, p4qso, p4mb, p4qb, p4msl = p4
    p5sc, p5me, p5qe, p5msi, p5qsi, p5mso, p5qso, p5mb, p5qb, p5msl = p5
    p55sc, p55me, p55qe, p55msi, p55qsi, p55mso, p55qso, p55mb, p55qb, p55msl = p55
    p6sc, p6me, p6qe, p6msi, p6qsi, p6mso, p6qso, p6mb, p6qb, p6msl = p6
    p7sc, p7me, p7qe, p7msi, p7qsi, p7mso, p7qso, p7mb, p7qb, p7msl = p7
    p8sc, p8me, p8qe, p8msi, p8qsi, p8mso, p8qso, p8mb, p8qb, p8msl = p8

    plot(p2qso[:,-1], p2mb, lw = 3, label = "Policy 2", color = '#FF1493')
    plot(p3qso[:,-1], p3mb, lw = 3, label = "Policy 3", color = "#FF0000") 
    plot(p4qso[:,-1], p4mb, lw = 3, label = "Policy 4", color = "#FFA500") 
    plot(p5qso[:,-1], p5mb, lw = 3, label = "Policy 5", color = "#FFFF00")
    plot(p55qso[:,-1], p55mb, lw = 3, label = "Policy 5.5", color = "#008000")
    plot(p6qso[:,-1], p6mb, lw = 3, label = "Policy 6", color = "#0000FF")
    plot(p7qso[:,-1], p7mb, lw = 3, label = "Policy 7", color = "#000000") 
    plot(p8qso[:,-1], p8mb, lw = 3, label = "Policy 8", color = "#800080") 

def plot_E_B_vs_K_S_in(p2,p3,p4,p5,p55,p6,p7,p8):
    p2sc, p2me, p2qe, p2msi, p2qsi, p2mso, p2qso, p2mb, p2qb, p2msl = p2
    p3sc, p3me, p3qe, p3msi, p3qsi, p3mso, p3qso, p3mb, p3qb, p3msl = p3
    p4sc, p4me, p4qe, p4msi, p4qsi, p4mso, p4qso, p4mb, p4qb, p4msl = p4
    p5sc, p5me, p5qe, p5msi, p5qsi, p5mso, p5qso, p5mb, p5qb, p5msl = p5
    p55sc, p55me, p55qe, p55msi, p55qsi, p55mso, p55qso, p55mb, p55qb, p55msl = p55
    p6sc, p6me, p6qe, p6msi, p6qsi, p6mso, p6qso, p6mb, p6qb, p6msl = p6
    p7sc, p7me, p7qe, p7msi, p7qsi, p7mso, p7qso, p7mb, p7qb, p7msl = p7
    p8sc, p8me, p8qe, p8msi, p8qsi, p8mso, p8qso, p8mb, p8qb, p8msl = p8

    plot(p2qsi[:,-1], p2mb, lw = 3, label = "Policy 2", color = '#FF1493')
    plot(p3qsi[:,-1], p3mb, lw = 3, label = "Policy 3", color = "#FF0000") 
    plot(p4qsi[:,-1], p4mb, lw = 3, label = "Policy 4", color = "#FFA500") 
    plot(p5qsi[:,-1], p5mb, lw = 3, label = "Policy 5", color = "#FFFF00")
    plot(p55qsi[:,-1], p55mb, lw = 3, label = "Policy 5.5", color = "#008000")
    plot(p6qsi[:,-1], p6mb, lw = 3, label = "Policy 6", color = "#0000FF")
    plot(p7qsi[:,-1], p7mb, lw = 3, label = "Policy 7", color = "#000000") 
    plot(p8qsi[:,-1], p8mb, lw = 3, label = "Policy 8", color = "#800080") 

def plot_Q_B_vs_K_S_out(p2,p3,p4,p5,p55,p6,p7,p8):
    p2sc, p2me, p2qe, p2msi, p2qsi, p2mso, p2qso, p2mb, p2qb, p2msl = p2
    p3sc, p3me, p3qe, p3msi, p3qsi, p3mso, p3qso, p3mb, p3qb, p3msl = p3
    p4sc, p4me, p4qe, p4msi, p4qsi, p4mso, p4qso, p4mb, p4qb, p4msl = p4
    p5sc, p5me, p5qe, p5msi, p5qsi, p5mso, p5qso, p5mb, p5qb, p5msl = p5
    p55sc, p55me, p55qe, p55msi, p55qsi, p55mso, p55qso, p55mb, p55qb, p55msl = p55
    p6sc, p6me, p6qe, p6msi, p6qsi, p6mso, p6qso, p6mb, p6qb, p6msl = p6
    p7sc, p7me, p7qe, p7msi, p7qsi, p7mso, p7qso, p7mb, p7qb, p7msl = p7
    p8sc, p8me, p8qe, p8msi, p8qsi, p8mso, p8qso, p8mb, p8qb, p8msl = p8

    plot(p2qso[:,-1], p2qb[:,0], lw = 3, label = "Policy 2", color = '#FF1493')
    plot(p3qso[:,-1], p3qb[:,0], lw = 3, label = "Policy 3", color = "#FF0000") 
    plot(p4qso[:,-1], p4qb[:,0], lw = 3, label = "Policy 4", color = "#FFA500") 
    plot(p5qso[:,-1], p5qb[:,0], lw = 3, label = "Policy 5", color = "#FFFF00")
    plot(p55qso[:,-1], p55qb[:,0], lw = 3, label = "Policy 5.5", color = "#008000")
    plot(p6qso[:,-1], p6qb[:,0], lw = 3, label = "Policy 6", color = "#0000FF")
    plot(p7qso[:,-1], p7qb[:,0], lw = 3, label = "Policy 7", color = "#000000") 
    plot(p8qso[:,-1], p8qb[:,0], lw = 3, label = "Policy 8", color = "#800080") 

    plot(p2qso[:,-1], p2qb[:,1], '--', lw = 3, color = '#FF1493')
    plot(p3qso[:,-1], p3qb[:,1], '--', lw = 3, color = "#FF0000") 
    plot(p4qso[:,-1], p4qb[:,1], '--', lw = 3, color = "#FFA500") 
    plot(p5qso[:,-1], p5qb[:,1], '--', lw = 3, color = "#FFFF00")
    plot(p55qso[:,-1], p55qb[:,1], '--', lw = 3, color = "#008000")
    plot(p6qso[:,-1], p6qb[:,1], '--', lw = 3, color = "#0000FF")
    plot(p7qso[:,-1], p7qb[:,1], '--', lw = 3, color = "#000000") 
    plot(p8qso[:,-1], p8qb[:,1], '--', lw = 3, color = "#800080") 

    plot(p2qso[:,-1], p2qb[:,2], '-.', lw = 1, color = '#FF1493')
    plot(p3qso[:,-1], p3qb[:,2], '-.', lw = 1, color = "#FF0000") 
    plot(p4qso[:,-1], p4qb[:,2], '-.', lw = 1, color = "#FFA500") 
    plot(p5qso[:,-1], p5qb[:,2], '-.', lw = 1, color = "#FFFF00")
    plot(p55qso[:,-1], p55qb[:,2], '-.', lw = 1, color = "#008000")
    plot(p6qso[:,-1], p6qb[:,2], '-.', lw = 1, color = "#0000FF")
    plot(p7qso[:,-1], p7qb[:,2], '-.', lw = 1, color = "#000000") 
    plot(p8qso[:,-1], p8qb[:,2], '-.', lw = 1, color = "#800080") 

    plot(p2qso[:,-1], p2qb[:,3], ':', lw = 1, color = '#FF1493')
    plot(p3qso[:,-1], p3qb[:,3], ':', lw = 1, color = "#FF0000") 
    plot(p4qso[:,-1], p4qb[:,3], ':', lw = 1, color = "#FFA500") 
    plot(p5qso[:,-1], p5qb[:,3], ':', lw = 1, color = "#FFFF00")
    plot(p55qso[:,-1], p55qb[:,3], ':', lw = 1, color = "#008000")
    plot(p6qso[:,-1], p6qb[:,3], ':', lw = 1, color = "#0000FF")
    plot(p7qso[:,-1], p7qb[:,3], ':', lw = 1, color = "#000000") 
    plot(p8qso[:,-1], p8qb[:,3], ':', lw = 1, color = "#800080") 

#    plot(p2qso[:,-1], p2qb[:,4], '.', lw = 1, color = '#FF1493')
#    plot(p3qso[:,-1], p3qb[:,4], '.', lw = 1, color = "#FF0000") 
#    plot(p4qso[:,-1], p4qb[:,4], '.', lw = 1, color = "#FFA500") 
#    plot(p5qso[:,-1], p5qb[:,4], '.', lw = 1, color = "#FFFF00")
#    plot(p55qso[:,-1], p55qb[:,4], '.', lw = 1, color = "#008000")
#    plot(p6qso[:,-1], p6qb[:,4], '.', lw = 1, color = "#0000FF")
#    plot(p7qso[:,-1], p7qb[:,4], '.', lw = 1, color = "#000000") 
#    plot(p8qso[:,-1], p8qb[:,4], '.', lw = 1, color = "#800080") 

def plot_Q_B_vs_C_S(p2,p3,p4,p5,p55,p6,p7,p8):
    p2sc, p2me, p2qe, p2msi, p2qsi, p2mso, p2qso, p2mb, p2qb, p2msl = p2
    p3sc, p3me, p3qe, p3msi, p3qsi, p3mso, p3qso, p3mb, p3qb, p3msl = p3
    p4sc, p4me, p4qe, p4msi, p4qsi, p4mso, p4qso, p4mb, p4qb, p4msl = p4
    p5sc, p5me, p5qe, p5msi, p5qsi, p5mso, p5qso, p5mb, p5qb, p5msl = p5
    p55sc, p55me, p55qe, p55msi, p55qsi, p55mso, p55qso, p55mb, p55qb, p55msl = p55
    p6sc, p6me, p6qe, p6msi, p6qsi, p6mso, p6qso, p6mb, p6qb, p6msl = p6
    p7sc, p7me, p7qe, p7msi, p7qsi, p7mso, p7qso, p7mb, p7qb, p7msl = p7
    p8sc, p8me, p8qe, p8msi, p8qsi, p8mso, p8qso, p8mb, p8qb, p8msl = p8

    plot(p2sc, p2qb[:,0], lw = 3, label = "Policy 2", color = '#FF1493')
    plot(p3sc, p3qb[:,0], lw = 3, label = "Policy 3", color = "#FF0000") 
    plot(p4sc, p4qb[:,0], lw = 3, label = "Policy 4", color = "#FFA500") 
    plot(p5sc, p5qb[:,0], lw = 3, label = "Policy 5", color = "#FFFF00")
    plot(p55sc, p55qb[:,0], lw = 3, label = "Policy 5.5", color = "#008000")
    plot(p6sc, p6qb[:,0], lw = 3, label = "Policy 6", color = "#0000FF")
    plot(p7sc, p7qb[:,0], lw = 3, label = "Policy 7", color = "#000000") 
    plot(p8sc, p8qb[:,0], lw = 3, label = "Policy 8", color = "#800080") 

    plot(p2sc, p2qb[:,1], '--', lw = 3, color = '#FF1493')
    plot(p3sc, p3qb[:,1], '--', lw = 3, color = "#FF0000") 
    plot(p4sc, p4qb[:,1], '--', lw = 3, color = "#FFA500") 
    plot(p5sc, p5qb[:,1], '--', lw = 3, color = "#FFFF00")
    plot(p55sc, p55qb[:,1], '--', lw = 3, color = "#008000")
    plot(p6sc, p6qb[:,1], '--', lw = 3, color = "#0000FF")
    plot(p7sc, p7qb[:,1], '--', lw = 3, color = "#000000") 
    plot(p8sc, p8qb[:,1], '--', lw = 3, color = "#800080") 

    plot(p2sc, p2qb[:,2], '-.', lw = 1, color = '#FF1493')
    plot(p3sc, p3qb[:,2], '-.', lw = 1, color = "#FF0000") 
    plot(p4sc, p4qb[:,2], '-.', lw = 1, color = "#FFA500") 
    plot(p5sc, p5qb[:,2], '-.', lw = 1, color = "#FFFF00")
    plot(p55sc, p55qb[:,2], '-.', lw = 1, color = "#008000")
    plot(p6sc, p6qb[:,2], '-.', lw = 1, color = "#0000FF")
    plot(p7sc, p7qb[:,2], '-.', lw = 1, color = "#000000") 
    plot(p8sc, p8qb[:,2], '-.', lw = 1, color = "#800080") 

    plot(p2sc, p2qb[:,3], ':', lw = 1, color = '#FF1493')
    plot(p3sc, p3qb[:,3], ':', lw = 1, color = "#FF0000") 
    plot(p4sc, p4qb[:,3], ':', lw = 1, color = "#FFA500") 
    plot(p5sc, p5qb[:,3], ':', lw = 1, color = "#FFFF00")
    plot(p55sc, p55qb[:,3], ':', lw = 1, color = "#008000")
    plot(p6sc, p6qb[:,3], ':', lw = 1, color = "#0000FF")
    plot(p7sc, p7qb[:,3], ':', lw = 1, color = "#000000") 
    plot(p8sc, p8qb[:,3], ':', lw = 1, color = "#800080") 

#    plot(p2qso[:,-1], p2qb[:,4], '.', lw = 1, color = '#FF1493')
#    plot(p3qso[:,-1], p3qb[:,4], '.', lw = 1, color = "#FF0000") 
#    plot(p4qso[:,-1], p4qb[:,4], '.', lw = 1, color = "#FFA500") 
#    plot(p5qso[:,-1], p5qb[:,4], '.', lw = 1, color = "#FFFF00")
#    plot(p55qso[:,-1], p55qb[:,4], '.', lw = 1, color = "#008000")
#    plot(p6qso[:,-1], p6qb[:,4], '.', lw = 1, color = "#0000FF")
#    plot(p7qso[:,-1], p7qb[:,4], '.', lw = 1, color = "#000000") 
#    plot(p8qso[:,-1], p8qb[:,4], '.', lw = 1, color = "#800080") 

def plot_Q_S_out_vs_C_S(p2,p3,p4,p5,p55,p6,p7,p8):
    p2sc, p2me, p2qe, p2msi, p2qsi, p2mso, p2qso, p2mb, p2qb, p2msl = p2
    p3sc, p3me, p3qe, p3msi, p3qsi, p3mso, p3qso, p3mb, p3qb, p3msl = p3
    p4sc, p4me, p4qe, p4msi, p4qsi, p4mso, p4qso, p4mb, p4qb, p4msl = p4
    p5sc, p5me, p5qe, p5msi, p5qsi, p5mso, p5qso, p5mb, p5qb, p5msl = p5
    p55sc, p55me, p55qe, p55msi, p55qsi, p55mso, p55qso, p55mb, p55qb, p55msl = p55
    p6sc, p6me, p6qe, p6msi, p6qsi, p6mso, p6qso, p6mb, p6qb, p6msl = p6
    p7sc, p7me, p7qe, p7msi, p7qsi, p7mso, p7qso, p7mb, p7qb, p7msl = p7
    p8sc, p8me, p8qe, p8msi, p8qsi, p8mso, p8qso, p8mb, p8qb, p8msl = p8

    plot(p2sc, p2qso[:,0], lw = 3, label = "Policy 2", color = '#FF1493')
    plot(p3sc, p3qso[:,0], lw = 3, label = "Policy 3", color = "#FF0000") 
    plot(p4sc, p4qso[:,0], lw = 3, label = "Policy 4", color = "#FFA500") 
    plot(p5sc, p5qso[:,0], lw = 3, label = "Policy 5", color = "#FFFF00")
    plot(p55sc, p55qso[:,0], lw = 3, label = "Policy 5.5", color = "#008000")
    plot(p6sc, p6qso[:,0], lw = 3, label = "Policy 6", color = "#0000FF")
    plot(p7sc, p7qso[:,0], lw = 3, label = "Policy 7", color = "#000000") 
    plot(p8sc, p8qso[:,0], lw = 3, label = "Policy 8", color = "#800080") 

    plot(p2sc, p2qso[:,1], '--', lw = 3, color = '#FF1493')
    plot(p3sc, p3qso[:,1], '--', lw = 3, color = "#FF0000") 
    plot(p4sc, p4qso[:,1], '--', lw = 3, color = "#FFA500") 
    plot(p5sc, p5qso[:,1], '--', lw = 3, color = "#FFFF00")
    plot(p55sc, p55qso[:,1], '--', lw = 3, color = "#008000")
    plot(p6sc, p6qso[:,1], '--', lw = 3, color = "#0000FF")
    plot(p7sc, p7qso[:,1], '--', lw = 3, color = "#000000") 
    plot(p8sc, p8qso[:,1], '--', lw = 3, color = "#800080") 

    plot(p2sc, p2qso[:,2], '-.', lw = 1, color = '#FF1493')
    plot(p3sc, p3qso[:,2], '-.', lw = 1, color = "#FF0000") 
    plot(p4sc, p4qso[:,2], '-.', lw = 1, color = "#FFA500") 
    plot(p5sc, p5qso[:,2], '-.', lw = 1, color = "#FFFF00")
    plot(p55sc, p55qso[:,2], '-.', lw = 1, color = "#008000")
    plot(p6sc, p6qso[:,2], '-.', lw = 1, color = "#0000FF")
    plot(p7sc, p7qso[:,2], '-.', lw = 1, color = "#000000") 
    plot(p8sc, p8qso[:,2], '-.', lw = 1, color = "#800080") 

    plot(p2sc, p2qso[:,3], ':', lw = 1, color = '#FF1493')
    plot(p3sc, p3qso[:,3], ':', lw = 1, color = "#FF0000") 
    plot(p4sc, p4qso[:,3], ':', lw = 1, color = "#FFA500") 
    plot(p5sc, p5qso[:,3], ':', lw = 1, color = "#FFFF00")
    plot(p55sc, p55qso[:,3], ':', lw = 1, color = "#008000")
    plot(p6sc, p6qso[:,3], ':', lw = 1, color = "#0000FF")
    plot(p7sc, p7qso[:,3], ':', lw = 1, color = "#000000") 
    plot(p8sc, p8qso[:,3], ':', lw = 1, color = "#800080") 

#    plot(p2qso[:,-1], p2qb[:,4], '.', lw = 1, color = '#FF1493')
#    plot(p3qso[:,-1], p3qb[:,4], '.', lw = 1, color = "#FF0000") 
#    plot(p4qso[:,-1], p4qb[:,4], '.', lw = 1, color = "#FFA500") 
#    plot(p5qso[:,-1], p5qb[:,4], '.', lw = 1, color = "#FFFF00")
#    plot(p55qso[:,-1], p55qb[:,4], '.', lw = 1, color = "#008000")
#    plot(p6qso[:,-1], p6qb[:,4], '.', lw = 1, color = "#0000FF")
#    plot(p7qso[:,-1], p7qb[:,4], '.', lw = 1, color = "#000000") 
#    plot(p8qso[:,-1], p8qb[:,4], '.', lw = 1, color = "#800080") 

def plot_Q_S_in_vs_C_S(p2,p3,p4,p5,p55,p6,p7,p8):
    p2sc, p2me, p2qe, p2msi, p2qsi, p2mso, p2qso, p2mb, p2qb, p2msl = p2
    p3sc, p3me, p3qe, p3msi, p3qsi, p3mso, p3qso, p3mb, p3qb, p3msl = p3
    p4sc, p4me, p4qe, p4msi, p4qsi, p4mso, p4qso, p4mb, p4qb, p4msl = p4
    p5sc, p5me, p5qe, p5msi, p5qsi, p5mso, p5qso, p5mb, p5qb, p5msl = p5
    p55sc, p55me, p55qe, p55msi, p55qsi, p55mso, p55qso, p55mb, p55qb, p55msl = p55
    p6sc, p6me, p6qe, p6msi, p6qsi, p6mso, p6qso, p6mb, p6qb, p6msl = p6
    p7sc, p7me, p7qe, p7msi, p7qsi, p7mso, p7qso, p7mb, p7qb, p7msl = p7
    p8sc, p8me, p8qe, p8msi, p8qsi, p8mso, p8qso, p8mb, p8qb, p8msl = p8

    plot(p2sc, p2qsi[:,0], lw = 3, label = "Policy 2", color = '#FF1493')
    plot(p3sc, p3qsi[:,0], lw = 3, label = "Policy 3", color = "#FF0000") 
    plot(p4sc, p4qsi[:,0], lw = 3, label = "Policy 4", color = "#FFA500") 
    plot(p5sc, p5qsi[:,0], lw = 3, label = "Policy 5", color = "#FFFF00")
    plot(p55sc, p55qsi[:,0], lw = 3, label = "Policy 5.5", color = "#008000")
    plot(p6sc, p6qsi[:,0], lw = 3, label = "Policy 6", color = "#0000FF")
    plot(p7sc, p7qsi[:,0], lw = 3, label = "Policy 7", color = "#000000") 
    plot(p8sc, p8qsi[:,0], lw = 3, label = "Policy 8", color = "#800080") 

    plot(p2sc, p2qsi[:,1], '--', lw = 3, color = '#FF1493')
    plot(p3sc, p3qsi[:,1], '--', lw = 3, color = "#FF0000") 
    plot(p4sc, p4qsi[:,1], '--', lw = 3, color = "#FFA500") 
    plot(p5sc, p5qsi[:,1], '--', lw = 3, color = "#FFFF00")
    plot(p55sc, p55qsi[:,1], '--', lw = 3, color = "#008000")
    plot(p6sc, p6qsi[:,1], '--', lw = 3, color = "#0000FF")
    plot(p7sc, p7qsi[:,1], '--', lw = 3, color = "#000000") 
    plot(p8sc, p8qsi[:,1], '--', lw = 3, color = "#800080") 

    plot(p2sc, p2qsi[:,2], '-.', lw = 1, color = '#FF1493')
    plot(p3sc, p3qsi[:,2], '-.', lw = 1, color = "#FF0000") 
    plot(p4sc, p4qsi[:,2], '-.', lw = 1, color = "#FFA500") 
    plot(p5sc, p5qsi[:,2], '-.', lw = 1, color = "#FFFF00")
    plot(p55sc, p55qsi[:,2], '-.', lw = 1, color = "#008000")
    plot(p6sc, p6qsi[:,2], '-.', lw = 1, color = "#0000FF")
    plot(p7sc, p7qsi[:,2], '-.', lw = 1, color = "#000000") 
    plot(p8sc, p8qsi[:,2], '-.', lw = 1, color = "#800080") 

    plot(p2sc, p2qsi[:,3], ':', lw = 1, color = '#FF1493')
    plot(p3sc, p3qsi[:,3], ':', lw = 1, color = "#FF0000") 
    plot(p4sc, p4qsi[:,3], ':', lw = 1, color = "#FFA500") 
    plot(p5sc, p5qsi[:,3], ':', lw = 1, color = "#FFFF00")
    plot(p55sc, p55qsi[:,3], ':', lw = 1, color = "#008000")
    plot(p6sc, p6qsi[:,3], ':', lw = 1, color = "#0000FF")
    plot(p7sc, p7qsi[:,3], ':', lw = 1, color = "#000000") 
    plot(p8sc, p8qsi[:,3], ':', lw = 1, color = "#800080") 

#    plot(p2qso[:,-1], p2qb[:,4], '.', lw = 1, color = '#FF1493')
#    plot(p3qso[:,-1], p3qb[:,4], '.', lw = 1, color = "#FF0000") 
#    plot(p4qso[:,-1], p4qb[:,4], '.', lw = 1, color = "#FFA500") 
#    plot(p5qso[:,-1], p5qb[:,4], '.', lw = 1, color = "#FFFF00")
#    plot(p55qso[:,-1], p55qb[:,4], '.', lw = 1, color = "#008000")
#    plot(p6qso[:,-1], p6qb[:,4], '.', lw = 1, color = "#0000FF")
#    plot(p7qso[:,-1], p7qb[:,4], '.', lw = 1, color = "#000000") 
#    plot(p8qso[:,-1], p8qb[:,4], '.', lw = 1, color = "#800080") 

def plot_Q_X_vs_C_S(p2,p3,p4,p5,p55,p6,p7,p8):
    p2sc, p2me, p2qe, p2msi, p2qsi, p2mso, p2qso, p2mb, p2qb, p2msl = p2
    p3sc, p3me, p3qe, p3msi, p3qsi, p3mso, p3qso, p3mb, p3qb, p3msl = p3
    p4sc, p4me, p4qe, p4msi, p4qsi, p4mso, p4qso, p4mb, p4qb, p4msl = p4
    p5sc, p5me, p5qe, p5msi, p5qsi, p5mso, p5qso, p5mb, p5qb, p5msl = p5
    p55sc, p55me, p55qe, p55msi, p55qsi, p55mso, p55qso, p55mb, p55qb, p55msl = p55
    p6sc, p6me, p6qe, p6msi, p6qsi, p6mso, p6qso, p6mb, p6qb, p6msl = p6
    p7sc, p7me, p7qe, p7msi, p7qsi, p7mso, p7qso, p7mb, p7qb, p7msl = p7
    p8sc, p8me, p8qe, p8msi, p8qsi, p8mso, p8qso, p8mb, p8qb, p8msl = p8

    plot(p2sc, p2qe[:,0], lw = 3, label = "Policy 2", color = '#FF1493')
    plot(p3sc, p3qe[:,0], lw = 3, label = "Policy 3", color = "#FF0000") 
    plot(p4sc, p4qe[:,0], lw = 3, label = "Policy 4", color = "#FFA500") 
    plot(p5sc, p5qe[:,0], lw = 3, label = "Policy 5", color = "#FFFF00")
    plot(p55sc, p55qe[:,0], lw = 3, label = "Policy 5.5", color = "#008000")
    plot(p6sc, p6qe[:,0], lw = 3, label = "Policy 6", color = "#0000FF")
    plot(p7sc, p7qe[:,0], lw = 3, label = "Policy 7", color = "#000000") 
    plot(p8sc, p8qe[:,0], lw = 3, label = "Policy 8", color = "#800080") 

    plot(p2sc, p2qe[:,1], '--', lw = 3, color = '#FF1493')
    plot(p3sc, p3qe[:,1], '--', lw = 3, color = "#FF0000") 
    plot(p4sc, p4qe[:,1], '--', lw = 3, color = "#FFA500") 
    plot(p5sc, p5qe[:,1], '--', lw = 3, color = "#FFFF00")
    plot(p55sc, p55qe[:,1], '--', lw = 3, color = "#008000")
    plot(p6sc, p6qe[:,1], '--', lw = 3, color = "#0000FF")
    plot(p7sc, p7qe[:,1], '--', lw = 3, color = "#000000") 
    plot(p8sc, p8qe[:,1], '--', lw = 3, color = "#800080") 

    plot(p2sc, p2qe[:,2], '-.', lw = 1, color = '#FF1493')
    plot(p3sc, p3qe[:,2], '-.', lw = 1, color = "#FF0000") 
    plot(p4sc, p4qe[:,2], '-.', lw = 1, color = "#FFA500") 
    plot(p5sc, p5qe[:,2], '-.', lw = 1, color = "#FFFF00")
    plot(p55sc, p55qe[:,2], '-.', lw = 1, color = "#008000")
    plot(p6sc, p6qe[:,2], '-.', lw = 1, color = "#0000FF")
    plot(p7sc, p7qe[:,2], '-.', lw = 1, color = "#000000") 
    plot(p8sc, p8qe[:,2], '-.', lw = 1, color = "#800080") 

    plot(p2sc, p2qe[:,3], ':', lw = 1, color = '#FF1493')
    plot(p3sc, p3qe[:,3], ':', lw = 1, color = "#FF0000") 
    plot(p4sc, p4qe[:,3], ':', lw = 1, color = "#FFA500") 
    plot(p5sc, p5qe[:,3], ':', lw = 1, color = "#FFFF00")
    plot(p55sc, p55qe[:,3], ':', lw = 1, color = "#008000")
    plot(p6sc, p6qe[:,3], ':', lw = 1, color = "#0000FF")
    plot(p7sc, p7qe[:,3], ':', lw = 1, color = "#000000") 
    plot(p8sc, p8qe[:,3], ':', lw = 1, color = "#800080") 

#    plot(p2qso[:,-1], p2qb[:,4], '.', lw = 1, color = '#FF1493')
#    plot(p3qso[:,-1], p3qb[:,4], '.', lw = 1, color = "#FF0000") 
#    plot(p4qso[:,-1], p4qb[:,4], '.', lw = 1, color = "#FFA500") 
#    plot(p5qso[:,-1], p5qb[:,4], '.', lw = 1, color = "#FFFF00")
#    plot(p55qso[:,-1], p55qb[:,4], '.', lw = 1, color = "#008000")
#    plot(p6qso[:,-1], p6qb[:,4], '.', lw = 1, color = "#0000FF")
#    plot(p7qso[:,-1], p7qb[:,4], '.', lw = 1, color = "#000000") 
#    plot(p8qso[:,-1], p8qb[:,4], '.', lw = 1, color = "#800080") 

def plot_Q_B_vs_K_S_in(p2,p3,p4,p5,p55,p6,p7,p8):
    p2sc, p2me, p2qe, p2msi, p2qsi, p2mso, p2qso, p2mb, p2qb, p2msl = p2
    p3sc, p3me, p3qe, p3msi, p3qsi, p3mso, p3qso, p3mb, p3qb, p3msl = p3
    p4sc, p4me, p4qe, p4msi, p4qsi, p4mso, p4qso, p4mb, p4qb, p4msl = p4
    p5sc, p5me, p5qe, p5msi, p5qsi, p5mso, p5qso, p5mb, p5qb, p5msl = p5
    p55sc, p55me, p55qe, p55msi, p55qsi, p55mso, p55qso, p55mb, p55qb, p55msl = p55
    p6sc, p6me, p6qe, p6msi, p6qsi, p6mso, p6qso, p6mb, p6qb, p6msl = p6
    p7sc, p7me, p7qe, p7msi, p7qsi, p7mso, p7qso, p7mb, p7qb, p7msl = p7
    p8sc, p8me, p8qe, p8msi, p8qsi, p8mso, p8qso, p8mb, p8qb, p8msl = p8

    plot(p2qsi[:,-1], p2qb[:,0], lw = 3, label = "Policy 2", color = '#FF1493')
    plot(p3qsi[:,-1], p3qb[:,0], lw = 3, label = "Policy 3", color = "#FF0000") 
    plot(p4qsi[:,-1], p4qb[:,0], lw = 3, label = "Policy 4", color = "#FFA500") 
    plot(p5qsi[:,-1], p5qb[:,0], lw = 3, label = "Policy 5", color = "#FFFF00")
    plot(p55qso[:,-1], p55qb[:,0], lw = 3, label = "Policy 5.5", color = "#008000")
    plot(p6qsi[:,-1], p6qb[:,0], lw = 3, label = "Policy 6", color = "#0000FF")
    plot(p7qsi[:,-1], p7qb[:,0], lw = 3, label = "Policy 7", color = "#000000") 
    plot(p8qsi[:,-1], p8qb[:,0], lw = 3, label = "Policy 8", color = "#800080") 

    plot(p2qsi[:,-1], p2qb[:,1], '--', lw = 3, color = '#FF1493')
    plot(p3qsi[:,-1], p3qb[:,1], '--', lw = 3, color = "#FF0000") 
    plot(p4qsi[:,-1], p4qb[:,1], '--', lw = 3, color = "#FFA500") 
    plot(p5qsi[:,-1], p5qb[:,1], '--', lw = 3, color = "#FFFF00")
    plot(p55qsi[:,-1], p55qb[:,1], '--', lw = 3, color = "#008000")
    plot(p6qsi[:,-1], p6qb[:,1], '--', lw = 3, color = "#0000FF")
    plot(p7qsi[:,-1], p7qb[:,1], '--', lw = 3, color = "#000000") 
    plot(p8qsi[:,-1], p8qb[:,1], '--', lw = 3, color = "#800080") 

    plot(p2qsi[:,-1], p2qb[:,2], '-.', lw = 1, color = '#FF1493')
    plot(p3qsi[:,-1], p3qb[:,2], '-.', lw = 1, color = "#FF0000") 
    plot(p4qsi[:,-1], p4qb[:,2], '-.', lw = 1, color = "#FFA500") 
    plot(p5qsi[:,-1], p5qb[:,2], '-.', lw = 1, color = "#FFFF00")
    plot(p55qsi[:,-1], p55qb[:,2], '-.', lw = 1, color = "#008000")
    plot(p6qsi[:,-1], p6qb[:,2], '-.', lw = 1, color = "#0000FF")
    plot(p7qsi[:,-1], p7qb[:,2], '-.', lw = 1, color = "#000000") 
    plot(p8qsi[:,-1], p8qb[:,2], '-.', lw = 1, color = "#800080") 

    plot(p2qsi[:,-1], p2qb[:,3], ':', lw = 1, color = '#FF1493')
    plot(p3qsi[:,-1], p3qb[:,3], ':', lw = 1, color = "#FF0000") 
    plot(p4qsi[:,-1], p4qb[:,3], ':', lw = 1, color = "#FFA500") 
    plot(p5qsi[:,-1], p5qb[:,3], ':', lw = 1, color = "#FFFF00")
    plot(p55qsi[:,-1], p55qb[:,3], ':', lw = 1, color = "#008000")
    plot(p6qsi[:,-1], p6qb[:,3], ':', lw = 1, color = "#0000FF")
    plot(p7qsi[:,-1], p7qb[:,3], ':', lw = 1, color = "#000000") 
    plot(p8qsi[:,-1], p8qb[:,3], ':', lw = 1, color = "#800080") 

#    plot(p2qso[:,-1], p2qb[:,4], '.', lw = 1, color = '#FF1493')
#    plot(p3qso[:,-1], p3qb[:,4], '.', lw = 1, color = "#FF0000") 
#    plot(p4qso[:,-1], p4qb[:,4], '.', lw = 1, color = "#FFA500") 
#    plot(p5qso[:,-1], p5qb[:,4], '.', lw = 1, color = "#FFFF00")
#    plot(p55qso[:,-1], p55qb[:,4], '.', lw = 1, color = "#008000")
#    plot(p6qso[:,-1], p6qb[:,4], '.', lw = 1, color = "#0000FF")
#    plot(p7qso[:,-1], p7qb[:,4], '.', lw = 1, color = "#000000") 
#    plot(p8qso[:,-1], p8qb[:,4], '.', lw = 1, color = "#800080") 

def plot_C_S_vs_K_B(p2,p3,p4,p5,p55,p6,p7,p8):
    p2sc, p2me, p2qe, p2msi, p2qsi, p2mso, p2qso, p2mb, p2qb, p2msl = p2
    p3sc, p3me, p3qe, p3msi, p3qsi, p3mso, p3qso, p3mb, p3qb, p3msl = p3
    p4sc, p4me, p4qe, p4msi, p4qsi, p4mso, p4qso, p4mb, p4qb, p4msl = p4
    p5sc, p5me, p5qe, p5msi, p5qsi, p5mso, p5qso, p5mb, p5qb, p5msl = p5
    p55sc, p55me, p55qe, p55msi, p55qsi, p55mso, p55qso, p55mb, p55qb, p55msl = p55
    p6sc, p6me, p6qe, p6msi, p6qsi, p6mso, p6qso, p6mb, p6qb, p6msl = p6
    p7sc, p7me, p7qe, p7msi, p7qsi, p7mso, p7qso, p7mb, p7qb, p7msl = p7
    p8sc, p8me, p8qe, p8msi, p8qsi, p8mso, p8qso, p8mb, p8qb, p8msl = p8

    plot(p2qb[:,-1], p2sc, lw = 3, label = "Policy 2", color = '#FF1493')
    plot(p3qb[:,-1], p3sc, lw = 3, label = "Policy 3", color = "#FF0000") 
    plot(p4qb[:,-1], p4sc, lw = 3, label = "Policy 4", color = "#FFA500") 
    plot(p5qb[:,-1], p5sc, lw = 3, label = "Policy 5", color = "#FFFF00")
    plot(p55qb[:,-1], p55sc, lw = 3, label = "Policy 5.5", color = "#008000")
    plot(p6qb[:,-1], p6sc, lw = 3, label = "Policy 6", color = "#0000FF")
    plot(p7qb[:,-1], p7sc, lw = 3, label = "Policy 7", color = "#000000") 
    plot(p8qb[:,-1], p8sc, lw = 3, label = "Policy 8", color = "#800080") 

def plot_Q_S_out_vs_K_B(p2,p3,p4,p5,p55,p6,p7,p8):
    p2sc, p2me, p2qe, p2msi, p2qsi, p2mso, p2qso, p2mb, p2qb, p2msl = p2
    p3sc, p3me, p3qe, p3msi, p3qsi, p3mso, p3qso, p3mb, p3qb, p3msl = p3
    p4sc, p4me, p4qe, p4msi, p4qsi, p4mso, p4qso, p4mb, p4qb, p4msl = p4
    p5sc, p5me, p5qe, p5msi, p5qsi, p5mso, p5qso, p5mb, p5qb, p5msl = p5
    p55sc, p55me, p55qe, p55msi, p55qsi, p55mso, p55qso, p55mb, p55qb, p55msl = p55
    p6sc, p6me, p6qe, p6msi, p6qsi, p6mso, p6qso, p6mb, p6qb, p6msl = p6
    p7sc, p7me, p7qe, p7msi, p7qsi, p7mso, p7qso, p7mb, p7qb, p7msl = p7
    p8sc, p8me, p8qe, p8msi, p8qsi, p8mso, p8qso, p8mb, p8qb, p8msl = p8

    plot(p2qb[:,-1], p2qso[:,0], lw = 3, label = "Policy 2", color = '#FF1493')
    plot(p3qb[:,-1], p3qso[:,0], lw = 3, label = "Policy 3", color = "#FF0000") 
    plot(p4qb[:,-1], p4qso[:,0], lw = 3, label = "Policy 4", color = "#FFA500") 
    plot(p5qb[:,-1], p5qso[:,0], lw = 3, label = "Policy 5", color = "#FFFF00")
    plot(p55qb[:,-1], p55qso[:,0], lw = 3, label = "Policy 5.5", color = "#008000")
    plot(p6qb[:,-1], p6qso[:,0], lw = 3, label = "Policy 6", color = "#0000FF")
    plot(p7qb[:,-1], p7qso[:,0], lw = 3, label = "Policy 7", color = "#000000") 
    plot(p8qb[:,-1], p8qso[:,0], lw = 3, label = "Policy 8", color = "#800080") 

    plot(p2qb[:,-1], p2qso[:,1], '--', lw = 3, color = '#FF1493')
    plot(p3qb[:,-1], p3qso[:,1], '--', lw = 3, color = "#FF0000") 
    plot(p4qb[:,-1], p4qso[:,1], '--', lw = 3, color = "#FFA500") 
    plot(p5qb[:,-1], p5qso[:,1], '--', lw = 3, color = "#FFFF00")
    plot(p55qb[:,-1], p55qso[:,1], '--', lw = 3, color = "#008000")
    plot(p6qb[:,-1], p6qso[:,1], '--', lw = 3, color = "#0000FF")
    plot(p7qb[:,-1], p7qso[:,1], '--', lw = 3, color = "#000000") 
    plot(p8qb[:,-1], p8qso[:,1], '--', lw = 3, color = "#800080") 

    plot(p2qb[:,-1], p2qso[:,2], '-.', lw = 1, color = '#FF1493')
    plot(p3qb[:,-1], p3qso[:,2], '-.', lw = 1, color = "#FF0000") 
    plot(p4qb[:,-1], p4qso[:,2], '-.', lw = 1, color = "#FFA500") 
    plot(p5qb[:,-1], p5qso[:,2], '-.', lw = 1, color = "#FFFF00")
    plot(p55qb[:,-1], p55qso[:,2], '-.', lw = 1, color = "#008000")
    plot(p6qb[:,-1], p6qso[:,2], '-.', lw = 1, color = "#0000FF")
    plot(p7qb[:,-1], p7qso[:,2], '-.', lw = 1, color = "#000000") 
    plot(p8qb[:,-1], p8qso[:,2], '-.', lw = 1, color = "#800080") 

    plot(p2qb[:,-1], p2qso[:,3], ':', lw = 1, color = '#FF1493')
    plot(p3qb[:,-1], p3qso[:,3], ':', lw = 1, color = "#FF0000") 
    plot(p4qb[:,-1], p4qso[:,3], ':', lw = 1, color = "#FFA500") 
    plot(p5qb[:,-1], p5qso[:,3], ':', lw = 1, color = "#FFFF00")
    plot(p55qb[:,-1], p55qso[:,3], ':', lw = 1, color = "#008000")
    plot(p6qb[:,-1], p6qso[:,3], ':', lw = 1, color = "#0000FF")
    plot(p7qb[:,-1], p7qso[:,3], ':', lw = 1, color = "#000000") 
    plot(p8qb[:,-1], p8qso[:,3], ':', lw = 1, color = "#800080") 

#    plot(p2qso[:,-1], p2qb[:,4], '.', lw = 1, color = '#FF1493')
#    plot(p3qso[:,-1], p3qb[:,4], '.', lw = 1, color = "#FF0000") 
#    plot(p4qso[:,-1], p4qb[:,4], '.', lw = 1, color = "#FFA500") 
#    plot(p5qso[:,-1], p5qb[:,4], '.', lw = 1, color = "#FFFF00")
#    plot(p55qso[:,-1], p55qb[:,4], '.', lw = 1, color = "#008000")
#    plot(p6qso[:,-1], p6qb[:,4], '.', lw = 1, color = "#0000FF")
#    plot(p7qso[:,-1], p7qb[:,4], '.', lw = 1, color = "#000000") 
#    plot(p8qso[:,-1], p8qb[:,4], '.', lw = 1, color = "#800080") 

def get_gammarange_data(gammarange, alpha, eta_in = 1., eta_out = 1., policy = 2, quantile = .99, storage_capacity = NaN):
    N = len(gammarange)
    quantiles = array(quantile)
    L = quantiles.size
    print L
    sc = empty(N)
    me = empty(N)
    qe = empty([N,L])
    msi = empty(N)
    qsi = empty([N,L])
    mso = empty(N)
    qso = empty([N,L])
    mb = empty(N)
    qb = empty([N,L])
    msl = empty(N)
    for i in arange(len(gammarange)):
        sc[i], me[i], qe[i], msi[i], qsi[i], mso[i], qso[i], mb[i], qb[i], msl[i] = get_storage_data(gammarange[i], alpha, storage_capacity, eta_in, eta_out, policy, quantile)
    return sc, me, qe, msi, qsi, mso, qso, mb, qb, msl
    
def get_contours(gammarange, alpharange, load, wind, solar, eta_in = 1., eta_out = 1.):
    N = len(gammarange)
    L = len(alpharange)
    sc_contour = empty([N,L])
    b_contour = empty([N,L])
    length = len(load)
    for i in arange(N):
        for j in arange(L):
            sc_contour[i,j], a, b, c, d, e, f, b_contour[i,j], g, h = get_storage_data(NaN, 1, NaN, eta_in, eta_out, 2, 1., gammarange[i]*(alpharange[j]*wind + (1-alpharange[j])*solar) - load)
    return sc_contour, b_contour

def get_fixed_gamma_contours(gamma, storagerange, alpharange, load, wind, solar, eta_in = 1., eta_out = 1.):
    N = len(storagerange)
    L = len(alpharange)
    sc_contour = empty([N,L])
    b_contour = empty([N,L])
    length = len(load)
    for i in arange(N):
        for j in arange(L):
            sc_contour[i,j], a, b, c, d, e, f, b_contour[i,j], g, h = get_storage_data(NaN, 1, storagerange[i], eta_in, eta_out, 2, 1., gamma*(alpharange[j]*wind + (1-alpharange[j])*solar) - load)
    return sc_contour, b_contour

def get_contours_given_max_storage(gammarange, alpharange, maxstorage, load, wind, solar, eta_in = 1., eta_out = 1.):
    N = len(gammarange)
    L = len(alpharange)
    sc_contour = empty([N,L])
    b_contour = empty([N,L])
    length = len(load)
    for i in arange(N):
        for j in arange(L):
            sc_contour[i,j], a, b, c, d, e, f, b_contour[i,j], g, h = get_storage_data(NaN, 1, maxstorage, eta_in, eta_out, 2, 1., gammarange[i]*(alpharange[j]*wind + (1-alpharange[j])*solar) - load)
    return sc_contour, b_contour

def get_contours_for_limited_in_and_output_capacities(mismatch, eta_in = .6, eta_out = .6, N = 30, maxstorage = NaN, L = None):
    if L is None:
        L = N
    b_contour = empty([N,L])
    sc_contour = empty([N,L])
    inrange = linspace(0, mismatch.max(), N)
    outrange = linspace(0, -mismatch.min(), L)
    for i in arange(N):
        print "i: ", i, inrange[i]
        constant.policy_5_in_capacity = inrange[i]
        sc_contour[i,0] = 0
        b_contour[i,0] = neg(mismatch).mean()
        for j in arange(1,L):
            print "j: ", j, outrange[i]
            constant.policy_5_out_capacity = outrange[j]
            sc_contour[i,j], a, b, c, d, e, f, b_contour[i,j], g, h = get_storage_data(NaN, 1, maxstorage, eta_in, eta_out, 5, 1., mismatch)
    return inrange, outrange, b_contour, sc_contour

def loadLoadSolarWind(state='CA', txtlabel='', path='./'):

    filename = 'LoadWindSolar_' + state + txtlabel
    table = np.load(path+filename+'.npy')

    return table[0]

def get_CA_mismatch(gamma, alpha):
    data = loadLoadSolarWind()
    load = data['load']
    solar = data['solar']
    wind = data['wind']
 
    return gamma * load.mean() * (alpha * wind / wind.mean() + (1 - alpha) * solar / solar.mean()) - load

def get_CA_load_wind_solar():
    data = loadLoadSolarWind()
    load = data['load']
    solar = data['solar']
    wind = data['wind']
    
    return load, wind * load.mean() / wind.mean(), solar * load.mean() / solar.mean()

def find_alpha_min_balancing(gamma, storage_capacity, eta_in = 1., eta_out = 1., load = get_Load(countries.Region), wind = get_Wind(countries.Region), solar = get_PV(countries.Region)):
    def find_balancing(alpha):
        s, a, b, c, d, e, f, balancing, g, h = get_storage_data(NaN, 1, storage_capacity, eta_in, eta_out, 2, 1., gamma*(alpha*wind + (1-alpha)*solar) - load)
        return balancing
    return fminbound(find_balancing, 0, 1)

def find_alpha_min_storage(gamma, eta_in = 1., eta_out = 1., load = get_Load(countries.Region), wind = get_Wind(countries.Region), solar = get_PV(countries.Region)):
    def find_storage(alpha):
        dmmy, storage = get_policy_2_storage(gamma*(alpha*wind + (1-alpha)*solar) - load, eta_in, eta_out, NaN)
        return storage
    return fminbound(find_storage, 0, 1)

def find_K_B_given_storage(gamma, storage_size, load, wind, solar):
    def find_storage(K_B):
        constant.policy_7_capacity = K_B
        dummy, storage = get_policy_7_storage(gamma * (.6 * wind + .4 * solar) - load, 1., 1., NaN)
        print mquantiles(-dummy, 1.)
        return storage - storage_size
    return brentq(find_storage, 0, 1.4)

def find_K_B_given_storage2(gamma, storage_size, load, wind, solar):
    mismatch = gamma * (.6 * wind + .4 * solar) - load
    def find_storage(K_B):
        constant.policy_7_capacity = K_B
        dummy, storage = get_policy_7_storage(gamma * (.6 * wind + .4 * solar) - load, 1., 1., NaN)
        print storage, -min(mismatch), K_B
        print storage_size - storage
        return abs(storage_size - storage)
    return fminbound(find_storage, 0, -min(mismatch))

def find_K_B_that_balances(gamma, load, wind, solar):
    mismatch = gamma * (.6 * wind + .4 * solar) - load
    def find_equilibrium(K_B):
        storage_part = pos(mismatch) - pos(neg(mismatch) - K_B)
        return storage_part.sum()
    return brentq(find_equilibrium, 0, -min(mismatch))

def find_K_B_given_K_B(gamma, storage_size, load, wind, solar):
    def find_K_B(K_B):
        constant.policy_7_capacity = K_B
        #d1,d2,d3,d4,d5,d6,d7,d8,K_B2, d9 = get_storage_data(gamma, .6, storage_size, 1.,1.,7, 1.)
        mismatch = gamma * (.6 * wind + .4 * solar) - load
        new_mismatch, storage = get_policy_7_storage(mismatch, 1., 1., storage_size)
        K_B2 = mquantiles(-new_mismatch, 1.)
        print storage, storage_size, K_B, K_B2
        return K_B - K_B2
    return fminbound(find_K_B, 0, 1.4)

def plot_graf_med_indsats(mini, over, under, fast, figno, sc, steps):
    figure(figno)
    clf()
    ax1 = axes()
    ax2 = axes([0.58,0.6,0.275,0.275])
    plot(linspace(0,sc,steps)/8766,mini,'k', lw = 1)
    plot(linspace(0,sc,steps)/8766,over,'k--', lw = 1)
    plot(linspace(0,sc,steps)/8766,under,'k--', lw = 1)
    plot(linspace(0,sc,steps)/8766,fast * ones(steps),'k:', lw = .5)
    axis("tight")
    axis(ymin = 0, ymax = 1)
    xlabel(r"[av.y.l.]")
    ylabel(r"$\alpha_W$")
    oom = int(log10(sc/8766))
    size = int(sc/8766 * 10**(1-oom))
    xticks((linspace(0,size * 10**(oom-1),3)))
    axes(ax1)
    plot(linspace(0,sc,steps),mini,'k', lw = 1)
    plot(linspace(0,sc,steps),over,'k--', lw = 1)
    plot(linspace(0,sc,steps),under,'k--', lw = 1)
    plot(linspace(0,sc,steps),fast * ones(steps),'k:', lw = .5)
    axis(xmax = 24)
    xticks(arange(5)*6)
    xlabel("Storage energy capacity [av.h.l.]")
    ylabel(r"Wind fraction $\alpha_W$")

def get_policy_2_3_4_data(mismatch, N, eta_in, eta_out, quantile = .99):
    #policy 2:
    p2sc, p2me, p2qe, p2msi, p2qsi, p2mso, p2qso, p2mb, p2qb, p2msl = get_storages_data(NaN, 1, N, eta_in, eta_out, 2, quantile, concatenate([[0],exp(linspace(-31/10.,(N-1)/10.,(N-1)))/exp((N-1)/10.)]), mismatch)
    #policy 3:
    p3sc, p3me, p3qe, p3msi, p3qsi, p3mso, p3qso, p3mb, p3qb, p3msl = get_storages_data(NaN, 1, N, eta_in, eta_out, 3, quantile, concatenate([[0],exp(linspace(-31/10.,(N-1)/10.,(N-1)))/exp((N-1)/10.)]), mismatch)
    #policy 4:
    p4sc, p4me, p4qe, p4msi, p4qsi, p4mso, p4qso, p4mb, p4qb, p4msl = get_storages_data(NaN, 1, N, eta_in, eta_out, 4, quantile, concatenate([[0],exp(linspace(-31/10.,(N-1)/10.,(N-1)))/exp((N-1)/10.)]), mismatch)
    p2 = (p2sc, p2me, p2qe, p2msi, p2qsi, p2mso, p2qso, p2mb, p2qb, p2msl)
    p3 = (p3sc, p3me, p3qe, p3msi, p3qsi, p3mso, p3qso, p3mb, p3qb, p3msl)
    p4 = (p4sc, p4me, p4qe, p4msi, p4qsi, p4mso, p4qso, p4mb, p4qb, p4msl)
    return p2, p3, p4

def plot_MP_fig3a(alphagammacontour, alphastoragecontour, alphacontour, figno = 3):
    fig = figure(figno)
    fig.set_size_inches([3.425, 6.125*3.425/8.125])
    clf()
    ax1 = gca()
    ax1.set_position([.16,.15, .69, .8])
    contourf(alphagammacontour, alphastoragecontour + 1e-5, alphacontour, 300, cmap = cm.jet_r)
    contourf(alphagammacontour, alphastoragecontour + 1e-5, alphacontour, 300, cmap = cm.jet_r)
    contourf(alphagammacontour, alphastoragecontour + 1e-5, alphacontour, 300, cmap = cm.jet_r)
    axis("tight")
    setp(gca(), yscale = "log")
    axis(ymin = .5)
    cb = colorbar(ticks = linspace(.25,.8,12))
    cb.set_label(r"Wind fraction $\alpha_W$")
    xlabel(r"Renewable power generation $\gamma$")
    ylabel(r"Storage energy capacity [av.h.l.]")

def plot_MP_fig3b(alphagammacontour, alphastoragecontour, alphabalancingcontour, g075bcont, figno = 3):
    fig = figure(figno)
    fig.set_size_inches([3.425, 6.125*3.425/8.125])
    clf()
    ax1 = gca()
    ax1.set_position([.16,.15, .69, .8])
    contourf(alphagammacontour, alphastoragecontour + 1e-5, alphabalancingcontour, 300)
    clim([0,g075bcont.max()])
    contourf(alphagammacontour, alphastoragecontour + 1e-5, alphabalancingcontour, 300)
    clim([0,g075bcont.max()])
    contourf(alphagammacontour, alphastoragecontour + 1e-5, alphabalancingcontour, 300)
    clim([0,g075bcont.max()])
    axis("tight")
    axis(ymin = .5)
    setp(ax1, yscale = 'log')
    cb = colorbar()
    cb.set_label(r"Balancing energy $E_B$")
    xlabel(r"Renewable power generation $\gamma$")
    ylabel("Storage energy capacity [av.h.l.]")

def plot_MP_fig3d(alphagammacontour, alphastoragecontour, fixedbalancingcontour, alphabalancingcontour, g075bcont, figno = 3):
    fig = figure(figno)
    fig.set_size_inches([3.425, 6.125*3.425/8.125])
    clf()
    ax1 = gca()
    ax1.set_position([.16,.15, .69, .8])
    contourf(alphagammacontour, alphastoragecontour + 1e-5, fixedbalancingcontour - alphabalancingcontour, 300, cmap = cm.PuRd)
    clim([0,g075bcont.max()])
    contourf(alphagammacontour, alphastoragecontour + 1e-5, fixedbalancingcontour - alphabalancingcontour, 300, cmap = cm.PuRd)
    clim([0,g075bcont.max()])
    contourf(alphagammacontour, alphastoragecontour + 1e-5, fixedbalancingcontour - alphabalancingcontour, 300, cmap = cm.PuRd)
    clim([0,g075bcont.max()])
    axis("tight")
    axis(ymin = .5)
    setp(ax1, yscale = 'log')
    cb = colorbar()
    cb.set_label(r"Balancing energy difference")
    xlabel(r"Renewable power generation $\gamma$")
    ylabel("Storage energy capacity [av.h.l.]")
    contour(alphagammacontour,alphastoragecontour,fixedbalancingcontour - alphabalancingcontour, [.01, .005, .0025], colors = ('k','.5', '.75'))


    
def plot_MP_fig9b(mini1, mini2, mini3, over1, over2, over3, under1, under2, under3, sc1, sc2, sc3, fast = .6, figno = 9):
    steps1 = mini1.size
    steps2 = mini2.size
    steps3 = mini3.size
    fig = figure(figno)
    clf()
    fig.set_size_inches([3.425, 6.125*3.425/8.125])
    ax1 = axes()
    ax1.set_position([.15,.17,.84,.8])
#    ax2 = axes([0.28,0.3,0.3,0.25])
#    axes(ax1)
    matplotlib.rcParams['font.size'] = 10
    plot(linspace(0,sc1,steps1),mini1,'r', lw = 2)
    plot(linspace(0,sc1,steps1),over1, ':', color = 'r', lw = 1)
    plot(linspace(0,sc1,steps1),under1,':', color = 'r', lw = 1)
    plot(linspace(0,sc2,steps2),mini2,'g', lw = 2)
    plot(linspace(0,sc2,steps2),over2, ':', color = 'g', lw = 1)
    plot(linspace(0,sc2,steps2),under2,':', color = 'g', lw = 1)
    plot(linspace(0,sc3,steps3),mini3,'b', lw = 2)
    plot(linspace(0,sc3,steps3),over3, ':', color = 'b', lw = 1)
    plot(linspace(0,sc3,steps3),under3,':', color = 'b', lw = 1)
    plot(linspace(0,sc2,steps2),fast * ones(steps2), 'k', lw = 1, alpha = .5)
#    axis("tight")
    xlabel(r"Storage energy capacity [av.y.l.]")
    ylabel(r"Wind fraction $\alpha_W$")
    setp(ax1, xscale = 'log')
    axis([9,1250,0,1])
    #oom = int(log10(sc/8766))
    #size = int(sc/8766 * 10**(1-oom))
    #xticks((linspace(0,size * 10**(oom-1),4)))
#    axes(ax2)
#    plot(linspace(0,sc,steps),mini,'k', lw = 1)
#    plot(linspace(0,sc,steps),over, ':', color = ".25", lw = 1)
#    plot(linspace(0,sc,steps),under, ':', color = ".25", lw = 1)
#    plot(linspace(0,sc,steps),fast * ones(steps),'k', lw = 1, alpha = .5)
#    axis(ymin = .3, xmax = 24)
#    yticks(arange(.4,.9,.2))
#    xticks(arange(5)*6)
#    xlabel("[av.h.l.]")

def plot_MP_fig9c(mini1, mini2, mini3, over1, over2, over3, under1, under2, under3, sc1, sc2, sc3, fast = .6, figno = 9):
    steps1 = mini1.size
    steps2 = mini2.size
    steps3 = mini3.size
    fig = figure(figno)
    clf()
    fig.set_size_inches([3.425, 6.125*3.425/8.125])
    ax1 = axes()
    ax1.set_position([.17,.15,.79,.81])
#    ax2 = axes([0.28,0.3,0.3,0.25])
#    axes(ax1)
    matplotlib.rcParams['font.size'] = 10
    plot(linspace(0,sc1,steps1),mini1,'r', lw = 2, label = r"$\gamma = 0.75$")
    plot(linspace(0,sc1,steps1),over1, ':', color = 'r', lw = 1)
    plot(linspace(0,sc1,steps1),under1,':', color = 'r', lw = 1)
    plot(linspace(0,sc2,steps2),mini2,'g', lw = 2, label = r"$\gamma = 1.00$")
    plot(linspace(0,sc2,steps2),over2, ':', color = 'g', lw = 1)
    plot(linspace(0,sc2,steps2),under2,':', color = 'g', lw = 1)
    plot(linspace(0,sc3,steps3),mini3,'b', lw = 2, label = r"$\gamma = 1.25$")
    plot(linspace(0,sc3,steps3),over3, ':', color = 'b', lw = 1)
    plot(linspace(0,sc3,steps3),under3,':', color = 'b', lw = 1)
    plot(linspace(0,sc2,steps2),fast * ones(steps2), 'k', lw = 1, alpha = .5)
#    axis("tight")
    xlabel(r"Storage energy capacity $C_S$ [av.h.l.]")
    ylabel(r"Wind fraction $\alpha_W$")
    xticks((0,6,12,18,24))
    leg = legend(loc = 3, labelspacing = -.1)#, bbox_to_anchor = (0.06,0,1,1))
    setp(leg.get_texts(), fontsize = "small")
    axvline(x=6, ls='--', color = 'k')
    
#    setp(ax1, xscale = 'log')
    axis([0,24,0,1])
    #oom = int(log10(sc/8766))
    #size = int(sc/8766 * 10**(1-oom))
    #xticks((linspace(0,size * 10**(oom-1),4)))
#    axes(ax2)
#    plot(linspace(0,sc,steps),mini,'k', lw = 1)
#    plot(linspace(0,sc,steps),over, ':', color = ".25", lw = 1)
#    plot(linspace(0,sc,steps),under, ':', color = ".25", lw = 1)
#    plot(linspace(0,sc,steps),fast * ones(steps),'k', lw = 1, alpha = .5)
#    axis(ymin = .3, xmax = 24)
#    yticks(arange(.4,.9,.2))
#    xticks(arange(5)*6)
#    xlabel("[av.h.l.]")
   
def plot_MP_fig1(gammarange, so06, so08, sooptimal, figno = 1):
    matplotlib.rcParams['font.size'] = 10
    fig = figure(figno)
    clf()
    fig.set_size_inches([3.425, 6.125*3.425/8.125])
    ax1 = axes()
    ax1.set_position([.17,.15,.79,.81])
#    ax2 = axes([0.675,0.605,0.3,0.3])
#    axes(ax1)
    plot(gammarange, so06/8766, 'b', alpha = 1, lw = 1, label = r"$\alpha_W =\, 0.60$")
    plot(gammarange, so08/8766, 'r', alpha = .75, lw = 1, label = r"$\alpha_W =\, 0.80$")
    plot(gammarange, sooptimal/8766, '.5', lw = 1, label = r"Optimal $\alpha_W$")
    xlabel(r"Average RES power generation factor $\gamma$")
    ax1.set_ylabel(r"Storage energy capacity $C_S$ [av.y.l.]", y = .45)
#    ylabel("Storage energy capacity [av.y.l.]")
    leg = legend(loc = 1, labelspacing = -.1, bbox_to_anchor = (0.06,0,1,1))
    setp(leg.get_texts(), fontsize = "small")
#    axes(ax2)
#    plot(gammarange, sooptimal, 'k', lw = 1, alpha = .15, label = r"optimal $\alpha_W$")
    plot(gammarange, so06/8766, 'b', lw = 1, alpha = 1, label = r"$\alpha_W = 0.60$")
    plot(gammarange, so08/8766, 'r', lw = 1, alpha = .75, label = r"$\alpha_W = 0.80$")
#    xlabel(r"$\gamma$")
#    ylabel(r"[av.h.l.]")
#    setp(ax2, yscale = "log")
#    axis("tight")
#    axis(xmin = .85,xmax = 1.15)#, ymin = 20)
#    xticks(arange(0.9,1.1,.1))
#    yticks(arange(0,1210,400))

def plot_MP_fig2a(gammarange, bo06, bo08, booptimal, basebalancing, figno = 1):
    matplotlib.rcParams['font.size'] = 10
    fig = figure(figno)
    clf()
    fig.set_size_inches([3.425, 6.125*3.425/8.125])
    ax1 = axes()
    ax1.set_position([.17,.15,.79,.81])
    plot(gammarange, bo06, 'k', alpha = 1, lw = 1, label = r"$\alpha_W =\, 0.60$")
    plot(gammarange, bo08, 'k', alpha = .5, lw = 1, label = r"$\alpha_W =\, 0.80$")
    plot(gammarange, booptimal, 'k', lw = 1, color = ".5", label = r"Optimal $\alpha_W$")
    plot(gammarange, bo06, 'k', alpha = 1, lw = 1)
    plot(gammarange, bo08, 'k', alpha = .5, lw = 1)
    plot(gammarange, basebalancing, 'k:', lw = 1)
    xlabel(r"Average RES power generation $\gamma$")
    ylabel("Balancing energy fraction $E_B^{\textup{excess}}$")
    leg = legend(loc = 0, labelspacing = -.1)
    setp(leg.get_texts(), fontsize = "small")
    axis(ymax = 2*6.125/8.125*81/79)
    yticks(arange(0,1.1,.5))

def plot_MP_fig2b(gammarange, bo06, bo08, booptimal, bo6h, basebalancing, figno = 1):
    matplotlib.rcParams['font.size'] = 10
    fig = figure(figno)
    clf()
    fig.set_size_inches([3.425, 6.125*3.425/8.125])
    ax1 = axes()
    ax1.set_position([.17,.15,.79,.81])
    plot(gammarange, bo06 - basebalancing, 'b', lw = 1, alpha = 1, label = r"$\alpha_W = 0.60$")
    plot(gammarange[:1], (bo08 - basebalancing)[:1], 'r', lw = 1, alpha = .75, label = r"$\alpha_W = 0.80$")
    plot(gammarange, booptimal - basebalancing, 'k', lw = 1, color = ".5", label = r"Optimal $\alpha_W$")
    plot(gammarange, bo08 - basebalancing, 'r', lw = 1, alpha = .75)
    plot(gammarange, bo6h - basebalancing, 'b--', lw = 1, label = r"6 h storage")
    xlabel(r"Average RES power generation factor $\gamma$")
    ylabel(r"Excess balancing fraction $E_B^{\mathrm{excess}}$")
#    setp(ax2, yscale = "log")
    axis([0,2,0,.22])
    leg = legend(loc = 1, labelspacing = -.1, bbox_to_anchor = (.07,0.06,1.005,1.018))
    setp(leg.get_texts(), fontsize = "small")

def plot_MP_fig11b(gammarange, bo06, bo6h, basebalancing, figno = 1):
    matplotlib.rcParams['font.size'] = 10
    fig = figure(figno)
    clf()
    fig.set_size_inches([3.425, 6.125*3.425/8.125])
    ax1 = axes()
    ax1.set_position([.17,.15,.79,.81])
#    ax2 = axes([0.63,0.5,0.33,0.33])
#    axes(ax1)
    plot(gammarange, bo06 - basebalancing, 'k', alpha = 1, lw = 1, label = r"no storage")
    plot(gammarange, bo6h - basebalancing, 'k', lw = 1, color = ".5", label = r"6h storage")
#    plot(gammarange, bo06, 'k', alpha = 1, lw = 1)
#    plot(gammarange, basebalancing, 'k:', lw = 1)
    xlabel(r"Renewable power generation $\gamma$")
    ylabel("Excess balancing")
    leg = legend(loc = 1, labelspacing = -.1, bbox_to_anchor = (.07,.06,1,1))
    setp(leg.get_texts(), fontsize = "small")
#    axis(ymax = 2*6.125/8.125)
#    yticks(arange(0,1.1,.5))
#    axes(ax2)
#    plot(gammarange, bo06 - bo6h, 'k', lw = 1, label = r"optimal $\alpha_W$")
#    plot(gammarange, bo6h - basebalancing, lw = 1, color = ".5")
#    yticks(arange(0,.11,.05))
#    xlabel(r"$\gamma$")
#    ylabel(r"Balancing difference")
#    setp(ax2, yscale = "log")
#    axis("tight")
    axis([0,2,0,.22])

def plot_MP_fig11c(gammarange, p2, p3, p4, figno = 1):
    matplotlib.rcParams['font.size'] = 10
    fig = figure(figno)
    clf()
    fig.set_size_inches([3.425, 6.125*3.425/8.125])
    ax1 = axes()
    ax1.set_position([.17,.15,.79,.81])
    plot(gammarange, p2[:,2], 'k', alpha = 1, lw = 1, label = r"Policy 1")
    plot(gammarange, p3[:,2], 'k', lw = 1, alpha = .5, label = r"Policy 2")
    plot(gammarange, p4[:,2], 'k', alpha = .25, lw = 1, label = "Policy 3")
    xlabel(r"Average RES power generation factor $\gamma$")
    ylabel("Balancing quantiles ($Q_{0.999}$) [av.h.l.]      ")
    leg = legend(loc = 0, labelspacing = -.1)
    setp(leg.get_texts(), fontsize = "small")
    axis(ymin = 0, ymax = 1.5)

def plot_MP_wind_solar_gen_fig(gammarange, alphas, figno = 4):
    fig = figure(figno)
    clf()
    fig.set_size_inches([3.425, 6.125*3.425/8.125])
    ax1 = axes()
    ax1.set_position([.15,.15,.81,.79])
    plot(gammarange, gammarange * alphas, 'b', lw = 1, label = "wind")
    plot(gammarange, gammarange * (1 - alphas), 'y', lw = 1, label = "solar")
    xlabel(r"Total renewable power generation $\gamma$")
    ylabel(r"Specific power generation")
    leg = legend(loc = 2, labelspacing = -.1)
    setp(leg.get_texts(), fontsize = "small")
    axis("tight")
    axis(ymax = 1.2)

def plot_MP_fig4a(bX, bY, bcont, alphas, ulim, olim, sc, figno = 4):
    fig = figure(figno)
    clf()
    fig.set_size_inches([3.425,6.125*3.425/8.125])
    ax1 = axes()
    ax1.set_position([.16,.15,.72,.8])
    contourf(bX, bY + 1e-5, bcont, 300)
    clim([0,bcont.max()])
    contourf(bX, bY + 1e-5, bcont, 300)
    clim([0,bcont.max()])
    contourf(bX, bY + 1e-5, bcont, 300)
    clim([0,bcont.max()])
    plot(alphas, linspace(1e-5, sc, alphas.size), 'k', lw = 1)
    plot(ulim/100., linspace(1e-5, sc, alphas.size), 'k--', lw = 1)
    plot(olim/100., linspace(1e-5, sc, alphas.size), 'k--', lw = 1)
    axis("tight")
    axis(ymin = .5)
    setp(ax1, yscale = 'log')
    xticks(arange(0,1.1,.2))
    cb = colorbar()
    cb.set_label(r"Balancing energy $E_B$")
    xlabel(r"Wind fraction $\alpha_W$")
    ylabel(r"Storage energy capacity [av.h.l.]")

def plot_MP_fig4b(bX, bY, bcont, alphas, ulim, olim, sc, g075bcont, figno = 4):
    fig = figure(figno)
    clf()
    fig.set_size_inches([3.425,6.125*3.425/8.125])
    ax1 = axes()
    ax1.set_position([.16,.15,.72,.8])
    contourf(bX, bY + 1e-5, bcont, 300)
    clim([0,g075bcont.max()])
    contourf(bX, bY + 1e-5, bcont, 300)
    clim([0,g075bcont.max()])
    contourf(bX, bY + 1e-5, bcont, 300)
    clim([0,g075bcont.max()])
    plot(alphas, linspace(1e-5, sc, alphas.size), 'k', lw = 1)
    plot(ulim/100., linspace(1e-5, sc, alphas.size), 'k--', lw = 1)
    plot(olim/100., linspace(1e-5, sc, alphas.size), 'k--', lw = 1)
    axis("tight")
    axis(ymin = .5)
    setp(ax1, yscale = 'log')
    xticks(arange(0,1.1,.2))
    cb = colorbar()
    cb.set_label(r"Balancing energy $E_B$")
    xlabel(r"Wind fraction $\alpha_W$")
    ylabel(r"Storage energy capacity [av.h.l.]")

def plot_MP_fig5a(g1sr, g1br, g2sr, g2br, g3sr, g3br, g1ob, g2ob, g3ob, dayscale, figno = 5):
    fig = figure(figno)
    clf()
    fig.set_size_inches([3.425,6.125*3.425/8.125])
    ax1 = axes()
    ax1.set_position([.17, .15, .79, .81])
    plot(g1sr, g1br/70128 - .25, 'k', lw = 1, label = "$\gamma = 0.75$")
    plot(g2sr, g2br/70128, 'k:', lw = 1, label = "$\gamma = 1.00$")
    plot(g3sr, g3br/70128, 'k--', lw = 1, label = "$\gamma = 1.25$")
    axis("tight")
    axis(xmax = 24, ymax = .22, ymin = 0)
    xticks(arange(0,4.1)*6)
    plot(dayscale, g1ob/70128 - .25, 'k', lw = 1, alpha = .25)
    plot(dayscale, g2ob/70128, 'k', lw = 1, alpha = .25)
    plot(dayscale, g3ob/70128, 'k', lw = 1, alpha = .25)
    xlabel(r"Storage energy capacity $C_S$ [av.h.l.]")
    ylabel(r"Excess balancing fraction $E_B^{\mathrm{excess}}$")
#    axes(ax1)
#def plot_MP_fig5b(g1sr,g1br,g2sr,g2br,g3sr,g3br, figno = 5):
#    fig = figure(figno)
#    clf()
#    fig.set_size_inches([3.425,6.125*3.425/8.125])
    axvline(x=6, ls='--', color = 'k')
    ax3 = axes([.5, .53, .8-(.6-.22), .8 - (.6 - .22)])
    plot(g1sr, g1br/70128 - .25, 'k', lw = 1, label = "$\gamma = 0.75$")
    plot(g2sr, g2br/70128, 'k:', lw = 1, label = "$\gamma = 1.00$")
    plot(g3sr, g3br/70128, 'k--', lw = 1, label = "$\gamma = 1.25$")
    yticks(arange(-0.1,.31,.05))
#    xticks(arange(0,.091,.045))
    setp(ax3, xscale = "log")
    axis(ymin = 0, xmin = 1., xmax = 1250)
    leg = legend(loc = 1, labelspacing = -.1, bbox_to_anchor = (.07,0.06,1.15,1.1))
    setp(leg.get_texts(), fontsize = "small")
    axvline(x=6, ls='--', color = 'k')


def plot_MP_storage_usage(p2,mm, tekst = "Policy 1", figno = 6):
    fig = figure(figno)
    clf()
    p2bund = (mm-p2).cumsum().min()
    fig.set_size_inches([3.425,(6.125*3.425/8.125)*2/3])
    ax1 = axes([.17,((1-.14*2/3)-(.15*2/3))/2+(.15*2/3)+.0333333,.79,((1-.14*2/3)-(.15*2/3))/2])
    fill_between(linspace(0,70128,70128)[1680:1802],(mm-p2).cumsum()[1680:1802]-p2bund,0, color = "g", alpha = .7)
    plot(linspace(0,70128,2),array([0,0]), 'r', lw = 1)
    plot(linspace(0,70128,2),array([6,6]), 'b', lw = 1)
    plot((mm-p2).cumsum()-p2bund, 'g', lw = 1)
    xticks(arange(1680,1800,24), ['','','','',''])
    yticks(arange(0,7,2))
    ylabel("[av.h.l.]")
    text((1680+1800)/2,5,tekst, ha = "center", va="center")
    axis([1680, 1800, -1, 7])
    ax2 = axes([.17,.15*2/3,.79,((1-.14*2/3)-(.15*2/3))/2])
    fill_between(linspace(0,70128,70128)[1680:1801], p2[1680:1801], mm[1680:1801], p2[1680:1801] <= mm[1680:1801], color = "green", alpha = .7)
    fill_between(linspace(0,70128,70128)[1680:1801], p2[1680:1801], mm[1680:1801], p2[1680:1801] >= mm[1680:1801], color = "green", alpha = .3)
    fill_between(linspace(0,70128,70128)[1680:1801], p2[1680:1801], 0, p2[1680:1801] <= 0, color = "red")
    fill_between(linspace(0,70128,70128)[1680:1801], p2[1680:1801], 0, p2[1680:1801] >= 0, color = "blue")
    plot(p2, 'k', lw = 1)
    plot(mm, 'k', lw = 1)
    plot(linspace(0,70128, 2), array([0,0]), 'k', lw = 1)
    xticks(arange(1680,1800,24), ['','','','',''])
    ax1.set_ylabel("[av.h.l.]", labelpad=18)
    axis([1680, 1800, -.9, 1.1])

#    ax2 = axes([.17,.375,.825,.30])
#    p3bund = (mm-p3).cumsum().min()
#    fill_between(linspace(0,70128,70128)[1680:1802],(mm-p3).cumsum()[1680:1802]-p3bund,0, color = "g", alpha = .7)
#    plot(linspace(0,70128,2),array([0,0]), 'r', lw = 1)
#    plot(linspace(0,70128,2),array([6,6]), 'b', lw = 1)
#    plot((mm-p3).cumsum()-p3bund, 'g', lw = 1)
#    xticks(arange(1680,1800,24), ['','','','',''])
#    yticks(arange(0,7,2))
#    ylabel("[av.h.l.]")
#    axis([1680, 1800, -1, 7])
#    ax3 = axes([.17,.065,.825,.30])
#    p4bund = (mm-p4).cumsum().min()
#    fill_between(linspace(0,70128,70128)[1680:1802],(mm-p4).cumsum()[1680:1802]-p4bund,0, color = "g", alpha = .7)
#    plot(linspace(0,70128,2),array([0,0]), 'r', lw = 1)
#    plot(linspace(0,70128,2),array([6,6]), 'b', lw = 1)
#    plot((mm-p4).cumsum()-p4bund, 'g', lw = 1)
    xticks(arange(1680,1800,24), ['            Fri','            Sat','           Sun','            Mon','            Tue'],rotation = 0)
    yticks(arange(-1,1.2,.5))
    ax2.set_ylabel("[av.h.l.]", labelpad = 1)
    axis([1680, 1800, -.9, 1.1])
    
def plot_MP_fig6(p2,p3,p4,mm, figno = 6):
    fig = figure(figno)
    clf()
    fig.set_size_inches([3.425,6.125*3.425/8.125])
    ax1 = axes([.17,.685,.825,.30])
    fill_between(linspace(0,70128,70128)[1680:1801], p2[1680:1801], mm[1680:1801], p2[1680:1801] <= mm[1680:1801], color = "green", alpha = .7)
    fill_between(linspace(0,70128,70128)[1680:1801], p2[1680:1801], mm[1680:1801], p2[1680:1801] >= mm[1680:1801], color = "green", alpha = .3)
    fill_between(linspace(0,70128,70128)[1680:1801], p2[1680:1801], 0, p2[1680:1801] <= 0, color = "red")
    fill_between(linspace(0,70128,70128)[1680:1801], p2[1680:1801], 0, p2[1680:1801] >= 0, color = "blue")
    plot(p2, 'k', lw = 1)
    plot(mm, 'k', lw = 1)
    plot(linspace(0,70128, 2), array([0,0]), 'k', lw = 1)
    xticks(arange(1680,1800,24), ['','','','',''])
    ylabel("[av.h.l.]")
    axis([1680, 1800, -.9, 1.1])
    ax2 = axes([.17,.375,.825,.30])
    fill_between(linspace(0,70128,70128)[1680:1801], p3[1680:1801], mm[1680:1801], p3[1680:1801] <= mm[1680:1801], color = "green", alpha = .7)
    fill_between(linspace(0,70128,70128)[1680:1801], p3[1680:1801], mm[1680:1801], p3[1680:1801] >= mm[1680:1801], color = "green", alpha = .3)
    fill_between(linspace(0,70128,70128)[1680:1801], p3[1680:1801], 0, p3[1680:1801] <= 0, color = "red")
    fill_between(linspace(0,70128,70128)[1680:1801], p3[1680:1801], 0, p3[1680:1801] >= 0, color = "blue")
    plot(p3, 'k', lw = 1)
    plot(mm, 'k', lw = 1)
    plot(linspace(0,70128, 2), array([0,0]), 'k', lw = 1)
    xticks(arange(1680,1800,24), ['','','','',''])
    ylabel("[av.h.l.]")
    axis([1680, 1800, -.9, 1.1])
    ax3 = axes([.17,.065,.825,.30])
    fill_between(linspace(0,70128,70128)[1680:1801], p4[1680:1801], mm[1680:1801], p4[1680:1801] <= mm[1680:1801], color = "green", alpha = .7)
    fill_between(linspace(0,70128,70128)[1680:1801], p4[1680:1801], mm[1680:1801], p4[1680:1801] >= mm[1680:1801], color = "green", alpha = .3)
    fill_between(linspace(0,70128,70128)[1680:1801], p4[1680:1801], 0, p4[1680:1801] <= 0, color = "red")
    fill_between(linspace(0,70128,70128)[1680:1801], p4[1680:1801], 0, p4[1680:1801] >= 0, color = "blue")
    plot(p4, 'k', lw = 1)
    plot(mm, 'k', lw = 1)
    plot(linspace(0,70128, 2), array([0,0]), 'k', lw = 1)
    xticks(arange(1680,1800,24), ['            Fri','            Sat','           Sun','            Mon','            Tue'],rotation = 0)
    ylabel("[av.h.l.]")
    axis([1680, 1800, -.9, 1.1])
    
def plot_MP_sufficience_contours(gammas, alphas, balafsec, bndbalafsec, bal6h, percent = 10., ballevel = (150000./(370*8766)), figno = 14, tekst = "With 6-hour storage"):
    fr = percent/100.
    fig = figure(figno)
    clf()
    fig.set_size_inches([3.425,6.125*3.425/8.125])
    ax1 = axes([.17,.15,.79,.81])
    contourf(gammas, alphas, balafsec.transpose()/70128, [0,ballevel], colors = [(.8,0,0,.25)], alpha = .25)
    contourf(gammas, alphas, bndbalafsec.transpose()/70128, [0,ballevel], colors = [(0,.8,0,.25)], alpha = .25)
    contourf(gammas, alphas, bal6h.transpose()/70128, [0,ballevel], colors = [(0,0,.8,.25)], alpha = .25)
    contour(gammas, alphas, balafsec.transpose()/70128, array([1.-fr,1.,1+fr])*ballevel, linestyles = ['--', '-', '-.'], colors = 'r', linewidths = 1)
    contour(gammas, alphas, bndbalafsec.transpose()/70128, array([1.-fr,1.,1+fr])*ballevel, linestyles = ['--', '-', '-.'], colors = ['g','g','g'], linewidths = 1)
    contour(gammas, alphas, bal6h.transpose()/70128, array([1.-fr,1.,1+fr])*ballevel, linestyles = ['--', '-', '-.'], colors = 'b', linewidths = 1)
    xlabel(r"Average RES power generation factor $\gamma$")
    ylabel(r"Wind fraction $\alpha_W$")
    text(1.02,.9,tekst, horizontalalignment="left", weight = "bold")

#Legend
    pp_unl = Line2D((0,0),[1],[1],color='r',lw=1)
    pp_lim = Line2D((0,0), [1],[1], color = 'g', lw = 1)
    pp_no = Line2D((0,0), [1], [1], color = 'b', lw = 1)

    pp = (pp_unl, pp_lim, pp_no)
    pp_text = (r'Unlimited H$_\mathsf{2}$ storage', 'Limited H$_2$ storage', 'No H$_2$ storage')
    leg = legend(pp,pp_text,loc = 'center right', labelspacing = -.1, bbox_to_anchor = (1, .6))#, bbox_transform=gcf().transFigure,ncol=4,columnspacing=.3,handlelength=2,handletextpad=0.3,borderaxespad=0.0)
    ltext  = leg.get_texts();
    setp(ltext, fontsize='small')    # the legend text fontsize


    legend()
    axis(xmin = 1)

def plot_SLIDES_sufficience_contours(gammas, alphas, bndbalafsec, bndbalafsec6h, ballevel = (150000./(370*8766)), figno = 14, tekst = "With 6-hour storage"):
    fig = figure(figno)
    clf()
    fig.set_size_inches([3.425,6.125*3.425/8.125])
    ax1 = axes([.17,.17,.78,.79])
    contourf(gammas, alphas, bndbalafsec6h.transpose()/70128, [0,ballevel], colors = [(.75,.75,1)])    
    contourf(gammas, alphas, bndbalafsec.transpose()/70128, [0,ballevel], colors = [(.75,1,.75)])
    contour(gammas, alphas, bndbalafsec.transpose()/70128, [ballevel], linestyle = '-',  linewidths = 1, colors = 'g')
    contour(gammas, alphas, bndbalafsec6h.transpose()/70128, [ballevel], linestyle = '-', linewidths = 1, colors = 'b')
    xlabel(r"Average RES power generation $\gamma$")
    ylabel(r"Wind fraction $\alpha_W$")
    text(1.4,.94,tekst, horizontalalignment="left", weight = "ultralight", size = "xx-small")

#Legend
    pp_lim = Line2D((0,0), [1],[1], color = 'b', lw = 1)
    pp_no = Line2D((0,0), [1], [1], color = 'g', lw = 1)

    pp = (pp_lim, pp_no)
    pp_text = ('With 6-hour storage', 'Without 6-hour storage')
    leg = legend(pp,pp_text,loc = 'lower left', labelspacing = -.1)#, bbox_to_anchor = (1, .6))#, bbox_transform=gcf().transFigure,ncol=4,columnspacing=.3,handlelength=2,handletextpad=0.3,borderaxespad=0.0)
    ltext  = leg.get_texts();
    setp(ltext, fontsize='small')    # the legend text fontsize


    legend()
    axis(xmin = 1)

def plot_MP_storage_contours(gammas, alphas, secsto, stolim = 25000/(370.*8766), percent = 10., figno = 15):
    fr = percent/100.
    fig = figure(figno)
    clf()
    fig.set_size_inches([3.425,6.125*3.425/8.125])
    ax1 = axes([.145,.15,.74, .82])
    contourf(gammas, alphas, secsto.transpose()/8766,300)
    contourf(gammas, alphas, secsto.transpose()/8766,300)
    contourf(gammas, alphas, secsto.transpose()/8766,300)
    cb = colorbar(ticks = arange(0,.25,.025), format = '%.3f')
    xlabel(r"Average RES power generation $\gamma$")
    ylabel(r"Wind fraction $\alpha_W$")
    contour(gammas, alphas, secsto.transpose()/8766, array([1.-fr,1.,1.+fr]) * stolim, colors = "1", linestyles = ['--','-','-.'], linewidths = 1)
    yticks(arange(0,1.1,.2))
    contourf(gammas, alphas, secsto.transpose()/8766,300)
    axis([0,2,0,1])
    ax2 = cb.ax
    ax2.yaxis.set_label_coords(5.8,.485)
    cb.set_label("Storage energy capacity [av.y.l.]")

def plot_MP_fig13(gammarange, sc, no6hsc, limit = 25000./(370*8766), figno = 13):
    fig = figure(figno)
    clf()
    fig.set_size_inches([3.425,6.125*3.425/8.125])
    ax1 = axes([.17,.15,.79,.81])
#    print gammarange.size, sc.size
    plot(gammarange, sc/8766, 'k', lw = 1, label = "With 6 h storage")
    plot(gammarange, limit * ones(gammarange.size), 'k', alpha = .75, lw = 1)
    plot(gammarange, no6hsc/8766, 'k', alpha = .25, lw = 1, label = "No 6 h storage")
    axvline(x=1, ls='--', color = 'k')
    xlabel(r"Average RES power generation factor $\gamma$")
    ax1.set_ylabel(r"Storage energy capacity $C_S$ [av.y.l.]", y = .45)
    leg = legend(loc = 2, labelspacing = -.1)
    axis(ymax = .135)
    setp(leg.get_texts(), fontsize = "small")
    text(.2,limit+.002,"25 TWh", rotation = 0, size = "small")


def plot_MP_fig17(gammarange, fullsto, bndsto, only6h, hydropotential = 150000./(150*8766), figno = 17):
    matplotlib.rcParams['font.size'] = 10
    fig = figure(figno)
    clf()
    fig.set_size_inches([3.425, 6.125*3.425/8.125])
    ax1 = axes()
    ax1.set_position([.17,.15,.79,.81])
    ax2 = axes([0.35,0.61,0.58,0.33])
    axes(ax1)
    plot(gammarange, fullsto/70128, 'k:', alpha = 1, lw = 1, label = r"Unlimited H$_2$ storage")
    plot(gammarange, bndsto/70128, 'k', alpha = 1, lw = 1, label = r"Limited H$_2$ storage")
    plot(gammarange, only6h/70128, 'k--', lw = 1, label = r"Only 6 h storage")
    plot(gammarange, hydropotential*ones(101), 'k', alpha = 1, lw = 1)
    xlabel(r"Average RES power generation factor $\gamma$")
    ylabel("Balancing energy fraction $E_B$")
    leg = legend(loc = 7, labelspacing = -.1, bbox_to_anchor = (0.05,-.25,1,1))
    setp(leg.get_texts(), fontsize = "small")
    axis(ymax = 2*6.125/8.125*81/79)
    yticks(arange(0,1.6,.25))
    axes(ax2)
    plot(gammarange, fullsto/70128, 'k:', lw = 1)
    plot(gammarange, bndsto/70128, 'k', lw = 1, alpha = 1)
    plot(gammarange, only6h/70128, 'k--', lw = 1)
    plot(gammarange, hydropotential*ones(101), 'k', lw = 1)
    xlabel(r"$\gamma$")
    ax2.xaxis.set_label_coords(.5,-.1)
    ylabel(r"$E_B$")
#    setp(ax2, yscale = "log")
    axis([1.01, 1.135, 0, .1])
    yticks(arange(0,.11,.03))
    yplaces = [-interp(-hydropotential,-fullsto/70128,-gammarange),-interp(-hydropotential,-bndsto/70128,-gammarange),-interp(-hydropotential,-only6h/70128,-gammarange)]
    xticks(yplaces, ('%.3f'%yplaces[0],'%.3f'%yplaces[1],'%.3f'%yplaces[2]), rotation = 45)
    print yplaces

def get_afvigelse(x,y,sr,br,z0 = None, z1 = None):
    if z0 is None:
        z0 = br[-1]
    if z1 is None:
        z1 = br[0]
    bid1 = interp(sr[sr<x],[0,x],[z1,y])
    bid2 = interp(sr[sr>x],[x,sr[-1]],[y,z0])
    #plot(sr,br)
    #plot(sr[sr<=x],bid1)
    #plot(sr[sr>x],bid2)
    samlet = concatenate([bid1,[y],bid2])
    xsr = concatenate([sr[sr<x],[x],sr[sr>x]])
    xbr = concatenate([br[sr<x],[interp(x,sr,br)],br[sr>x]])
    #plot(xsr,xbr)
    #plot(xsr,samlet)
    #plot(xsr,square(samlet-xbr))
    return trapz(square(samlet-xbr),xsr)
    
def find_min_afvigelse(sr,br):
    return fminbound(least_y_afvigelse,sr.min(),sr.max(),(sr,br))

def find_least_y_afvigelse(x,sr,br):
    def get_y_afvigelse(y):
        return find_mindste_afvigelse(x,y,sr,br)
    return fminbound(get_y_afvigelse,br.min(),br.max())

def least_y_afvigelse(x,sr,br):
    return find_mindste_afvigelse(x,find_least_y_afvigelse(x,sr,br),sr,br)

def find_min_constrained_afvigelse(sr,br):
    def tmp(x):
        return get_afvigelse(x,interp(x,sr,br),sr,br)
    return fminbound(tmp,sr.min(),sr.max())

def find_optimal_z0(x,y,sr,br):
    def tmp(z):
        return get_afvigelse(x,y,sr,br,z)
    return None# fminbound(tmp,-square(br.max()),square(br.max()))

def find_bedste_start(xr,yr):
    def pris(y):
        return trapz(square(interp(xr,[xr[0],xr[-1]],[y,yr[-1]])-yr),xr)
    return fminbound(pris,yr.min()-2*(yr.max()-yr.min()),yr.max()+2*(yr.max()-yr.min()))

def find_bedste_slut(xr,yr):
    def pris(y):
        return trapz(square(interp(xr,[xr[0],xr[-1]], [yr[0], y])-yr),xr)
    return fminbound(pris,yr.min()-2*(yr.max()-yr.min()),yr.max()+2*(yr.max()-yr.min()))

def find_mindste_afvigelse(x,y,sr,br):
    xbid1 = concatenate([sr[sr<x],[x]])
    xbid2 = concatenate([[x],sr[sr>x]])
    ybid1 = interp(xbid1,sr,br)
    ybid2 = interp(xbid2,sr,br)
    y0 = find_bedste_start(xbid1,ybid1)
    y1 = find_bedste_slut(xbid2,ybid2)
    xr = concatenate([xbid1[:-1],xbid2])
    return trapz(square(interp(xr, [sr[0],x,sr[-1]],[y0,y,y1])-interp(xr,sr,br)),xr)
    
