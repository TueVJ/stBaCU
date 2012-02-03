from pylab import *
import numpy
from scipy.optimize import *
from scipy.integrate import *
from scipy.stats.mstats import mquantiles
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

###
#	@return: indices: 	indices in ts where ts changes sign.
#			integrals:	The integral of ts in the same-sign intervals.
###

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
    
##	-> cooppolicies.py:
##	get_storage
##	get_storage_after_policy_4

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

##	-> cooppolicies.py:
##	get_policy_2_storage
##	get_policy_1_storage
##	get_policy_5_storage
##	get_policy_7_storage
##	get_policy_8_storage
##	get_policy_9_storage
##	old_storage
##

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
    
##
##	get_mismatch_after_double_storage
##	get_mismatch_after_policy_N block
##

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
    
