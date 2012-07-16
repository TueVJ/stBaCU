from pylab import NaN
from numpy import array
import numpy
#from scipy.optimize import *
#from scipy.integrate import *
#from scipy.stats.mstats import mquantiles
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import countries
import constant
from copy import deepcopy
from mfunc import *
from mpolicies import *
from mismatch import *
from mutils import pos,neg


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













