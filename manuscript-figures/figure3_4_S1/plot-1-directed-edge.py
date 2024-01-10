# A network with two nodes and one directed edge.
# Generate figures 3 and S1(a)-(c) in the manuscript.
# figure 3 if if_robustness == 0
# figure S1(a)-(c) if if_robustness == 1
#
# The output figure is saved in fig.pdf

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy import linspace

type=2
# (a) type==0
# (b) type==1
# (c) type==2

# whether or not to run the robustness test with respect to the two $r$ values
if_robustness = 1

def dist_2_gaussian(mu1, std1, mu2, std2):

    # symmetrized KL divergence
#    return ( (std1**2 - std2**2)**2 + (std1**2 + std2**2)*(mu1-mu2)**2 ) / (2 * std1**2 * std2**2)

    # t-test based index
    return ( np.absolute(mu1-mu2) / np.sqrt(std1**2 + std2**2) )

def signal_quality(r, Deltar, w, D1, D2):
    x1 = -np.sqrt(-r)
    x2 = -np.sqrt(-r + Deltar + w * (np.sqrt(-r) - 1))
    mean_signal = np.empty(3)
    std_signal = np.empty(3)
    # {mean,std}_signal[0]: var(x1)
    # {mean,std}_signal[1]: var(x2)
    # {mean,std}_signal[2]: (var(x1) + var(x2))/2
    C11 = - D1**2 / (4*x1) # expected variance of x1
    C12 = w * D1**2 / (8*x1*(x1+x2)) # expected covariance between x1 and x2
    C22 = -w**2 * D1**2 / (16*x1*x2*(x1+x2)) - D2**2 / (4*x2) # expected variance of x2
    mean_signal[0] = C11    
    mean_signal[1] = C22
    mean_signal[2] = (C11 + C22) / 2
    std_signal[0] = np.sqrt(2/(L-1)) * C11
    std_signal[1] = np.sqrt(2/(L-1)) * C22
    std_signal[2] = np.sqrt( (C11**2 + 2 * C12**2 + C22**2) / (2*(L-1)) )
    return mean_signal, std_signal

Nsamples = 100
r_min = -0.5 # -1
r_max = -0.008
rvalues = np.linspace(r_min, r_max, Nsamples)
D1 = 0.1

if type==0:
    Deltar = 1.0
    D2 = 0.1
elif type==1:
    Deltar = 0.5
    D2 = 0.1
elif type==2:
    Deltar = 1 
    D2 = 0.2
else:
    raise ValueError('0 <= type <= 2 violated')

L = 100 # number of samples

w = 0.5 # Coupling strength. w < Deltar must be satisfied
mean_signal_across_r = np.empty((Nsamples, 3))
std_signal_across_r = np.empty((Nsamples, 3))
for i in range(Nsamples):
    r = rvalues[i]
    mean_signal, std_signal = signal_quality(r, Deltar, w, D1, D2)
    mean_signal_across_r[i, :] = mean_signal
    std_signal_across_r[i, :] = std_signal

# necessary to get the fonts right
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"]})

if if_robustness == 0: # no robustness test

    # plot \hat{V}_1
    plt.plot(rvalues, mean_signal_across_r[:,0], 'm-', linewidth=1.5, markerfacecolor='none', label=r'$\hat{V}_1$')
    #plt.plot(rvalues, mean_signal_across_r[:,0] + std_signal_across_r[:,0], 'm-', linewidth=0.5, markerfacecolor='none')
    #plt.plot(rvalues, mean_signal_across_r[:,0] - std_signal_across_r[:,0], 'm-', linewidth=0.5, markerfacecolor='none')
    plt.fill_between(rvalues, mean_signal_across_r[:,0] - std_signal_across_r[:,0], mean_signal_across_r[:,0] + std_signal_across_r[:,0], color='magenta', alpha=.3) # transparency option only works for pdf, not eps

    # plot \hat{V}_2
    plt.plot(rvalues, mean_signal_across_r[:,1], 'r-', linewidth=2, markerfacecolor='none', label=r'$\hat{V}_2$')
    #plt.plot(rvalues, mean_signal_across_r[:,1] + std_signal_across_r[:,1], 'r-', linewidth=0.5, markerfacecolor='none')
    #plt.plot(rvalues, mean_signal_across_r[:,1] - std_signal_across_r[:,1], 'r-', linewidth=0.5, markerfacecolor='none')
    plt.fill_between(rvalues, mean_signal_across_r[:,1] - std_signal_across_r[:,1], mean_signal_across_r[:,1] + std_signal_across_r[:,1], color='red', alpha=.3)

    # plot \hat{V}_{1, 2}
    plt.plot(rvalues, mean_signal_across_r[:,2], 'b-', linewidth=2, markerfacecolor='none', label=r'$\hat{V}_{\{1,2\}}$')
    #plt.plot(rvalues, mean_signal_across_r[:,2] + std_signal_across_r[:,2], 'b-', linewidth=0.5, markerfacecolor='none')
    #plt.plot(rvalues, mean_signal_across_r[:,2] - std_signal_across_r[:,2], 'b-', linewidth=0.5, markerfacecolor='none')
    plt.fill_between(rvalues, mean_signal_across_r[:,2] - std_signal_across_r[:,2], mean_signal_across_r[:,2] + std_signal_across_r[:,2], color='blue', alpha=.3)

    # calculate signal quality index
    rcheck = np.empty(2)
    rcheck[0] = -0.3
    rcheck[1] = -0.1
    print('checkpoint r =', rcheck[0], ',', rcheck[1])
    mean_signal_r0, std_signal_r0 = signal_quality(rcheck[0], Deltar, w, D1, D2)
    mean_signal_r1, std_signal_r1 = signal_quality(rcheck[1], Deltar, w, D1, D2)
    q = np.empty(3)
    for i in range(3):
        q[i] = dist_2_gaussian(mean_signal_r0[i], std_signal_r0[i], mean_signal_r1[i], std_signal_r1[i])
    print('quality index =', q)

    plt.xlabel(r'$r$', fontsize=24)
    #if data >= 1: # Fig. 2
    plt.ylabel('sample covariance', fontsize=24)
    plt.xlim(r_min, 0)
    plt.ylim(0, 0.03)
    plt.xticks(np.linspace(-0.5, 0, 6, endpoint=True), ('\N{MINUS SIGN}0.5', '\N{MINUS SIGN}0.4', '\N{MINUS SIGN}0.3', '\N{MINUS SIGN}0.2', '\N{MINUS SIGN}0.1', '0'), fontsize=20,fontname='Helvetica')
    #    plt.ylim(0, 1.01)
    plt.yticks(np.linspace(0, 0.03, 4, endpoint=True), ('0', '0.01', '0.02', '0.03'), fontsize=20)

    if type==0:
        plt.title('(a)', fontsize=28, x=-0.08, y=1.05)
        plt.legend(loc = 'upper left', numpoints = 1, labelspacing=0.25, frameon=False, fontsize=20)
    elif type==1:
        plt.title('(b)', fontsize=28, x=-0.08, y=1.05)
    elif type==2:
        plt.title('(c)', fontsize=28, x=-0.08, y=1.05)

    plt.subplots_adjust(bottom=0.15)
    plt.subplots_adjust(left=0.16)

else: # if_robustness == 1, so carry out the robustness test

    which_is_best = np.empty((Nsamples, Nsamples))
    q = np.empty(3)
    for i in range(Nsamples): # for all r values
        which_is_best[i, i] = -1
        for j in range(i):
            for k in range(3):
                q[k] = dist_2_gaussian(mean_signal_across_r[i, k], std_signal_across_r[i, k], mean_signal_across_r[j, k], std_signal_across_r[j, k])
            which_is_best[i, j] = np.argmax(q)
            which_is_best[j, i] = which_is_best[i, j]

#    print(which_is_best)            
    plt.xlabel(r'$r$', fontsize=24)
    #if data >= 1: # Fig. 2
    plt.ylabel(r'$r$', fontsize=24)
    plt.xticks(np.linspace(-0.5, 0, 6, endpoint=True), ('\N{MINUS SIGN}0.5', '\N{MINUS SIGN}0.4', '\N{MINUS SIGN}0.3', '\N{MINUS SIGN}0.2', '\N{MINUS SIGN}0.1', '0'), fontsize=20,fontname='Helvetica')
    plt.yticks(np.linspace(-0.5, 0, 6, endpoint=True), ('\N{MINUS SIGN}0.5', '\N{MINUS SIGN}0.4', '\N{MINUS SIGN}0.3', '\N{MINUS SIGN}0.2', '\N{MINUS SIGN}0.1', '0'), fontsize=20,fontname='Helvetica')
    plt.xlim(r_min, r_max)
    plt.ylim(r_min, r_max)

    cmap = matplotlib.colors.ListedColormap(['white','plum','tomato','dodgerblue','cyan', 'orange'])
    bounds = np.array([-1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5])
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
    plt.imshow(which_is_best, interpolation='nearest', origin='lower', extent=[r_min, r_max, r_min, r_max], aspect=1.0, cmap=cmap, norm=norm) # no smoothing. For a technical reason of Python, interpolation='none' produces a smoothed fig when the figure is saved as .eps
#    plt.colorbar()
#    plt.subplots_adjust(right=0.95)

#        im = ax.imshow(nR.T, interpolation='nearest', origin='lower', extent=[r_min, -0.008, r_min, -0.008nf_rate_max, 0.0, b_max], cmap='jet') # no smoothing. For a technical reason of Python, interpolation='none' produces a smoothed fig when the figure is saved as .eps

    if type==0:
        plt.title('(a)', fontsize=28, x=-0.08, y=1.05)
        plt.text(-0.4, -0.15, r'$\hat{V}_{\{1, 2 \}}$', fontsize=24)
    elif type==1:
        plt.title('(b)', fontsize=28, x=-0.08, y=1.05)
        plt.text(-0.4, -0.15, r'$\hat{V}_{\{1, 2 \}}$', fontsize=24)
        plt.text(-0.12, 0.02, r'$\hat{V}_{2}$', fontsize=24)
        plt.arrow(-0.08, 0.02, 0.05, -0.026, linewidth=0.7, color='black', head_width=0.007, clip_on=False)
    elif type==2:
        plt.title('(c)', fontsize=28, x=-0.08, y=1.05)
        plt.text(-0.4, -0.15, r'$\hat{V}_{1}$', fontsize=24)
        plt.text(-0.39, 0.02, r'$\hat{V}_{\{1, 2\}}$', fontsize=24)
        plt.arrow(-0.39, 0.02, -0.04, -0.027, linewidth=0.7, color='black', head_width=0.007, clip_on=False)

#plt.subplots_adjust(bottom=0.16, left=0.15, right=0.95)
#plt.tick_params(labelsize=20)

plt.savefig("fig.pdf", bbox_inches='tight')
# https://stackoverflow.com/questions/4042192/reduce-left-and-right-margins-in-matplotlib-plot