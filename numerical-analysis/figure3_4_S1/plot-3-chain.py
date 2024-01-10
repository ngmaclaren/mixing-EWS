# A chain with 3 nodes.
# Generate figures 4 and S1(d)-(f) in the manuscript.
# figure 4 if if_robustness == 0
# figure S1(d)-(f) if if_robustness == 1
#
# The output figure is saved in fig.pdf

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy import linspace

type=1
# (a) type==0
# (b) type==1
# (c) type==2

# whether or not to run the robustness test with respect to the two $r$ values
if_robustness = 1

def dist_2_gaussian(mu1, std1, mu2, std2):
    return ( np.absolute(mu1-mu2) / np.sqrt(std1**2 + std2**2) )

def signal_quality(r, w, D1, D2, x1, x2):
    tolerance = 1e-6
    error = 1.0
    while error > tolerance:
        x1_next = -np.sqrt(-r - w*(x2+1))
        x2_next = -np.sqrt(-r - 2*w*(x1+1))
        error = np.absolute(x1_next - x1) + np.absolute(x2_next - x2)
        x1 = x1_next
        x2 = x2_next
    # equilibrium x1 and x2 obtained
#    print(r, x1, x2)

    mean_signal = np.empty(5)
    std_signal = np.empty(5)
    # {mean,std}_signal[0]: var(x1) (= var(x3))
    # {mean,std}_signal[1]: var(x2)
    # {mean,std}_signal[2]: (var(x1) + var(x2))/2
    # {mean,std}_signal[3]: (var(x1) + var(x3))/2
    # {mean,std}_signal[4]: (var(x1) + var(x2) + var(x3))/3
    C11 = ( - 4*x1*x2*(x1+x2)*(D1**2) + (w**2) * ( (2*x1+x2)*(D1**2) - x1*(D2**2) ) ) / ( 8*x1*(2*x1*x2-(w**2))*(x1+x2) )
    C12 = w * (x2*(D1**2)+x1*(D2**2)) / ( 4*(2*x1*x2-(w**2))*(x1+x2) )
    C13 = - (w**2)*(x2*(D1**2) + x1*(D2**2)) / ( 8*x1*(2*x1*x2-(w**2))*(x1+x2) )
    C22 = ( -2*x1*x2*(x1+x2)*(D2**2) + (w**2)*x2*(-D1**2+D2**2) ) / ( 4*x2*(2*x1*x2-(w**2))*(x1+x2) )
    mean_signal[0] = C11 # expected variance of x1
    mean_signal[1] = C22
    mean_signal[2] = (C11 + C22) / 2
    mean_signal[3] = C11 # (x1+x3)/2
    mean_signal[4] = (2*C11 + C22) / 3
    std_signal[0] = np.sqrt(2/(L-1)) * C11
    std_signal[1] = np.sqrt(2/(L-1)) * C22
    std_signal[2] =np.sqrt( (C11**2 + 2 * C12**2 + C22**2) / (2*(L-1)) )
    std_signal[3] =np.sqrt( (2 * C11**2 + 2 * C13**2) / (2*(L-1)) )
    std_signal[4] = np.sqrt( 2 * (2 * C11**2 + 4 * C12**2 + 2 * C13**2 + C22**2) / (9*(L-1)) )
    return mean_signal, std_signal, x1, x2

Nsamples = 100
D2 = 0.1

if type==0:
    w = 0.05 # Coupling strength
    D1 = 0.1
elif type==1:
    w = 0.05
    D1 = 0.7
elif type==2:
    w = 0.05
    D1 = 0.015
else:
    raise ValueError('0 <= type <= 2 violated')

r_min = -0.5 # -1
r_max = -2*w* ( (w+1) - np.sqrt(w*(w+1)) )
print('r_max =', r_max)
rvalues = np.linspace(r_min, r_max - 0.01, Nsamples)

L = 100 # number of samples

mean_signal_across_r = np.empty((Nsamples, 5))
std_signal_across_r = np.empty((Nsamples, 5))
x1 = -1.0
x2 = -1.0
for i in range(Nsamples):
    r = rvalues[i]
    mean_signal, std_signal, x1, x2 = signal_quality(r, w, D1, D2, x1, x2)
    mean_signal_across_r[i, :] = mean_signal
    std_signal_across_r[i, :] = std_signal

# necessary to get the fonts right
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"]})

if if_robustness == 0: # no robustness test

    # x1 only
    # plt.plot(rvalues, mean_signal_across_r[:,0], 'k-', linewidth=1.5, markerfacecolor='none', label=r'$\hat{V}_1$')
    # plt.plot(rvalues, mean_signal_across_r[:,0] + std_signal_across_r[:,0], 'k-', linewidth=0.5, markerfacecolor='none')
    # plt.plot(rvalues, mean_signal_across_r[:,0] - std_signal_across_r[:,0], 'k-', linewidth=0.5, markerfacecolor='none')

    # x2 only
    plt.plot(rvalues, mean_signal_across_r[:,1], 'r-', linewidth=2, markerfacecolor='none', label=r'$\hat{V}_2$')
    #plt.plot(rvalues, mean_signal_across_r[:,1] + std_signal_across_r[:,1], 'k-', linewidth=0.5, markerfacecolor='none')
    #plt.plot(rvalues, mean_signal_across_r[:,1] - std_signal_across_r[:,1], 'k-', linewidth=0.5, markerfacecolor='none')
    plt.fill_between(rvalues, mean_signal_across_r[:,1] - std_signal_across_r[:,1], mean_signal_across_r[:,1] + std_signal_across_r[:,1], color='red', alpha=.3) # transparency option only works for pdf, not eps

    # x1 and x2
    plt.plot(rvalues, mean_signal_across_r[:,2], 'b-', linewidth=2, markerfacecolor='none', label=r'$\hat{V}_{\{1,2\}}$')
    #plt.plot(rvalues, mean_signal_across_r[:,2] + std_signal_across_r[:,2], 'b-', linewidth=0.5, markerfacecolor='none')
    #plt.plot(rvalues, mean_signal_across_r[:,2] - std_signal_across_r[:,2], 'b-', linewidth=0.5, markerfacecolor='none')
    plt.fill_between(rvalues, mean_signal_across_r[:,2] - std_signal_across_r[:,2], mean_signal_across_r[:,2] + std_signal_across_r[:,2], color='blue', alpha=.3) # transparency option only works for pdf, not eps

    # x1 and x3
    plt.plot(rvalues, mean_signal_across_r[:,3], 'c-', linewidth=2, markerfacecolor='none', label=r'$\hat{V}_{\{1,3\}}$')
    #plt.plot(rvalues, mean_signal_across_r[:,3] + std_signal_across_r[:,3], 'c-', linewidth=0.5, markerfacecolor='none')
    #plt.plot(rvalues, mean_signal_across_r[:,3] - std_signal_across_r[:,3], 'c-', linewidth=0.5, markerfacecolor='none')
    plt.fill_between(rvalues, mean_signal_across_r[:,3] - std_signal_across_r[:,3], mean_signal_across_r[:,3] + std_signal_across_r[:,3], color='cyan', alpha=.3) # transparency option only works for pdf, not eps

    # all 3 nodes
    plt.plot(rvalues, mean_signal_across_r[:,4], color='orange', linestyle='solid', linewidth=2, markerfacecolor='none', label=r'$\hat{V}_{\{1,2,3\}}$')
    #plt.plot(rvalues, mean_signal_across_r[:,4] + std_signal_across_r[:,4], 'b-', linewidth=0.5, markerfacecolor='none')
    #plt.plot(rvalues, mean_signal_across_r[:,4] - std_signal_across_r[:,4], 'b-', linewidth=0.5, markerfacecolor='none')
    plt.fill_between(rvalues, mean_signal_across_r[:,4] - std_signal_across_r[:,4], mean_signal_across_r[:,4] + std_signal_across_r[:,4], color='orange', alpha=.3) # transparency option only works for pdf, not eps

    # calculate signal quality index
    rcheck = np.empty(2)
    rcheck[0] = -0.3
    rcheck[1] = -0.1
    print('checkpoint r =', rcheck[0], ',', rcheck[1])
    mean_signal_r0, std_signal_r0, x1_tmp, x2_tmp = signal_quality(rcheck[0], w, D1, D2, -0.8, -0.8)
    mean_signal_r1, std_signal_r1, x1_tmp, x2_tmp = signal_quality(rcheck[1], w, D1, D2, -0.8, -0.8)
    q = np.empty(5)
    for i in range(5):
        q[i] = dist_2_gaussian(mean_signal_r0[i], std_signal_r0[i], mean_signal_r1[i], std_signal_r1[i])
    print('quality index =', q)

    plt.xlabel(r'$r$', fontsize=24)
    #if data >= 1: # Fig. 2
    plt.ylabel('sample covariance', fontsize=24)
    #plt.xlim(r_min, r_max)
    plt.xlim(-0.4, r_max-0.005)
    plt.xticks(np.linspace(-0.4, -0.1, 4, endpoint=True), ('\N{MINUS SIGN}0.4', '\N{MINUS SIGN}0.3', '\N{MINUS SIGN}0.2', '\N{MINUS SIGN}0.1'), fontsize=20)

    if type==0:
        plt.title('(a)', fontsize=28, x=-0.08, y=1.05)
        plt.legend(loc = 'upper left', numpoints = 1, labelspacing=0.25, frameon=False, fontsize=20)
        plt.yticks(np.linspace(0, 0.03, 4, endpoint=True), ('0', '0.01', '0.02', '0.03'), fontsize=20)
        plt.ylim(0, 0.026)
    elif type==1:
        plt.title('(b)', fontsize=28, x=-0.08, y=1.05)
        plt.yticks(np.linspace(0, 0.3, 7, endpoint=True), ('0', '0.05', '0.1', '0.15', '0.2', '0.25', '0.3'), fontsize=20)
        plt.ylim(0, 0.31)
    elif type==2:
        plt.title('(c)', fontsize=28, x=-0.08, y=1.05)
        plt.yticks(np.linspace(0, 0.03, 4, endpoint=True), ('0', '0.01', '0.02', '0.03'), fontsize=20)
        plt.ylim(0, 0.026)

else: # if_robustness == 1, so carry out the robustness test

    which_is_best = np.empty((Nsamples, Nsamples))
    q = np.empty(5)
    for i in range(Nsamples): # for all r values
        which_is_best[i, i] = -1
        for j in range(i):
            for k in range(5):
                q[k] = dist_2_gaussian(mean_signal_across_r[i, k], std_signal_across_r[i, k], mean_signal_across_r[j, k], std_signal_across_r[j, k])
            which_is_best[i, j] = np.argmax(q)
            which_is_best[j, i] = which_is_best[i, j]

#    print(which_is_best)            
    plt.xlabel(r'$r$', fontsize=24)
    #if data >= 1: # Fig. 2
    plt.ylabel(r'$r$', fontsize=24)
    plt.xticks(np.linspace(-0.5, 0, 6, endpoint=True), ('\N{MINUS SIGN}0.5', '\N{MINUS SIGN}0.4', '\N{MINUS SIGN}0.3', '\N{MINUS SIGN}0.2', '\N{MINUS SIGN}0.1', '0'), fontsize=20,fontname='Helvetica')
    plt.yticks(np.linspace(-0.5, 0, 6, endpoint=True), ('\N{MINUS SIGN}0.5', '\N{MINUS SIGN}0.4', '\N{MINUS SIGN}0.3', '\N{MINUS SIGN}0.2', '\N{MINUS SIGN}0.1', '0'), fontsize=20,fontname='Helvetica')
    plt.xlim(r_min, r_max-0.01) # this must come after plt.xtics
    plt.ylim(r_min, r_max-0.01)

#    cmap = matplotlib.cm.get_cmap('Pastel1', 6)
#    cmap = matplotlib.cm.get_cmap('Paired_r', 10)
    # choosing a colormap:
    # https://matplotlib.org/stable/users/explain/colors/colormaps.html
    cmap = matplotlib.colors.ListedColormap(['white','plum','tomato','dodgerblue','cyan', 'orange'])

    bounds = np.array([-1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5])
#    norm = matplotlib.colors.BoundaryNorm(boundaries=bounds, ncolors=256)
#    norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
#    vmin, vmax = min(bounds), max(bounds)
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

    plt.imshow(which_is_best, interpolation='nearest', origin='lower', extent=[r_min, r_max-0.01, r_min, r_max-0.01], aspect=1.0, cmap=cmap, norm=norm) # no smoothing. For a technical reason of Python, interpolation='none' produces a smoothed fig when the figure is saved as .eps
#    plt.imshow(which_is_best, interpolation='nearest', origin='lower', extent=[r_min, r_max-0.01, r_min, r_max-0.01], aspect='auto', cmap=cmap, norm=norm, vmin=vmin, vmax=vmax) # no smoothing. For a technical reason of Python, interpolation='none' produces a smoothed fig when the figure is saved as .eps
    # https://hydro.iis.u-tokyo.ac.jp/~akira/page/Python/contents/plot/color/colormap.html

#    plt.colorbar()

    if type==0:
        plt.title('(d)', fontsize=28, x=-0.08, y=1.05)
        plt.text(-0.43, -0.23, r'$\hat{V}_{\{1, 2, 3\}}$', fontsize=24)
    elif type==1:
        plt.title('(e)', fontsize=28, x=-0.08, y=1.05)
        plt.text(-0.43, -0.25, r'$\hat{V}_{\{1, 2, 3\}}$', fontsize=24)
        plt.text(-0.25, -0.15, r'$\hat{V}_{2}$', fontsize=24)
    elif type==2:
        plt.title('(f)', fontsize=28, x=-0.08, y=1.05)
        plt.text(-0.43, -0.23, r'$\hat{V}_{\{1, 3\}}$', fontsize=24)

plt.subplots_adjust(bottom=0.15)
plt.subplots_adjust(left=0.16)

plt.savefig("fig.pdf", bbox_inches='tight')
# https://stackoverflow.com/questions/4042192/reduce-left-and-right-margins-in-matplotlib-plot