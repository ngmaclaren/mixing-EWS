# Generate figure S29 in the manuscript.
# The output figure is saved in fig.pdf

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from scipy import linspace
import math # sqrt

Nsamples = 200

type=2

w = 0.1 # coupling strength

if type==0:
    r = -0.4
elif type==1:
    r = 2*w*(-w-1-np.sqrt(w*(w+1)))
    print('r = r_c^prime =', r) # -0.286332495807108
elif type==2:
    r = 2*w*(-w-1+np.sqrt(w*(w+1))) # -0.15366750419289202
    print('r_p =', r)
elif type==3:
    r = -0.1535

aspect_ratio = 3.5/2.2
if type==0 or type==1:
    xmax = 2.2
    xmin = -xmax
    ymax = xmax * aspect_ratio
    ymin = -ymax
elif type==2:
    xmax = 0.7
    xmin = -xmax
    ymax = xmax * aspect_ratio
    ymin = -ymax
elif type==3:
    dx = 0.024
    xmax = -1-r/(2*w) + dx
    xmin = -1-r/(2*w) - dx 
    ymax = dx * aspect_ratio
    ymin = -ymax

plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)

x1_a = np.linspace(xmin, xmax, Nsamples)
x2_a = - x1_a**2 / w  - 1 - r / w
x2_b = np.linspace(ymin, ymax, Nsamples)
x1_b = - x2_b**2 / (2*w) - 1 - r / (2*w)

# necessary to get the fonts right
#plt.rcParams.update({
#     "text.usetex": True,
#"font.family": "sans-serif",
#"font.sans-serif": ["Times"]})

plt.rcParams["text.usetex"] = True # then the font for text labels also change (not x/ylabel or x/yticks)
plt.rcParams["font.family"] = 'serif'
plt.rcParams["font.serif"] = ['Times New Roman']

# plt.rcParams["font.family"] = "sans-serif"
# plt.rcParams["font.sans-serif"] = "DejaVu Sans"

# rc('font',**{'family':'serif','serif':['Times']})
# rc('text', usetex=True)

plt.plot(x1_a, x2_a, 'b-', linewidth=1.5, markerfacecolor='none')
plt.plot(x1_b, x2_b, 'm-', linewidth=1.5, markerfacecolor='none')

plt.plot(0, -1-r/w, marker="o", markersize=7, markerfacecolor="black", markeredgewidth=0) # , markeredgecolor="black", markerfacecolor="black")
plt.plot(np.sqrt(-w-r), 0, marker="o", markersize=7, markerfacecolor="black", markeredgewidth=0) # , markeredgecolor="black", markerfacecolor="black")
plt.plot(-np.sqrt(-w-r), 0, marker="o", markersize=7, markerfacecolor="black", markeredgewidth=0) # , markeredgecolor="black", markerfacecolor="black")
plt.plot(-1-r/(2*w), 0, marker="o", markersize=7, markerfacecolor="black", markeredgewidth=0) # , markeredgecolor="black", markerfacecolor="black")
if -2*w-r >= 0:
    plt.plot(0, np.sqrt(-2*w-r), marker="o", markersize=7, markerfacecolor="black", markeredgewidth=0) # , markeredgecolor="black", markerfacecolor="black")
    plt.plot(0, -np.sqrt(-2*w-r), marker="o", markersize=7, markerfacecolor="black", markeredgewidth=0) # , markeredgecolor="black", markerfacecolor="black")

if type==0:
    plt.text(0.05, -1-r/w, r'${\rm V}_1$', fontsize=14)
    plt.text(np.sqrt(-w-r)-0.23, 0.16, r'${\rm P}_1^+$', fontsize=14)
    plt.text(-1-r/(2*w)-0.05, 0.12, r'${\rm V}_2$', fontsize=14)
    plt.text(0.03, -np.sqrt(-2*w-r) - 0.29, r'${\rm P}_2^-$', fontsize=14) # -2*w-r >=0 is satisfied with this r value
    plt.text(-np.sqrt(-w-r)-0.21, 0.17, r'${\rm P}_1^-$', fontsize=14)
    plt.text(2/2.2*xmin, 2.95/3.5*ymax, "(a)", fontsize=28)
    plt.xticks(np.arange(-2, 2.01, 1.0), fontsize=18, fontname='Times New Roman')
    plt.yticks(np.arange(-3, 3.01, 1.0), fontsize=18, fontname='Times New Roman')
elif type==1:
    plt.text(0.05, -1-r/w, r'${\rm V}_1$', fontsize=14)
    plt.text(-1-r/(2*w)-0.03, 0.14, r'${\rm V}_2 (= {\rm P}_1^+)$', fontsize=14)
    plt.text(0.03, -np.sqrt(-2*w-r) - 0.32, r'${\rm P}_2^-$', fontsize=14) # -2*w-r >=0 is satisfied with this r value
    plt.text(-np.sqrt(-w-r)+0.05, 0.05, r'${\rm P}_1^-$', fontsize=14)
    plt.text(2/2.2*xmin, 2.95/3.5*ymax, "(b)", fontsize=28)
    plt.xticks(np.arange(-2, 2.01, 1.0), fontsize=18, fontname='Times New Roman')
    plt.yticks(np.arange(-3, 3.01, 1.0), fontsize=18, fontname='Times New Roman')
elif type==2:
    plt.text(0.01, -1-r/w + 0.03, r'${\rm V}_1$', fontsize=14)
    plt.text(np.sqrt(-w-r), 0.07, r'${\rm P}_1^+$', fontsize=14)
    plt.text(-1-r/(2*w)+0.01, -0.13, r'${\rm V}_2$', fontsize=14)
    plt.text(-1-r/(2*w)+0.01, -0.26, r'$(= {\rm P}_1^-)$', fontsize=14)
    plt.text(2/2.2*xmin, 2.95/3.5*ymax, "(c)", fontsize=28)
    plt.xticks(np.arange(-0.5, 0.501, 0.5), ('\N{MINUS SIGN}0.5', '0', '0.5'), fontsize=18, fontname='Times New Roman')
    plt.yticks(np.linspace(-1, 1, 5, endpoint=True), ('\N{MINUS SIGN}1', '\N{MINUS SIGN}0.5', '0', '0.5', '1'), fontsize=18, fontname='Times New Roman')
elif type==3:
    plt.text(-1-r/(2*w)-0.0025, 0.002, r'${\rm V}_2$', fontsize=14)
    plt.text(-1-r/(2*w)+0.0015, -0.005, r'${\rm P}_1^-$', fontsize=14)
    plt.text(xmin + 0.2/4.4 * (xmax-xmin), 2.95/3.5*ymax, "(d)", fontsize=28)
    plt.yticks(np.arange(-0.03, 0.0301, 0.01), ('\N{MINUS SIGN}0.03', '\N{MINUS SIGN}0.02', '\N{MINUS SIGN}0.01', '0', '0.01', '0.02', '0.03'), fontsize=18, fontname='Times New Roman')
    plt.xticks(np.arange(-0.25, -0.21, 0.01), ('\N{MINUS SIGN}0.25', '\N{MINUS SIGN}0.24', '\N{MINUS SIGN}0.23', '\N{MINUS SIGN}0.22', '\N{MINUS SIGN}0.21'), fontsize=18, fontname='Times New Roman')


if -2*w-r >= 0:
    plt.text(0.03, np.sqrt(-2*w-r) + 0.1, r'${\rm P}_2^+$', fontsize=14)

if type==0 or type==1 or type==2:
    plt.arrow(0, ymin, 0, 1.95*ymax, head_width=xmax*0.02, head_length=ymax*0.02, linewidth=1, color='black')
    plt.arrow(xmin, 0, 1.95*xmax, 0, head_width=ymax*0.02, head_length=xmax*0.02, linewidth=1, color='black')
elif type==3:
    plt.arrow(0, ymin, 0, 1.95/2.0 * (ymax-ymin), head_width=dx*0.02, head_length=ymax*0.02, linewidth=1, color='black')
    plt.arrow(xmin, 0, 1.95/2.0 *(xmax -xmin), 0, head_width=ymax*0.02, head_length=dx*0.02, linewidth=1, color='black')


plt.yticks(fontsize=18)
plt.xlabel(r'$x_1$', fontsize=24, fontname='Times New Roman')

if type==2:
    plt.ylabel(r'$x_2$', fontsize=24, labelpad=-8) # default value of labelpad = 4
if type==3:
    plt.ylabel(r'$x_2$', fontsize=24, labelpad=-18, fontname='Times New Roman') # default value of labelpad = 4
else:
    plt.ylabel(r'$x_2$', fontsize=24)

plt.subplots_adjust(bottom=0.145)
plt.subplots_adjust(left=0.14)

# remove plot box borders
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)

plt.savefig("fig.pdf")