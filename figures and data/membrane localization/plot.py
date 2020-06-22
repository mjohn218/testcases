# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 15:47:21 2020

@author: phenning
"""
import numpy as np
import matplotlib.pyplot as plt

from collections import OrderedDict

linestyles = OrderedDict(
    [('solid',               (0, ())),
     ('loosely dotted',      (0, (1, 10))),
     ('dotted',              (0, (1, 5))),
     ('densely dotted',      (0, (1, 1))),

     ('loosely dashed',      (0, (5, 10))),
     ('dashed',              (0, (5, 5))),
     ('densely dashed',      (0, (5, 1))),

     ('loosely dashdotted',  (0, (3, 10, 1, 10))),
     ('dashdotted',          (0, (3, 5, 1, 5))),
     ('densely dashdotted',  (0, (3, 1, 1, 1))),

     ('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),
     ('dashdotdotted',         (0, (3, 5, 1, 5, 1, 5))),
     ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1)))])

import matplotlib.ticker as mticker
f = mticker.ScalarFormatter(useMathText=True)
f.set_powerlimits((-3,3))
"${}$".format(f.format_data(0.0001))

default_x = 5*1.4
default_y = 3*1.4

plt.gcf().set_size_inches(default_x, default_y)

#%%
## load data ##
ode_mem = np.loadtxt('new_run/ODE.txt', skiprows=1)
ode_mem[:,1:9] *= 0.2209

pde_mem = np.loadtxt('new_run/PDE.txt', skiprows=1)

ssa_mem = np.loadtxt('new_run/SSA_MA_MB.txt', delimiter=',')

smoldyn_mem = np.loadtxt('Smoldyn_MA_MB_large.txt', delimiter=',')

fpr_mem = np.loadtxt('Re__Data_etc_/Avg10trajAMandBM_FPR.txt')

## plot figure ##
fig, ax = plt.subplots()
ax.plot(ode_mem[:,0], ode_mem[:,4], '--', label='ODE', color = 'tab:blue')
ax.plot(pde_mem[:,0], pde_mem[:,4], '-.', label = 'PDE', color = 'tab:green')
ax.plot(fpr_mem[:,0], fpr_mem[:,1], linestyle=linestyles['densely dashdotdotted'], label='NERDSS', color = 'tab:red')
ax.plot(ssa_mem[0,:], ssa_mem[1,:], ':', label='Gillespie', color = 'tab:purple')
ax.plot(smoldyn_mem[:,0], smoldyn_mem[:,1], linestyle=linestyles['densely dashed'], label = 'Smoldyn', color = 'tab:orange')
plt.xlim(0, 2)
plt.legend(loc=1, labelspacing = 0.3)
plt.xlabel('time [s]');
plt.ylabel('MA(t)');

axins = ax.inset_axes([0.25, 0.6, 0.37, 0.37])
axins.set_xlim(1.2, 1.75)
axins.set_ylim(50, 65)
ax.indicate_inset_zoom(axins)
axins.plot(ode_mem[:,0], ode_mem[:,4], '--', label='ODE', color = 'tab:blue')
axins.plot(pde_mem[:,0], pde_mem[:,4], '-.', label = 'PDE', color = 'tab:green')
axins.plot(fpr_mem[:,0], fpr_mem[:,1], linestyle=linestyles['densely dashdotdotted'], label='NERDSS', color = 'tab:red')
axins.plot(ssa_mem[0,:], ssa_mem[1,:], ':', label='Gillespie', color = 'tab:purple')
axins.plot(smoldyn_mem[:,0], smoldyn_mem[:,1], linestyle=linestyles['densely dashed'], label = 'Smoldyn', color = 'tab:orange')

#plt.savefig("plot_membrane_location.pdf",bbox_inches='tight', dpi = 400)
#plt.savefig("plot_membrane_location.svg",bbox_inches='tight', dpi = 400)




