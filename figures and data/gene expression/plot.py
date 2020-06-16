# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 15:43:39 2020

@author: phenning
"""

import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict
import scipy.io as sio

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

########## Plot gene expression ##########

gene_data = sio.loadmat('Original/MCellData/matlab.mat')
time = gene_data["excelt_osci"]
ode_a = gene_data["ode_osci_a"]
ode_r = gene_data["ode_osci_r"]
smoldyn_a = gene_data["smol_oci_a"]
smoldyn_r = gene_data["smol_osci_r"]
mcell_r = gene_data["mcelldat_r"]

fig, ax = plt.subplots()
ax.plot(time, ode_r, label='R, ODE', color='tab:orange')
ax.plot(time, ode_a, label='A, ODE',color='tab:blue',)
ax.plot(time, smoldyn_r, label='R, Smoldyn', color='tab:orange', linestyle=linestyles['densely dotted'],)
ax.plot(time, smoldyn_a,color='tab:blue', label='A, Smoldyn', linestyle=linestyles['densely dotted'],)
plt.legend(loc=9, labelspacing = 0.3, ncol =2)
plt.xlabel('time [s]');
plt.ylabel('N(t)');
plt.xlim(0, 100)
plt.ylim(0, 1900)

plt.savefig("gene_expression.pdf",bbox_inches='tight', dpi = 400)
plt.savefig("gene_expression.svg",bbox_inches='tight', dpi = 400)


#%%

########## Plot gene expression ##########

gene_data2 = sio.loadmat('GRAPH_AND_DATA/matlab.mat')
ode_AR = gene_data2["ARode"]
time_ode = ode_AR[:,0]
ode_a = ode_AR[:,1]
ode_r = ode_AR[:,2]

smoldyn_AR = gene_data2["ARstoch"]
time_smoldyn = smoldyn_AR[:,0]
smoldyn_a = np.mean(smoldyn_AR[:,1:21:2],1)
smoldyn_r = np.mean(smoldyn_AR[:,2:21:2],1)

fig, ax = plt.subplots()
ax.plot(time_ode, ode_r, label='R, ODE', color='tab:blue')
ax.plot(time_ode, ode_a, label='A, ODE',color='tab:blue', linestyle=linestyles['densely dotted'],)
ax.plot(time_smoldyn, smoldyn_r, label='R, Smoldyn', color='tab:orange')
ax.plot(time_smoldyn, smoldyn_a,color='tab:orange', label='A, Smoldyn', linestyle=linestyles['densely dotted'],)
plt.legend(loc=9, labelspacing = 0.3, ncol =2)
plt.xlabel('time [s]');
plt.ylabel('N(t)');
plt.xlim(0, 100)
plt.ylim(0, 1900)