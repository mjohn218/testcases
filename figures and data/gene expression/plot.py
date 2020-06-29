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
smoldyn2 = np.loadtxt('smoldyn_200s_uniform.dat')

fig, ax = plt.subplots(2)
ax[0].plot(time[0:390], ode_r[0:390], label='R', color='tab:orange')
ax[0].plot(time, ode_a, label='A',color='tab:blue')
ax[0].legend(loc=9, labelspacing = 0.3, ncol =2)
ax[0].set_xlabel('time [s]');
ax[0].set_ylabel('N(t)');
ax[0].set_xlim(0, 200)
ax[0].set_ylim(0, 1900)

#ax[1].plot(time, smoldyn_r, label='R', color='tab:orange')
#ax[1].plot(time, smoldyn_a,color='tab:blue', label='A')
ax[1].plot(smoldyn2[:,0], smoldyn2[:,2], label='R', color='tab:orange')
ax[1].plot(smoldyn2[:,0], smoldyn2[:,1],color='tab:blue', label='A')
plt.legend(loc=9, labelspacing = 0.3, ncol =2)
ax[1].legend(loc=9, labelspacing = 0.3, ncol =2)
ax[1].set_xlabel('time [s]');
ax[1].set_ylabel('N(t)');
ax[1].set_xlim(0, 200)
ax[1].set_ylim(0, 1900)

#plt.savefig("gene_expression.pdf",bbox_inches='tight', dpi = 400)
#plt.savefig("gene_expression.svg",bbox_inches='tight', dpi = 400)


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
plt.xlim(0, 200)
plt.ylim(0, 1900)

#%%

pde_fix = np.genfromtxt('FixedPromoter_DA2_DR50_ARC.csv', delimiter  = ',',skip_header=9)

dist_fix = pde_fix[0,3:27]
dist_fix[0] = 0
time_fix = pde_fix[1:202,1]
Aout = pde_fix[1:202,4]
Aout = Aout*602.2*4
Acenter = pde_fix[1:202,14]
Acenter = Acenter*602.2*4

fig, ax = plt.subplots()
ax.plot(time, ode_a, label='ODE', color='tab:blue')
ax.plot(time_fix, Aout, label='PDE, periphery', color='tab:green')
ax.plot(time_fix, Acenter, label='PDE, center', color='tab:red')
plt.xlabel('time [s]');
plt.ylabel('A(t)');
plt.legend(loc=9, labelspacing = 0.3, ncol =3)
plt.xlim(0, 200)
plt.ylim(0, 1700)

default_x = 5*1.4
default_y = 1.5*1.4
plt.gcf().set_size_inches(default_x, default_y)

plt.savefig("gene_expression_space.pdf",bbox_inches='tight', dpi = 400)
plt.savefig("gene_expression_space.svg",bbox_inches='tight', dpi = 400)--

