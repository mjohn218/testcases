# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 17:01:00 2020

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
plt.rcParams.update({'font.size': 20})


#%%
########## Plot MinCDE PDE ##########
default_x = 5*1.4*2
default_y = 3*1.4

#minDt_b_pde = np.genfromtxt('minDt_WallOnly__Curves_vars(1)_times(301).csv', delimiter  = ',',skip_header=9)
#time_pde = minDt_b_pde[1:302,1]
#space_pde = minDt_b_pde[0,2:113]
#minDt_b_pde = minDt_b_pde[1:,2:]


minDt_b_pde1 = np.loadtxt('MinCDE_PDE\pde_minDt_1.txt')
time_pde = np.loadtxt('MinCDE_PDE\minD_time.txt')

minDt_b_pde2 = np.loadtxt('MinCDE_PDE\pde_minDt_2.txt')

minDt_b_pde1_fine = np.loadtxt('MinCDE_PDE\minDt_1_fine.txt', skiprows=1)
time_pde_fine = minDt_b_pde1_fine[:,0]
minDt_b_pde1_fine = minDt_b_pde1_fine[:,2]

minDt_b_pde2_fine = np.loadtxt('MinCDE_PDE\minDt_2_fine.txt', skiprows=1)
minDt_b_pde2_fine = minDt_b_pde2_fine[:,2]

fig, ax = plt.subplots(2)
ax[0].plot(time_pde, minDt_b_pde1)
ax[0].plot(time_pde, minDt_b_pde2)
ax[0].set_xlim(0, 600)
ax[0].set_ylim(0, 400)
ax[0].get_xaxis().set_visible(False)
#ax[0].set_xlabel('time [s]');
ax[0].set_ylabel(r'$\rho_{MinD-ATP} [1/\mu m^2]$');

ax[1].plot(time_pde_fine, minDt_b_pde1_fine)
ax[1].plot(time_pde_fine, minDt_b_pde2_fine)
ax[1].set_xlim(0, 600)
ax[1].set_ylim(0, 400)
ax[1].set_xlabel('time [s]');
#ax[1].set_ylabel(r'$\rho_{MinE-MinD-ATP} [1/\mu m^2]$');
plt.gcf().set_size_inches(default_x, default_y)

plt.savefig("Min_pde.svg",bbox_inches='tight', dpi = 400)
plt.savefig("Min_pde.pdf",bbox_inches='tight', dpi = 400)

########## Plot MinCDE Smoldyn ##########
default_x = 5*1.4*4/9
default_y = 3*1.4



minDt_b_smoldyn = np.genfromtxt('run5/NDtb.txt', delimiter = ',')
EminDt_smoldyn = np.genfromtxt('run5/NEDt.txt', delimiter = ',')
time_smoldyn = np.genfromtxt('run5/time.txt', delimiter = ',')
space_smodyn = np.genfromtxt('run5/space.txt', delimiter = ',')


fig, ax = plt.subplots()
ax.plot(time_smoldyn, minDt_b_smoldyn[0, :])
ax.plot(time_smoldyn, minDt_b_smoldyn[-1, :])
ax.set_xlim(0, 200)
ax.set_ylim(0, 800)
ax.set_xticklabels([])

plt.ylabel(r'$\rho_{MinD-ATP} [1/\mu m^2]$');
plt.gcf().set_size_inches(default_x, default_y)

#plt.savefig("Min_smoldyn.svg",bbox_inches='tight', dpi = 400)
#plt.savefig("Min_smoldyn.pdf",bbox_inches='tight', dpi = 400)

########## Plot MinCDE Snapshots ##########
EminDt_b_pde = np.genfromtxt('Emin_WallOnly__Curves_vars(1)_times(301).csv', delimiter  = ',',skip_header=9)
EminDt_b_pde = EminDt_b_pde[1:,2:]

plt.rcParams.update({'font.size': 16})
fig, ax = plt.subplots(1,2)
ax[0].plot(space_pde, minDt_b_pde[71,:], label = 'MinD-ATP')
ax[0].plot(space_pde, EminDt_b_pde[71,:], label = 'E-MinD-ATP')
ax[0].set_xlim(0, 4)
ax[0].set_ylim(0, 550)
ax[0].set_xlabel(r'$z \;[\mu m]$');
ax[0].set_ylabel(r'$\rho \;[1/\mu m^2]$');
ax[0].legend(loc=0, labelspacing = 0.3)
ax[0].text(0.2, 500, 't < 250 s');

ax[1].plot([],[], ' ', label = 'PDE')
ax[1].plot(space_pde, minDt_b_pde[155,:], label = 'MinD-ATP', color = 'tab:blue')
ax[1].plot(space_pde, EminDt_b_pde[155,:], label = 'E-MinD-ATP', color = 'tab:orange')
ax[1].plot([],[], ' ', label = 'Smoldyn')
ax[1].plot(space_smodyn[0:-1], minDt_b_smoldyn[:,39], label = 'MinD-ATP', color = 'blue')
ax[1].plot(space_smodyn[0:-1], EminDt_smoldyn[:,39], label = 'E-MinD-ATP', color = 'red') 
ax[1].set_xlim(0, 4)
ax[1].set_ylim(0, 550)
ax[1].set_xlabel(r'$z \;[\mu m]$');
ax[1].legend(loc=0, labelspacing = 0.3)
ax[1].text(0.2, 500, 't > 250 s')
ax[1].set_yticklabels([])

default_x = 10*1.4
default_y = 3*1.4
plt.gcf().set_size_inches(default_x, default_y)
plt.subplots_adjust(wspace=0.05)

#plt.savefig("PAPER_Min_vs_z_1.svg",bbox_inches='tight', dpi = 400)
#plt.savefig("PAPER_Min_vs_z_1.pdf",bbox_inches='tight', dpi = 400)





