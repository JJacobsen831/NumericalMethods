# -*- coding: utf-8 -*-
"""
Created on Sat May 23 14:41:50 2020

@author: Jasen
"""
import os
os.chdir("/Users/Jasen/Documents/GitHub/NumericalMethods/Python")

import matplotlib.pyplot as plt
import numpy as np

from NPZ_MechalisMentonDeriv import MM_Deriv
from NumericalMethodFuns import EulerMethod as EM

#Run Mechalis-Menton NPZ model with Euler Method

#Set initial condition values
N_tot = 10.
Nute0 = N_tot/2
Phyt0 = N_tot/5
Zoo0 = N_tot - Nute0 - Phyt0

#set depth
nz = 10
z_bot = -65.
depth = np.linspace(z_bot, 0., nz)

#time step paramters
sec_per_day = 86400.
maxdays = 200.
tau = sec_per_day/500.   # <- Time Step
maxtime = maxdays*sec_per_day
n_step = int(maxtime/tau)


#save interval
dtsave_days = 0.1
dtsave_sec = dtsave_days*sec_per_day
dtsave_step = np.round(dtsave_sec/tau, decimals = 0)
nsaves = int(np.round(n_step/dtsave_step, decimals = 0))

#save variables
t_save = np.ones((1,nsaves), dtype = float)
N_save = np.ones((nz, nsaves), dtype = float)
P_save = N_save
Z_save = N_save

#model parameters
Vm = 2.0/sec_per_day     # max nutrient uptake rate
Ks = 0.1                 # Michaelis-menton half saturation value
Kext = 0.06              # e-folding scale of PAR with depth
Rm = 0.5/sec_per_day     # max grazing rate
Lam = 0.2                # level of saturated grazing
gam = 0.3                # percent of sloppy grazing (1-Gamma) is assimialtion eff.
m = 0.1/sec_per_day      # phytoplankton mortality
g = 0.2/sec_per_day      # zooplankton mortality

#define initial conditions
Nute = np.ones((nz, 1), dtype = float)*Nute0
Phyt = np.ones((nz, 1), dtype = float)*Phyt0
Zoo = np.ones((nz, 1), dtype = float)*Zoo0

#save initial condition to output array
N_save[:, 0] = Nute[:,0]
P_save[:, 0] = Phyt[:,0]
Z_save[:, 0] = Zoo[:,0]

#time stepping
isave = 0
tstep = np.linspace(0, int(maxtime), n_step)
for t_idx, istep in enumerate(tstep):
    time = istep*tau
    
    #compute derivatives
    (dNdt, dPdt, dZdt) = MM_Deriv(Nute, Phyt, Zoo, depth, Vm, Ks, Kext, Rm, Lam, gam, m, g)
    
    #Euler Method
    Nute = EM(Nute, dNdt, tau)
    Phyt = EM(Phyt, dPdt, tau)
    Zoo = EM(Zoo, dZdt, tau)
    
    if np.mod(istep ,dtsave_step) == 0:
        t_save[:,isave] = time
        N_save[:,isave]= Nute.T
        P_save[:,isave] = Phyt.T
        Z_save[:,isave] = Zoo.T
        isave = isave+1
        print((istep, isave, dNdt - dPdt))

plt.figure()
plt.plot(t_save.T, N_save[5,:])
plt.plot(t_save.T, P_save[5,:])
plt.plot(t_save.T, Z_save[5,:])
plt.show

