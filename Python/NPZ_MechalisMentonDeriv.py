# -*- coding: utf-8 -*-
"""
Created on Sat May 23 14:41:50 2020

@author: Jasen
"""
# One dimensional NPZ model computed with Euler Method
import numpy as np

def MM_Deriv(Nute, Phyt, Zoo, depth, Vm, Ks, Kext, Rm, Lam, gam, m, g):
   """Use Euler Method to compute 1D NPZ with Mechalis-Menton model """
  
   # Pre-allocate
   dP_dt = np.empty((depth.size,1))
   dZ_dt = dP_dt
   dN_dt = dP_dt
   
   for iz, z in enumerate(depth) :
      
       #compute terms
       NuteUp = (Vm*Nute[iz]*Phyt[iz]/(Ks+Nute[iz]))*np.exp(Kext*z)
       Graze = -Rm*Zoo[iz]*(Phyt[iz]**2)/((Lam**2)+Phyt[iz]**2)
       ZooAssim = (1 - gam)*Graze
       SlopFeed = gam*Graze
       PhytMort = m*Phyt[iz]
       ZooMort = g*Zoo[iz]
       
       #compute derivative
       dP_dt[iz] = NuteUp - Graze - PhytMort
       dZ_dt[iz] = ZooAssim - ZooMort
       dN_dt[iz] = -NuteUp + SlopFeed + PhytMort + ZooMort

   #store output 
   Derivs = (dN_dt, dP_dt, dZ_dt)
   
   return Derivs
