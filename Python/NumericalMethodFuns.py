# -*- coding: utf-8 -*-
"""
Created on Sat May 23 17:41:58 2020

@author: Jasen
"""
#compilation of numerical methods
def EulerMethod(Val, d_dt, tau):
    """Euler Numerical Method"""
    new_val = Val+d_dt*tau
    
    return new_val



def RK_2ndOrder_Heun(Val, d_dt, tau):
    """2nd order Runge-Kutta using Heun's method"""
    #set Runge-Kutta parameters
    a2 = 0.5
    a1 = 1 - a2
    p1 = 1/(2*a2)
    q11= 1/(2*a2)
    
    
