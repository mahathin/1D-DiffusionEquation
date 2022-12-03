# -*- coding: utf-8 -*-
"""
Created on Fri Jan 14 15:20:47 2022

@author: Mahathi
"""

from scipy.special import erfc
import numpy as np
import matplotlib.pyplot as plt


# initials
length = 0.05 #distance upto which we are considering the velocities
t_max = 0.8 #in m
v = 1.6*10**(-4) #kinematic viscosity (m2/s)(this value is for air)
v0 = 10 # initial velocity of particles - diffusion constant in m/s
dx = 0.0006 #our step sinze for x

#Calculated Parameters
dt = (dx**2)/(2*v) #calculating our step size for t - divided by 2 for system to be stable - from explicit method
gamma = (v*dt)/(dx**2) 

#finding numerical solution using fdm

def diffusion_numerical(dt,dx,t_max,length,v,v0):
    x = np.arange(0,length+dx,dx) 
    t = np.arange(0,t_max+dt,dt)
    V = np.zeros([len(t),len(x)]) #velocities
    V[:,0] = v0 #initial velocity
    for i in range(0,len(t)-1): 
        for j in range(1,len(x)-1): 
            V[i+1,j] = V[i,j] + gamma*(V[i,j-1] - 2*V[i,j] + V[i,j+1]) #using fdm
    return x,V

#finding analytical solutions using error function module erfc

def diffusion_analytic(t,length,V0,dx,v): 
    x = np.arange(0,length+dx,dx)
    e1 = length/(2*(t*v)**0.5)
    e = x/(2*(t*v)**0.5)
    s1 = 0
    s2 = 0
    for k in range(0,10000):
        s1 = s1 + erfc(2*k*e1 + e)
        s2 = s2 + erfc(2*(k+1)*e1-e)
    V_analytic = V0*(s1-s2)
    return V_analytic

#calling functions

x,V = diffusion_numerical(dt,dx,t_max,length,v,v0)

#plotting

fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(1,1,1)
ax.set_title('Comparison between numerical and analytic solutions')
ax.set_xlabel('Distance(m)')
ax.set_ylabel('Velocity of particles(m/s)')
    
fig1 = plt.figure(figsize=(10,8))
ax1 = fig1.add_subplot(1,1,1)
ax1.set_title('Error between numerical and analytical solutions')
ax1.set_xlabel('Distance(m)')
ax1.set_ylabel('Velocity difference of particles between numerical and analytical solutions(m/s)')

i = np.arange(0,1.0,0.06) #plotting intervals

for t in i:
    ax.plot(x,V[int(t/dt),:],label='numerical(different colours represent different times)')
    
    if t != 0:
        V_analytic = diffusion_analytic(t,length,v0,dx,v)
        ax.plot(x,V_analytic,'ok',label='analytic', markersize=2)
        
        error = V_analytic - V[int(t/dt),:]
        ax1.plot(x, error)
        
    if t == 0.06:
        ax.legend()

ax.show()
ax1.show()
