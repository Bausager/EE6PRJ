# -*- coding: utf-8 -*-
"""
Created on Sun Oct  3 12:23:24 2021

@author: Bausa
"""

import numpy as np
import matplotlib.pyplot as plt

error_old = 0.0
i_part_error_old = 0

Ts = 0.01
kp = 42
Ti = 10
Td = 0.001

def p_part(error):
    error_old = error * kp
    return (error_old)


def i_part_martin(error, error_old):
    error_old = error_old + (Ts/2.0) * (error - error_old)
    return (error_old/Ti)

def i_part(error, error_old):
    error_old = error_old + (Ts/(2.0*Ti)) * (error - error_old)
    return (error_old)



def d_part_martin(error, error_old):    
    error_old = -1*error_old + (2.0/Ts) * (error - error_old)
    return (error_old*Td)


def d_part(error, error_old):
    error_old = -1*error_old + ((2.0*Td)/Ts) * (error - error_old)
    return (error_old)


N = 100

n = int(N)

x =  np.sin(np.linspace(0, 2*np.pi, num=int(N)))
#x = np.append(np.append(np.zeros(N), np.ones(N)), np.zeros(N))
y = np.zeros(n)
y1 = np.zeros(n)


#plt.plot(x)
for i in range(0, n-1):
    error = x[i] - error_old
    error_old = d_part_martin(error, error_old)
    y[i] = error_old
    
error_old = 0

for i in range(0, n-1):
    error = x[i] - error_old
    #print(error)
    error_old = d_part(error, error_old)
    y1[i] = error_old
    
#%%
fig, ax = plt.subplots(1, 3)

ax[0].plot(y1-y, label='y1-y')
ax[0].grid()
ax[0].legend()

ax[1].plot(y, label='y')
ax[1].grid()
ax[1].legend()

ax[2].plot(y1, label='y1')
ax[2].grid()
ax[2].legend()
#%%
import numpy as np
import matplotlib.pyplot as plt


def G(frq, Tw):
    num = np.sin((frq*Tw)/2.0)
    den = (frq*Tw)/2.0
    gain = np.abs(num/den)
    
    ang = -((frq*Tw)/2.0)
    
    return gain, ang


def G1(frq, Tw):
    num = np.sin((frq*2*np.pi*Tw)/2.0)
    den = (frq*2*np.pi*Tw)/2.0
    gain = np.abs(num/den)
    
    #ang = -((frq*2*np.pi*Tw)/2.0)
    ang = -180.0 * f * Tw
    
    return gain, ang
    
    
Tw = 0.02


N = 1000
n = 100000
f = (np.linspace(start=0, stop=N, num=n))

Gain = np.zeros(N)
Angle = np.zeros(N)


Gain, Angle = G1(f, Tw)


Gain = 20*np.log10(Gain)
#Angle = Angle * 180/np.pi

f[Gain < -3]


print(-180.0 * 60 * 0.02)



#%%
#print(f[Gain < -3][0])
#print(Gain[f > 30][0])

#print(f'Phase Shift at: {f[Angle < -50][0]}')
#print(Angle[f > f[Angle > -50][-1]][0])

x = 0

for i in range(0, len(Angle)):
    Angle[i] = Angle[i] - x
    if Angle[i] <= -180:
        x = x - 180
        
        #print(f'Ang: {Angle[i]}, x {x}')



#plt.semilogx(f, Gain, label="Gain")
plt.semilogx(f, Angle, label="Angle")
plt.grid()
plt.legend()









#%%

x = (1/(2.363*10**-4) * 1/np.sqrt(2) * 2) / (2*np.pi)

print(x)




















