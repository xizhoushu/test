# -*- coding: utf-8 -*-
"""
Created on Wed Oct 20 09:47:14 2021

@author: Xi zhou 
@Modified 03/01/2022
"""
import numpy as np
import math
from math import sqrt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.fftpack import fft, fftshift, ifft
from scipy.fftpack import fftfreq
# from utility import *
######################################## Initi#################################
pi = np.pi
g = 9.81

n_freq = 20
n_dir = 11

amplitude = np.zeros((n_freq, 1))
direction = np.zeros((n_dir, 2))  # direction vector for star

###### Wave and ship Info #######
L = 230
model_scale = 31.599

H31 = 10.66 / model_scale  # sea 8
Tp = 0.8*sqrt(L/model_scale) # set as wave length identical to  ship length   0.8*sqrt(lamba)
# Tp = 3
omega_peak = 2*pi/Tp   # Peak frequency omega_0 (rad/s)
angle_mean = 0 # domiant wave direction 
freq_cutoff = 2  # Frequency cutoff factor [1.5, 5], default 3
dir_cutoff = 0  # Direction cutoff [0, 3pi/8], default 0

### Frequence step ###
# Cutoff frequency
omega_max = freq_cutoff * omega_peak
# Frequency step
delta_omega = omega_max / n_freq
# Frequency vector, starting at delta_omega
Omega_vec = np.linspace(1, omega_max, n_freq)
period_wave = 2*pi/Omega_vec
lambda_wave = 1.56*(2*pi/Omega_vec)**2
####  Direction step ####
delta_angle = (pi-2*dir_cutoff)/n_dir
   		
# Start direction
angle_start = angle_mean - pi/2 + delta_angle/2 + dir_cutoff
   	
# Max direction
angle_max = angle_mean + pi/2 - delta_angle/2 - dir_cutoff

angle_range = np.linspace(angle_start,angle_max,n_dir)

delta_theta = pi/n_dir

################################## Function ###################################
# directional spectrum function
# def directionalSpreading(n, theta):
#     k = 1 / np.sqrt(pi) * math.gamma(n / 2 + 1) / math.gamma((n + 1) / 2)
#     D = k * np.cos(theta) ** 2
#     return D


def directionalSpreading(spread,angle,angle_mean,Type):
    # Type 1 from fossen demo
    # Type 2 k*cos^2(angle)
    
    if(Type == 1):
        k = (2**(2*spread-1)*math.factorial(spread)*math.factorial(spread-1))/(pi*math.factorial(2*spread-1))
        D = k * np.cos(angle - angle_mean)**(2*spread)
    if(Type == 2):
        k = 1 / np.sqrt(pi) * math.gamma(spread / 2 + 1) / math.gamma((spread + 1) / 2)
        D = k * np.cos(angle) ** spread
    return D


def spectrum(omiga,H31,Tp,specType):
   if(specType==1):
        A = 0.0081 * g ** 2
        B = 3.11 / (H31 ** 2)
        S = (A / omiga ** 5) * np.exp(-B / omiga ** 4)
   if(specType==2):
       A=487*H31**2/Tp**4
       B=1949/Tp**4;  
       S = (A / omiga ** 5) * np.exp(-B / omiga ** 4)
   if(specType == 3):
    # p1 = Hs  --Significant wave height (Hs = 4 m0^1/2) [m]
    # p2 = T1   --Average wave period (T1 = 2 pi m0/m1)     [s]
       A=173*H31**2/Tp**4
       B=691/Tp**4
       S = A*omiga**(-5) * np.exp(-B / (omiga**4))
     
   return S



def directional_Sprectrum2(S,DirSpre):
    D = np.zeros((len(S), len(DirSpre)))
    for i in range(len(S)):
     for j in range(len(DirSpre)):
         D[i][j] = S[i]* DirSpre[j]
    return D
    

# def directional_Sprectrum(omiga,H31,Tp,specType, theta, n):
#     D = np.zeros((len(omiga), len(theta)))
#     for i in range(len(omiga)):
#         for j in range(len(theta)):
#             D[i][j] = spectrum(omiga[i],H31,Tp,specType) * directionalSpreading(n, theta[j])
#     return D



###############Prepare for plot and generating paramenters of JAVA Macro ##############
# -To compute parameter to be placed in javaMacro
tan_value = np.tan(angle_range)  	
for i in range(n_dir):
    direction[i] = -1, tan_value[i]

# caculate the frequency and directional distribution function
DirSpre = directionalSpreading(spread = 2, angle = angle_range, angle_mean = 0, Type = 2)  # Type 1 Fossem Type 2 k*cos(angle)^2
S = spectrum(Omega_vec,H31,Tp,specType = 1)


# caculate the wave H in [x,y,z]
Position = [0,0,0] # x,y,z
num = 1000
h = np.zeros(num)
time = np.zeros(num)

phase = np.array([np.random.randint(0, 2*pi) for _ in range(n_dir*n_freq) ])
phase = phase.reshape(n_freq,n_dir)
for t in range(num):
    zeta = 0
    time[t] = 0.1*t
    for i in range(n_freq):
        for j in range(n_dir):
            k = Omega_vec[i]**2/g
            
            zeta += np.sqrt(2 * S[i] * DirSpre[j] * delta_theta * delta_omega)*np.sin(Omega_vec[i]*(t-1)/10
                        - k*(np.sin(Position[1])+np.cos(Position[0])) + phase[i][j])
    
    h[t] = zeta
##################### Plot ################

plt.figure(num=1)
plt.subplot(3,1,1)
plt.plot(Omega_vec,S)
plt.subplot(3,1,2)
plt.plot(lambda_wave,S)
plt.subplot(3,1,3)
plt.plot(period_wave,S)
plt.title('frequency spectrum')
plt.figure(num=2)
plt.plot(angle_range,DirSpre)
plt.title('directional function')
### plot the spectrum
X, Y = np.meshgrid(Omega_vec, angle_range)
Z = directional_Sprectrum2(S,DirSpre)



fig = plt.figure(num=3)
ax1 = Axes3D(fig)
# ax1 = plt.axes(projection='3d')
# ax1.plot_surface(omiga.T,angle_range.T,Z)
ax1.plot_surface(X.T, Y.T, Z, rstride=1, cstride=1, cmap="viridis", edgecolor="none")
ax1.set_xlabel("omiga")
ax1.set_ylabel("anlge")
ax1.set_zlabel("S(w,/theta)")
plt.show()



plt.figure(num=4)
plt.subplot(1,2,1)
plt.plot(time,h)
plt.subplot(1,2,2)
plt.boxplot(h)
# plt.subplot(1,3,3)
# plt.plot(np.abs(fft(h)))
plt.show()

### generate java code
# Set the editable index
wave_index = 0
superpositionVofWave_Index = 3
unit_Index = 1
res = ""


#################### Code generation ###################
# with () as f:
with open("log_frequence.txt",'w') as f:

# f = open("log_frequence.txt", "w")
    for j in range(0, n_dir):

        for i in range(0, n_freq):

            S_directioanl = np.sqrt(2 * S[i] * DirSpre[j] * delta_theta * delta_omega)
            str1 = "FirstOrderSuperposingVofWave firstOrderSuperposingVofWave_%d = \n" % (
                wave_index
            )

            str2 = (
                '   superpositionVofWave_%d.getSuperposingVofWaveManager().createSuperposingVofWave(FirstOrderSuperposingVofWave.class, "FirstOrderSuperposing");\n'
                % superpositionVofWave_Index
            )

            str3 = (
                "  firstOrderSuperposingVofWave_%d.getAdvancingDirection().setComponents(%.5f, %.5f, 0.0);\n"
                % ((wave_index), direction[j][0], direction[j][1])
            )

            str4 = "  firstOrderSuperposingVofWave_%d.getAmplitude().setValue(%.5f);\n" % (
                (wave_index),
                S_directioanl,
            )

            str5 = "  firstOrderSuperposingVofWave_%d.getPhase().setUnits(units_%d);" % (
                wave_index,
                unit_Index,
            )

            str6 = "  firstOrderSuperposingVofWave_%d.getPhase().setValue(%d);\n" % (
                (wave_index),
                np.random.randint(0, 360)
            )

            str7 = (
                "  firstOrderSuperposingVofWave_%d.getSpecificationOption().setSelected(VofWaveSpecificationOption.Type.WAVE_PERIOD_SPECIFIED);\n"
                % wave_index
            )

            str8 = (
                "  ((VofWavePeriodSpecification) firstOrderSuperposingVofWave_%d.getVofWaveSpecification()).getWavePeriod().setValue(%.5f);\n"
                % ((wave_index), 2 * pi / Omega_vec[i])
            )

            res = str1 + str2 + str3 + str4 + str5 + str6 + str7 + str8
            # S_directioanl = am
            wave_index += 1
            # print(res, file=f)
            f.write(res)


print(wave_index, res)
