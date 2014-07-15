import matplotlib
matplotlib.use('Agg') ## TO WRITE PDF'S AND PNGS
from numpy import fft
import numpy as np
import matplotlib.pyplot as plt
fig = matplotlib.pyplot.gcf()
fig.set_size_inches(12.5,10.5)
##### LOAD THE ELECTRICAl SIGNAL YOU WANT ####
########  TO FOURIER TRANSFORM ###############
fx=np.loadtxt('dat1.dat')
n = 2048   		# Number of data points
dx = 89e-6 		# Sampling period (in meters)
x = dx*np.arange(0,n)   # x coordinates
w1 = 100.0 		# wavelength (meters)
w2 = 20.0 		# wavelength (meters)
Fk = fft.fft(fx)        #AMPLITUDE IN FREQUENCY DOMAIN 
nu = fft.fftfreq(n,dx)  # Natural frequencies
Fk = fft.fftshift(Fk)   # Shift zero freq to center
Fkabs=np.absolute(Fk)**2# Power spectrum
nu = fft.fftshift(nu)   # Shift zero freq to center
combined = np.vstack((nu,Fkabs)).T
np.savetxt('out1.dat',combined,delimiter=" ")
##################################################
plt.subplot(311)
plt.plot(x,fx,'r-')
plt.ylabel('Voltage V(t)')
plt.xlabel('Time (sec)')
plt.grid(True)

plt.subplot(312)
plt.plot(nu,Fk,'r-')
plt.xlim(-1000,1000)
plt.xlabel('Frequency (rad/sec)')
plt.ylabel('Voltage V(w)')
plt.grid(True)

plt.subplot(313)
plt.plot(nu,Fkabs,'r-')
plt.xlim(-1000,1000)
plt.xlabel('Frequency (rad/sec)')
plt.ylabel('Power Spectrum |V(w)|^2')
plt.grid(True)
plt.savefig('ffts.png',figsize=(8, 9))
#plt.show()
