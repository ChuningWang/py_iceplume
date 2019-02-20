import numpy as np
import matplotlib.pyplot as plt

data_zr = np.loadtxt('../outputs/iceplume_zr.txt', skiprows=1)
data_zw = np.loadtxt('../outputs/iceplume_zw.txt', skiprows=1)

zr = data_zr[:, 1]
zw = data_zw[:, 1]

rhoAm = data_zr[:, 8]
rho   = data_zw[:, 8]

plt.plot(rhoAm, zr, 'r')
plt.plot(rho,   zw, 'k')
