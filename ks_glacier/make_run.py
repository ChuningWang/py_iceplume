import numpy as np
import netCDF4 as nc
from scipy.io import loadmat
import subprocess

# ----------------------------------------------------------------------
ctddata = loadmat('KS2014_KS1_TS.mat')
salt = np.nanmean(ctddata['july2829'][0, 0][0], axis=1)
temp = np.nanmean(ctddata['july2829'][0, 0][1], axis=1)
zr = -1*ctddata['july2829'][0, 0][2].squeeze()

salt = salt[::-1][-300:]
temp = temp[::-1][-300:]
zr = zr[::-1][-300:]
zw = np.arange(-301.5, -1, 1)

ntracer = 1
N = len(zr)

dye = np.zeros((ntracer, N))
v = np.zeros(N)
w = np.zeros(N)

np.savetxt('../inputs/iceplume_zw.txt',
        np.array([zw]).T, fmt='%20.10e')
np.savetxt('../inputs/iceplume_zr.txt',
        np.vstack((temp, salt, v, w, dye)).T, fmt='%20.10e')
np.savetxt('../inputs/iceplume_scalar.txt',
        np.array([N, 100, 300, 15, 200, 0, 0, 1, 4, -260, 220]), fmt='%20.10e')

subprocess.call('cd ..; ./build.bash', shell=True)
subprocess.call('cd ..; ./iceplume_test.exe', shell=True)

data_zr = np.loadtxt('../outputs/iceplume_zr.txt', skiprows=1)
data_zw = np.loadtxt('../outputs/iceplume_zw.txt', skiprows=1)
data_dye = np.loadtxt('../outputs/iceplume_dye.txt')
zw = data_zw[:, 1]
f = data_zw[:, 2]
w = data_zw[:, 3]
t = data_zw[:, 4]
s = data_zw[:, 5]
a = data_zw[:, 6]
rho = data_zw[:, 8]
zr = data_zr[:, 1]
dz = data_zr[:, 2]
ent = data_zr[:, 3]
det = data_zr[:, 4]
rhoAm = data_zr[:, 9]
m = data_zr[:, 10]
dye = data_dye[2:]
