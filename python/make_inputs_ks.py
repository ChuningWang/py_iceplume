import numpy as np
from scipy.io import loadmat

ctddata = loadmat('../data/KS2014_KS1_TS.mat')
salt0 = np.nanmean(ctddata['july2829'][0, 0][0], axis=1)
temp0 = np.nanmean(ctddata['july2829'][0, 0][1], axis=1)
zr0 = -1*ctddata['july2829'][0, 0][2].squeeze()

salt0 = salt0[::-1]
temp0 = temp0[::-1]
zr0 = zr0[::-1]
salt0[np.isnan(salt0)] = salt0[498]
temp0[np.isnan(temp0)] = temp0[498]

zw = np.arange(-300, 1, 1)
zr = 0.5*(zw[1:] + zw[:-1])

salt = np.interp(zr, zr0, salt0)
temp = np.interp(zr, zr0, temp0)

# ------------ scalar inputs -------------------------------------------
N = len(zr)
dx = 100.
dy = 300.
dt = 20.
sg_runoff = 230.
sg_temp = 0.
sg_salt = 0.
sg_dye01 = 1.
sg_typ = 4
sg_dep = -260.
sg_len = 220.

# ------------ load profiles -------------------------------------------
dye01 = np.zeros(N)
v = np.zeros(N)
w = np.zeros(N)

np.savetxt('../inputs/iceplume_zw.txt', np.array([zw]).T, fmt='%20.10e')
np.savetxt('../inputs/iceplume_zr.txt',
        np.array([temp, salt, v, w, dye01]).T, fmt='%20.10e')
np.savetxt('../inputs/iceplume_scalar.txt',
        np.array([[N, dx, dy, dt,
                   sg_runoff, sg_temp, sg_salt, sg_dye01,
                   sg_typ, sg_dep, sg_len]]),
        fmt='%03d ' + 9*'%20.10e ' + '%20.10e')
