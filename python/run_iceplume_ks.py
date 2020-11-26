import subprocess
import numpy as np
import netCDF4 as nc
from scipy.io import loadmat
import pyroms
import gsw

ctd_file  = '../data/KS2014_KS1_TS.mat'
hist_file = '/Users/cw686/roms_outputs/fjord_ks/outputs/fjord_avg.nc'
grid_file = '/Users/cw686/roms_archive/fjord_ks/fjord_grid.nc'

grd = pyroms.grid.get_ROMS_grid('grd_temp', hist_file=hist_file,
                                grid_file=grid_file)
zr  = grd.vgrid.z_r[:][:, 0, 0]
zw  = grd.vgrid.z_w[:][:, 0, 0]

ctddata = loadmat('../data/KS2014_KS1_TS.mat')
salt0   = np.nanmean(ctddata['july2829'][0, 0][0], axis=1)
temp0   = np.nanmean(ctddata['july2829'][0, 0][1], axis=1)
zr0     = -1*ctddata['july2829'][0, 0][2].squeeze()

salt0 = salt0[::-1]
temp0 = temp0[::-1]
zr0   = zr0[::-1]
salt0[np.isnan(salt0)] = salt0[498]
temp0[np.isnan(temp0)] = temp0[498]

# zw = np.arange(-300, 1, 1)
# zr = 0.5*(zw[1:] + zw[:-1])

salt  = np.interp(zr, zr0, salt0)
temp  = np.interp(zr, zr0, temp0)
rhoAm = gsw.rho(salt, temp, np.abs(zr))

# ------------ scalar inputs -------------------------------------------
N           = len(zr)
dx          = 600.
dy          = 300.
dt          = 100.
sg_runoff   = 200.
sg_typ      = 4
sg_dep      = -260.
sg_len      = 220.
sg_temp     = 0.
sg_salt     = 0.
sg_dye      = 1.

# ------------ write profiles ------------------------------------------
dye = np.zeros(N)
v   = np.zeros(N)
w   = np.zeros(N)

np.savetxt('../inputs/iceplume_zw.txt', np.array([zw]).T, fmt='%30.15e')
np.savetxt('../inputs/iceplume_zr.txt',
        np.array([temp, salt, v, w, dye, temp, salt]).T, fmt='%30.15e')
np.savetxt('../inputs/iceplume_scalar.txt',
        np.array([[N, dx, dy, dt, sg_typ, sg_dep, sg_len,
                   sg_runoff, sg_temp, sg_salt, sg_dye]]), fmt='%30.15e')

subprocess.call('cd ..;. ./build.bash', shell=True)
subprocess.call('cd ..; ./iceplume_test.exe 1', shell=True)

data_zr = np.loadtxt('../outputs/iceplume_zr.txt', skiprows=1)
data_zw = np.loadtxt('../outputs/iceplume_zw.txt', skiprows=1)
