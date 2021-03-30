import subprocess
import numpy as np
import matplotlib.pyplot as plt
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
dz  = np.diff(zw)

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
sg_trs_list = [5, 10, 20, 50, 100, 200, 300]
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

ent = np.zeros((N, len(sg_trs_list)))
det = np.zeros((N, len(sg_trs_list)))
rho = np.zeros((N+1, len(sg_trs_list)))
rhoAm = np.zeros((N, len(sg_trs_list)))

subprocess.call('cd ..;. ./build.bash', shell=True)

for i, sg_trs in enumerate(sg_trs_list):
    np.savetxt('../inputs/iceplume_zw.txt', np.array([zw]).T, fmt='%30.15e')
    np.savetxt('../inputs/iceplume_zr.txt',
            np.array([temp, salt, v, w, dye, temp, salt]).T, fmt='%30.15e')
    np.savetxt('../inputs/iceplume_scalar.txt',
            np.array([[N, dx, dy, dt, sg_typ, sg_dep, sg_len,
                       sg_trs, sg_temp, sg_salt, sg_dye]]), fmt='%30.15e')

    subprocess.call('cd ..; ./iceplume_test.exe 1', shell=True)
    data_zr = np.loadtxt('../outputs/iceplume_zr.txt', skiprows=1)
    data_zw = np.loadtxt('../outputs/iceplume_zw.txt', skiprows=1)
    ent[:, i] = data_zr[:, 2]
    det[:, i] = data_zr[:, 3]
    rho[:, i] = data_zw[:, 8]
    rhoAm[:, i] = data_zr[:, 8]

ent = np.ma.masked_where(ent==0, ent)
det = np.ma.masked_where(det==0, det)
rho = np.ma.masked_where(rho==0, rho)

entv = ent.copy()
detv = det.copy()
for i in range(N):
    entv[i] = ent[i]/sg_len/dz[i]
    detv[i] = det[i]/sg_len/dz[i]

rho = rho - 1000
rhoAm = rhoAm - 1000

clist = [u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd',
         u'#8c564b', u'#e377c2', u'#7f7f7f', u'#bcbd22', u'#17becf']

fig, ax = plt.subplots(1, 2, sharey=True, figsize=(10, 6))
fig.subplots_adjust(wspace=0.1, hspace=0.1,
        left=0.08, right=0.98, bottom=0.08, top=0.98)
ax[0].set_xlim([22, 29])
ax[1].set_xlim([-0.2, 0.7])
ax[0].set_ylim([-260, 0])
ax[0].set_yticks([-200, -150, -100, -50, 0])
ax[0].set_yticklabels(['200', '150', '100', '50', '0'])
ax[0].set_ylabel('Depth [m]')
ax[0].set_xlabel(r'Density [kg$\cdot$m$^{-3}$]')
ax[1].set_xlabel(r'Entrainment/Detrainment Vel [m$^3$s$^{-1}$]')
ax[0].plot(rhoAm[:, 0], zr, 'k', lw=3, label='Ambient')
ax[1].plot([0, 0], [-400, 0], '--k', lw=3)
ax[1].plot([0, 0], [-400, -400], '-', color=clist[0], label='Detrainment')
ax[1].plot([0, 0], [-400, -400], '--', color=clist[0], label='Entrainment')

for i, sg_trs in enumerate(sg_trs_list):
    ax[0].plot(rho[:, i], zw, color = clist[i], label=r'$Q_{sg}=$%3d' % sg_trs)
    ax[1].plot(detv[:, i], zr, '->', ms=3, color=clist[i])
    ax[1].plot(entv[:, i], zr, '--', ms=3, color=clist[i])
ax[0].legend()
ax[1].legend()

fig.savefig('../figs/ks_out.png', dpi=300)
