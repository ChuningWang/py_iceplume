import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors
from cmocean import cm
import gsw
import pickle

matplotlib.rcParams.update({'font.size': 8})

data = pickle.load(open('./py_iceplume.pickle', 'rb'))

tGade = -90
ti = 0

time = data['time'][ti]
zr = data['zr'][ti, :, 0]
zw = data['zw'][ti, :, 0]
detI = data['detI'][ti, :, 0]
ent = data['ent'][ti, :, 0]
det = data['det'][ti, :, 0]
detF = data['detF'][ti, :, 0]
detE = data['detE'][ti, :, 0]
tAm = data['tAm'][ti, :, 0]
sAm = data['sAm'][ti, :, 0]
rhoAm = data['rhoAm'][ti, :, 0]
t = data['t'][ti, :, 0]
s = data['s'][ti, :, 0]
rho = data['rho'][ti, :, 0]
tDet = data['tDet'][ti, :, 0]
sDet = data['sDet'][ti, :, 0]
rhoDet = data['rhoDet'][ti, :, 0]

msk = detI==0
h1 = zr[~msk].max()
h2 = zr[~msk].min()

rhoDet = np.ma.masked_where(msk, rhoDet)
tDet = np.ma.masked_where(msk, tDet)
sDet = np.ma.masked_where(msk, sDet)

rhoP = rho[-1]

rhoAm0 = gsw.rho(sAm, tAm, 0)

rhoDet1 = rhoAm0[~msk].min()
rhoDet2 = rhoAm0[~msk].max()

msk2 = ent<0
rhoEnt1 = rhoAm0[msk2].min()
rhoEnt2 = rhoAm0[msk2].max()

s0 = np.arange(0, 35, 0.1)
t0 = np.arange(-5, 5, 0.1)

ss, tt = np.meshgrid(s0, t0)
rr = gsw.rho(ss, tt, 0)

fig, axs = plt.subplots(1, 3, sharey=True)
fig.subplots_adjust(wspace=0.1)
axs[0].set_ylim(-60, 0)
axs[0].set_xlim(24.5, 27.5)
axs[1].set_xlim(30.5, 33.5)
axs[2].set_xlim(0.5, 2.5)
axs[0].set_xticks([25, 26, 27])
axs[1].set_xticks([31, 32, 33])
axs[2].set_xticks([1, 2])
axs[0].set_ylabel('Depth [m]')
axs[0].set_xlabel(r'$\sigma$ [kg/m$^3$]')
axs[1].set_xlabel(r'S [PSU]')
axs[2].set_xlabel(r'T [$^{\circ}$C]')
axs[0].text(24.6, -58, 'a)')
axs[1].text(30.6, -58, 'b)')
axs[2].text(0.6,  -58, 'c)')

axs[0].fill_between([0, 100], h2, h1, color=[0.9, 0.9, 0.9], alpha=0.7)
axs[1].fill_between([0, 100], h2, h1, color=[0.9, 0.9, 0.9], alpha=0.7)
axs[2].fill_between([0, 100], h2, h1, color=[0.9, 0.9, 0.9], alpha=0.7)

# axs[0].plot(rhoP-1000,   zr, '--', color='grey', linewidth=1)
axs[0].plot(rhoAm-1000,  zr, '-k', linewidth=1)
axs[0].plot(rho-1000,    zw, '-r', linewidth=1)
axs[0].plot(rhoDet-1000, zr, '-og', markersize=3, linewidth=0.5)
axs[1].plot(sAm,         zr, '-k', linewidth=1)
axs[1].plot(s,           zw, '-r', linewidth=1)
axs[1].plot(sDet,        zr, '-og', markersize=3, linewidth=0.5)
axs[2].plot(tAm,         zr, '-k', linewidth=1)
axs[2].plot(t,           zw, '-r', linewidth=1)
axs[2].plot(tDet,        zr, '-og', markersize=3, linewidth=0.5)

axs[2].legend(['Ambient', 'Plume', 'Detrainment'])
fig.savefig('../figs/KS-profile_%03d.png' % int(time/3600), dpi=500)
plt.close()
