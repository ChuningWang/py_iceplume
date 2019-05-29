import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from cmocean import cm
import gsw

tGade = -90

data_zr = np.loadtxt('../outputs/iceplume_zr.txt', skiprows=1)
data_zw = np.loadtxt('../outputs/iceplume_zw.txt', skiprows=1)

zr = data_zr[:, 1]
zw = data_zw[:, 1]
detI = data_zr[:, 4]

ent = data_zr[:, 2]
det = data_zr[:, 3]
detF = data_zr[:, 9]
detE = data_zr[:, 10]

tAm = data_zr[:, 5]
sAm = data_zr[:, 6]
t = data_zw[:, 5]
s = data_zw[:, 6]

tF = -8.79e-3
tE = 1.023
sF = 0
sE = 33.6106

msk = detI==0
rho   = data_zw[:, 8]
rhoAm = data_zr[:, 8]
rhoDet = data_zr[:, 13]
rhoDet = np.ma.masked_where(msk, rhoDet)
rhoP = rho[-1]*np.ones(zr.shape)
rhoP = np.ma.masked_where(msk, rhoP)

tDet = data_zr[:, 11]
sDet = data_zr[:, 12]
tDet = np.ma.masked_where(msk, tDet)
sDet = np.ma.masked_where(msk, sDet)

rhoAm0 = gsw.rho(sAm, tAm, 0)
rhoDet1 = rhoAm0[~msk].min()
rhoDet2 = rhoAm0[~msk].max()
msk2 = ent<0
rhoEnt1 = rhoAm0[msk2].min()
rhoEnt2 = rhoAm0[msk2].max()

s0 = np.arange(0, 35, 0.1)
t0 = np.arange(-1, 5, 0.1)

ss, tt = np.meshgrid(s0, t0)
rr = gsw.rho(ss, tt, 0)

fig, axs = plt.subplots(1, 2, sharey=True)
fig.subplots_adjust(wspace=0.05)

axs[0].set_xlabel(r'T [$^{\circ}$C]')
axs[0].set_ylabel(r'S [PSU]')

axs[0].set_ylim(-0.5, 3.5)
axs[0].set_xlim(-0.5, 3)
axs[1].set_xlim(31, 34.5)
axs[0].spines['right'].set_visible(False)
axs[1].spines['left'].set_visible(False)
axs[0].yaxis.tick_left()
axs[1].yaxis.tick_right()
d = .015  # how big to make the diagonal lines in axes coordinates
kwargs = dict(transform=axs[0].transAxes, color='k', clip_on=False)
axs[0].plot((1-d, 1+d), (1-d, 1+d), **kwargs)  # top-left diagonal
axs[0].plot((1-d, 1+d), ( -d,  +d), **kwargs)  # top-right diagonal
kwargs.update(transform=axs[1].transAxes)  # switch to the bottom axes
axs[1].plot(( -d,  +d), (1-d, 1+d), **kwargs)  # bottom-left diagonal
axs[1].plot(( -d,  +d), ( -d,  +d), **kwargs)  # bottom-right diagonal

ct = axs[0].contour(ss, tt, rr-1000, np.linspace(0, 30, 61),
        colors='grey', linewidths=1)
axs[0].clabel(ct)
axs[0].plot(sF, tF, 'ok')
axs[0].plot([sF, sE], [tF, tE], '--', color='gray')
axs[0].text(sF, tF-0.2, r'$S_F$, $T_F$', horizontalalignment='center')

ct = axs[1].contour(ss, tt, rr-1000, np.linspace(0, 30, 61),
        colors='grey', linewidths=1)
axs[1].clabel(ct)
sct = axs[1].scatter(sAm, tAm, c=-zr, cmap=cm.amp, vmin=0, vmax=300)
axs[1].plot(s[-1], t[-1], 'ok')
axs[1].plot(sE, tE, 'ok')
axs[1].plot(sF, tF, 'ok')
axs[1].plot([sF, sE], [tF, tE], '--', color='gray')
axs[1].text(sE, tE-0.2, r'$S_E$, $T_E$', horizontalalignment='center')
axs[1].text(s[-1], t[-1]-0.2, r'$S_P$, $T_P$', horizontalalignment='center')

# cmap_array = np.array([[0.0, 0.0, 0.0, 0.0],
#                        [0.0, 0.0, 0.5, 0.2],
#                        [0.0, 0.0, 0.0, 0.0]])
# cmap = colors.ListedColormap(cmap_array)
# axs[1].contourf(ss, tt, rr, [rhoDet1, rhoDet2], cmap=cmap, linestyles='--')
# cmap_array = np.array([[0.0, 0.0, 0.0, 0.0],
#                        [0.0, 0.5, 0.0, 0.2],
#                        [0.0, 0.0, 0.0, 0.0]])
# cmap = colors.ListedColormap(cmap_array)
# axs[1].contourf(ss, tt, rr, [rhoEnt1, rhoEnt2], cmap=cmap, linestyles='--')

axs[1].contour(ss, tt, rr, [rhoDet1, rhoDet2], colors=[[0.0, 0.0, 0.5]], linestyles='--')
axs[1].contour(ss, tt, rr, [rhoEnt1, rhoEnt2], colors=[[0.0, 0.5, 0.0]], linestyles='--')

cbar_ax = fig.add_axes([0.59, 0.85, 0.3, 0.02])
cbar = fig.colorbar(sct, cax=cbar_ax, orientation='horizontal', ticks=[0, 100, 200, 300])
cbar.set_label('Depth [m]')

axs[1].plot([32.94975, 31.33634], [1.00272, 0.95318], 'r')
fig.savefig('../figs/KS-tsdiag.png', dpi=500)
plt.close()

fig, axs = plt.subplots(1, 3, sharey=True)
fig.subplots_adjust(wspace=0.1)
axs[0].set_ylim(-35, 0)
axs[0].set_xlim(25, 27)
axs[1].set_xlim(31, 33)
axs[2].set_xlim(0.5, 2.5)
axs[0].set_xticks([25, 26, 27])
axs[1].set_xticks([31, 32, 33])
axs[2].set_xticks([1, 2])

axs[0].fill_between([24.5, 27.5], -33.5, -3.5, color=[0.9, 0.9, 0.9])
axs[1].fill_between([31, 33],     -33.5, -3.5, color=[0.9, 0.9, 0.9])
axs[2].fill_between([0.5, 2.5],   -33.5, -3.5, color=[0.9, 0.9, 0.9])

axs[0].plot(rhoP-1000,   zr, '--', color='grey', linewidth=2)
axs[0].plot(rhoAm-1000,  zr, '-ok', markersize=3, linewidth=0.5)
axs[0].plot(rho-1000,    zw, '-or', markersize=3, linewidth=0.5)
axs[0].plot(rhoDet-1000, zr, '-og', markersize=3, linewidth=0.5)
axs[1].plot(sAm,         zr, '-ok', markersize=3, linewidth=0.5)
axs[1].plot(s,           zw, '-or', markersize=3, linewidth=0.5)
axs[1].plot(sDet,        zr, '-og', markersize=3, linewidth=0.5)
axs[2].plot(tAm,         zr, '-ok', markersize=3, linewidth=0.5)
axs[2].plot(t,           zw, '-or', markersize=3, linewidth=0.5)
axs[2].plot(tDet,        zr, '-og', markersize=3, linewidth=0.5)
fig.savefig('../figs/KS-profile.png', dpi=500)
