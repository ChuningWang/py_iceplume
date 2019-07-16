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
ti = -1

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

# sE = 33.635662600733035
# sF = 0.00096084596605349154
# sP = 32.498979130286109
# tE = 1.1075875150694721
# tF = -3.4976593488690160
# tP = 0.95195337192371021

# sE = 33.688802520056356
# sF = 9.5996854244683632e-4
# sP = 32.550765856890976
# tE = 1.1474499585835232
# tF = -3.5760401655695260
# tP = 0.98788186124497179

sE = 33.647260203477778
sF = 9.6145565125687390e-4
sP = 32.506309460908142
tE = 1.0156293331810291
tF = -3.4431957764179657
tP = 0.86442994397291506

msk = detI==0
rhoDet = np.ma.masked_where(msk, rhoDet)
rhoP = rho[-1]*np.ones(zr.shape)
rhoP = np.ma.masked_where(msk, rhoP)

tDet = np.ma.masked_where(msk, tDet)
sDet = np.ma.masked_where(msk, sDet)

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

fig, axs = plt.subplots(2, 2,
        gridspec_kw={'width_ratios': [1, 4], 'height_ratios': [4, 1]})
fig.subplots_adjust(wspace=0.025, hspace=0.025)

axs[1, 0].set_xlabel(150*r' '+r'S [PSU]')
axs[1, 0].set_ylabel(100*r' '+r'T [$^{\circ}$C]')

for i in range(2):
    axs[0, i].set_ylim(-0.5, 3.5)
    axs[1, i].set_ylim(-4, -3)
    axs[0, i].set_yticks(np.arange(0, 4, 0.5))
    axs[1, i].set_yticks([-4, -3.5])
for i in range(2):
    axs[i, 0].set_xlim(-0.5, 0.5)
    axs[i, 1].set_xlim(30.5, 34.5)
    axs[i, 0].set_xticks([-0.5, 0])
    axs[i, 1].set_xticks(np.arange(31, 35, 0.5))
for i in range(2):
    axs[i, 0].spines['right'].set_visible(False)
    axs[i, 1].spines['left'].set_visible(False)
    axs[i, 0].yaxis.tick_left()
    axs[i, 1].yaxis.tick_right()
for i in range(2):
    axs[0, i].spines['bottom'].set_visible(False)
    axs[1, i].spines['top'].set_visible(False)
    axs[0, i].xaxis.tick_top()
    axs[1, i].xaxis.tick_bottom()
d1 = .0075
d2 = .015
kwargs = dict(transform=axs[0, 0].transAxes, color='k', clip_on=False)
axs[0, 0].plot((1-d1*4, 1+d1*4), (1-d2,   1+d2  ), **kwargs)
axs[0, 0].plot(( -d1*4,  +d1*4), ( -d2/4,  +d2/4), **kwargs)
kwargs.update(transform=axs[1, 0].transAxes)
axs[1, 0].plot((1-d1*4, 1+d1*4), ( -d2*4,  +d2*4), **kwargs)
axs[1, 0].plot(( -d1*4,  +d1*4), (1-d2  , 1+d2  ), **kwargs)
kwargs.update(transform=axs[0, 1].transAxes)
axs[0, 1].plot((1-d1  , 1+d1  ), ( -d2/4,  +d2/4), **kwargs)
axs[0, 1].plot(( -d1  ,  +d1  ), (1-d2  , 1+d2  ), **kwargs)
kwargs.update(transform=axs[1, 1].transAxes)
axs[1, 1].plot((1-d1  , 1+d1  ), (1-d2  , 1+d2  ), **kwargs)
axs[1, 1].plot(( -d1  ,  +d1  ), ( -d2*4,  +d2*4), **kwargs)

# for i in range(2):
#     for j in range(2):
#         axs[i, j].plot([0, sE], [tGade, tE], '--', color='gray', linewidth=0.5)
#         axs[i, j].plot([0, sE], [0, tE], '--', color='gray', linewidth=0.5)
#         axs[i, j].fill_between([0, sE], [tGade, tE], [0, tE], color='gray', alpha=0.1)

ct = axs[0, 1].contour(ss, tt, rr-1000, np.linspace(24, 28, 21),
        colors='grey', linewidths=0.5, linestyles='--')
axs[0, 1].clabel(ct, [24, 25, 26, 27], fontsize=6)

axs[1, 0].plot([sF, sE], [tF, tE], '--', color='gray')
axs[1, 0].plot(sF, tF, 'ok')
axs[1, 0].text(sF, tF-0.3, r'$S_F$, $T_F$', horizontalalignment='center')

axs[0, 1].plot([sF, sE], [tF, tE], '--', color='gray')
# axs[0, 1].plot([s1, s2], [t1, t2], 'r')
axs[0, 1].plot([sDet.min(), sDet.max()], [tDet.min(), tDet.max()], '-or', markersize=3)
axs[0, 1].plot(sE, tE, 'ok')
axs[0, 1].text(sE, tE-0.3, r'$S_E$, $T_E$', horizontalalignment='center')
axs[0, 1].plot(sP, tP, 'ok')
axs[0, 1].text(sP, tP-0.3, r'$S_P$, $T_P$', horizontalalignment='center')

sct = axs[0, 1].scatter(sAm, tAm, c=-zr, cmap=cm.amp, vmin=0, vmax=300)
cbar_ax = fig.add_axes([0.56, 0.83, 0.3, 0.02])
cbar = fig.colorbar(sct, cax=cbar_ax, orientation='horizontal', ticks=[0, 100, 200, 300])
cbar.set_label('Depth [m]')

axs[0, 1].contour(ss, tt, rr, [rhoDet1, rhoDet2], colors=[[0.0, 0.0, 0.5]], linestyles='--', linewidths=1)
# axs[0, 1].contour(ss, tt, rr, [rhoEnt1, rhoEnt2], colors=[[0.0, 0.5, 0.0]], linestyles='--', linewidths=1)

fig.savefig('../figs/KS-tsdiag_%03d.png' % int(time/3600), dpi=500)
plt.close()

# fig, axs = plt.subplots(1, 3, sharey=True)
# fig.subplots_adjust(wspace=0.1)
# axs[0].set_ylim(-35, 0)
# axs[0].set_xlim(25, 27)
# axs[1].set_xlim(31, 33)
# axs[2].set_xlim(0.5, 2.5)
# axs[0].set_xticks([25, 26, 27])
# axs[1].set_xticks([31, 32, 33])
# axs[2].set_xticks([1, 2])
# 
# axs[0].fill_between([24.5, 27.5], -33.5, -3.5, color=[0.9, 0.9, 0.9])
# axs[1].fill_between([31, 33],     -33.5, -3.5, color=[0.9, 0.9, 0.9])
# axs[2].fill_between([0.5, 2.5],   -33.5, -3.5, color=[0.9, 0.9, 0.9])
# 
# axs[0].plot(rhoP-1000,   zr, '--', color='grey', linewidth=2)
# axs[0].plot(rhoAm-1000,  zr, '-ok', markersize=3, linewidth=0.5)
# axs[0].plot(rho-1000,    zw, '-or', markersize=3, linewidth=0.5)
# axs[0].plot(rhoDet-1000, zr, '-og', markersize=3, linewidth=0.5)
# axs[1].plot(sAm,         zr, '-ok', markersize=3, linewidth=0.5)
# axs[1].plot(s,           zw, '-or', markersize=3, linewidth=0.5)
# axs[1].plot(sDet,        zr, '-og', markersize=3, linewidth=0.5)
# axs[2].plot(tAm,         zr, '-ok', markersize=3, linewidth=0.5)
# axs[2].plot(t,           zw, '-or', markersize=3, linewidth=0.5)
# axs[2].plot(tDet,        zr, '-og', markersize=3, linewidth=0.5)
# fig.savefig('../figs/KS-profile_%03d.png' % int(time/3600), dpi=500)
# plt.close()
