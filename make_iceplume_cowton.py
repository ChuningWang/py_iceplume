import subprocess
import csv
import matplotlib.pyplot as plt
import numpy as np
from scipy import io
import csv
import pyroms

# ------------ vertical grid construction ------------------------------
dy = 200

fh = open('./data/outputs.txt', 'r')
csvR = csv.reader(fh, delimiter=' ', skipinitialspace=True)
data_cowton = []
for row in csvR:
    data_cowton.append([float(datai) for datai in row])
fh.close()
data_cowton = np.array(data_cowton)

depth_raw = data_cowton[::-1, 0]
temp_raw = data_cowton[::-1, 1]
salt_raw = data_cowton[::-1, 2]

zw = -np.arange(501)
zw = zw[::-1]
zr = 0.5*(zw[1:] + zw[:-1])
dz = np.diff(zw)

temp = np.interp(zr, depth_raw, temp_raw)
salt = np.interp(zr, depth_raw, salt_raw)
dye01 = np.zeros(temp.shape)
v = np.zeros(temp.shape)
w = np.zeros(temp.shape)

np.savetxt('./data/iceplume_zw.txt', zw)
np.savetxt('./data/iceplume_t.txt', temp)
np.savetxt('./data/iceplume_s.txt', salt)
np.savetxt('./data/iceplume_v.txt', v)
np.savetxt('./data/iceplume_w.txt', w)
np.savetxt('./data/iceplume_dye01.txt', dye01)

# ------------ build and run the executable ----------------------------
subprocess.call('./build.bash', shell=True)
subprocess.call('./iceplume_test.exe', shell=True)

# ------------ load results from txt files -----------------------------
f = open('./outputs/plume_out_zr.txt', 'rb')
header_zr = f.readline().split()
f.close()
data_zr = np.loadtxt('./outputs/plume_out_zr.txt', skiprows=1)
f = open('./outputs/plume_out_zw.txt', 'rb')
header_zw = f.readline().split()
f.close()
data_zw = np.loadtxt('./outputs/plume_out_zw.txt', skiprows=1)

data = {}
for (i, header) in enumerate(header_zr):
    data[header.decode('utf-8')] = data_zr[:, i]
for (i, header) in enumerate(header_zw):
    data[header.decode('utf-8')] = data_zw[:, i]

msk = data['s']==0
data['s'][msk] = np.NaN
data['t'][msk] = np.NaN
data['rho'][msk] = np.NaN
data['r'][msk] = np.NaN
data['r'][data['r'] == 0] = np.NaN

fig, axs = plt.subplots(2, 2)
fig.subplots_adjust(hspace=0.1, wspace=0.1)
for i in range(2):
    for j in range(2):
        axs[i, j].grid(True)
        axs[i, j].set_ylim([-500, 0])
        if i == 0:
            axs[i, j].xaxis.tick_top()
            axs[i, j].xaxis.set_label_position('top')
        if j == 1:
            axs[i, j].yaxis.tick_right()
            axs[i, j].yaxis.set_label_position('right')

axs[0, 0].set_ylabel('Depth')
axs[1, 0].set_ylabel('Depth')

axs[0, 0].set_xlabel(r'Salt [PSU]')
axs[0, 1].set_xlabel(r'Temp [$^{\circ}$C]')
axs[1, 0].set_xlabel(r'$\rho$ Anomaly [kg$\cdot$m$^{-3}$]')
axs[1, 1].set_xlabel(r'Radius [m]')

axs[0, 0].plot(data['s'], zw, '-C0', linewidth=4)
axs[0, 0].plot(data_cowton[:, 5], data_cowton[:, 0], '-C3')
axs[0, 0].plot(data_cowton[:, 2], data_cowton[:, 0], '-k')

axs[0, 1].plot(data['t'], zw, '-C0', linewidth=4)
axs[0, 1].plot(data_cowton[:, 4], data_cowton[:, 0], '-C3')
axs[0, 1].plot(data_cowton[:, 1], data_cowton[:, 0], '-k')

axs[1, 0].plot(data['rho'], zw, '-C0', linewidth=4)
axs[1, 0].plot(data_cowton[:, 6]-1000, data_cowton[:, 0], '-C3')
axs[1, 0].plot(data_cowton[:, 3]-1000, data_cowton[:, 0], '-k')

axs[1, 1].plot(data['r'], zw, '-C0', linewidth=4)
axs[1, 1].plot(data_cowton[:, 8], data_cowton[:, 0], '-C3')

axs[0, 0].legend(['ROMS', 'Cowton15', 'Ambient'])

plt.savefig('./figs/cowton15_out.png', dpi=600)

for i in range(2):
    for j in range(2):
        axs[i, j].set_ylim([-60, -40])
axs[0, 0].set_xlim([33.6, 34.0])
axs[0, 1].set_xlim([1.8, 2.1])
axs[1, 0].set_xlim([27.15, 27.35])
axs[1, 1].set_xlim([62, 66])

axs[0, 0].plot(data['s'], zw, 'oC0', linewidth=4, markerfacecolor='None', markersize=5)
axs[0, 0].plot(data_cowton[:, 5], data_cowton[:, 0], 'oC3', markerfacecolor='None', markersize=3)
axs[0, 0].plot(data_cowton[:, 2], data_cowton[:, 0], 'ok', markerfacecolor='None', markersize=3)

axs[0, 1].plot(data['t'], zw, 'oC0', linewidth=4, markerfacecolor='None', markersize=5)
axs[0, 1].plot(data_cowton[:, 4], data_cowton[:, 0], 'oC3', markerfacecolor='None', markersize=3)
axs[0, 1].plot(data_cowton[:, 1], data_cowton[:, 0], 'ok', markerfacecolor='None', markersize=3)

axs[1, 0].plot(data['rho'], zw, 'oC0', linewidth=4, markerfacecolor='None', markersize=5)
axs[1, 0].plot(data_cowton[:, 6]-1000, data_cowton[:, 0], 'oC3', markerfacecolor='None', markersize=3)
axs[1, 0].plot(data_cowton[:, 3]-1000, data_cowton[:, 0], 'ok', markerfacecolor='None', markersize=3)

axs[1, 1].plot(data['r'], zw, 'oC0', linewidth=4, markerfacecolor='None', markersize=5)
axs[1, 1].plot(data_cowton[:, 8], data_cowton[:, 0], 'oC3', markerfacecolor='None', markersize=3)

plt.savefig('./figs/cowton15_out2.png', dpi=600)
