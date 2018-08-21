import subprocess
import matplotlib.pyplot as plt
import numpy as np
from scipy import io
import csv
import pyroms

# ------------ basic grid parameters -----------------------------------
# vertical grid specs
theta_b = 1.0
theta_s = 1.0
Tcline = 10
N = 40

# ------------ vertical grid construction ------------------------------
h = np.ones((1, 1))*800
vgrd = pyroms.vgrid.s_coordinate_4(h, theta_b, theta_s, Tcline, N)
zr = vgrd.z_r[:]
zw = vgrd.z_w[:]

zw = np.array([0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 90,
               100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300,
               320, 340, 360, 380, 400, 420, 440, 460, 480, 500, 550,
               600, 650, 700, 750, 800])
# zw = np.array([0, 5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100,
#                120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
#                340, 360, 380, 400, 420, 440, 460, 480, 500, 550, 600,
#                650, 700, 750, 800, 850])
# zw = np.array([0, 5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100,
#                120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
#                340, 360, 380, 400, 420, 440, 460, 480, 500, 550, 600,
#                650, 700, 750, 800, 850, 900])

zw = -zw[::-1]
zr = 0.5*(zw[1:] + zw[:-1])
dz = np.diff(zw)
dy = 200

# ------------ load profiles -------------------------------------------
data = io.loadmat('../carroll_2017_JGR_oceans_initialTS.mat')
depth_raw = np.array(data['depth']).squeeze()
temp_raw = np.array(data['temperature']).squeeze()
salt_raw = np.array(data['salinity']).squeeze()
depth_raw = depth_raw[::-1]
temp_raw = temp_raw[::-1]
salt_raw = salt_raw[::-1]

temp = np.interp(zr, depth_raw, temp_raw)
salt = np.interp(zr, depth_raw, salt_raw)
dye01 = np.ones(N)
v = np.zeros(N)
w = np.zeros(N)

np.savetxt('../data/iceplume_zw.txt', zw)
np.savetxt('../data/iceplume_t.txt', temp)
np.savetxt('../data/iceplume_s.txt', salt)
np.savetxt('../data/iceplume_v.txt', v)
np.savetxt('../data/iceplume_w.txt', w)
np.savetxt('../data/iceplume_dye01.txt', dye01)

# ------------ build and run the executable ----------------------------
subprocess.call('../build.bash', shell=True)
subprocess.call('../iceplume_test.exe', shell=True)

# ------------ load results from txt files -----------------------------
f = open('../outputs/plume_out_zr.txt', 'rb')
header_zr = f.readline().split()
f.close()
data_zr = np.loadtxt('./outputs/plume_out_zr.txt', skiprows=1)
f = open('../outputs/plume_out_zw.txt', 'rb')
header_zw = f.readline().split()
f.close()
data_zw = np.loadtxt('../outputs/plume_out_zw.txt', skiprows=1)

data = {}
for (i, header) in enumerate(header_zr):
    data[header.decode('utf-8')] = data_zr[:, i]
for (i, header) in enumerate(header_zw):
    data[header.decode('utf-8')] = data_zw[:, i]

msk = data['s']==0
data['s'][msk] = np.NaN
data['t'][msk] = np.NaN
data['rho'][msk] = np.NaN

fig, axs = plt.subplots(1, 3)
axs[0].grid(True)
axs[1].grid(True)
axs[2].grid(True)
axs[0].set_ylim([-800, 0])
axs[1].set_ylim([-800, 0])
axs[2].set_ylim([-800, 0])
axs0 = axs[0].twiny()
axs[1].set_yticklabels([''])
axs[2].set_yticklabels([''])
axs[2].plot([0, 0], [-800, 0], '--', color='grey')
axs[0].set_xlabel(r'Salt [PSU]', color='b')
axs0.set_xlabel(r'Temp [$^{\circ}$C]', color='r')
axs[1].set_xlabel(r'$\rho$ Anomaly [kg$\cdot$m$^{-3}$]')
axs[2].set_xlabel(r'Ent/Det Vel [m$\cdot$s$^{-1}$]')

axs[0].plot(data['s'], zw, '-b')
axs[0].plot(data['sAm'], zr, '--.b')
axs0.plot(data['t'], zw, '-r')
axs0.plot(data['tAm'], zr, '--.r')
axs[1].plot(data['rho'], zw, '-k')
axs[1].plot(data['rhoAm'], zr, '--.k')
axs[2].plot(data['ent']/dz/dy, zr, '-.k')
axs[2].plot(data['det']/dz/dy, zr, '-.k')

axs[1].legend(['Plume', 'Ambient'])
plt.savefig('../figs/plume_out.png', dpi=600)
