import subprocess
import numpy as np
from scipy import io
import csv
import pyroms

# ------------ build the executable ------------------------------------
subprocess.call('../build.bash', shell=True)

# ------------ basic grid parameters -----------------------------------
# vertical grid specs
theta_b = 1.0
theta_s = 1.0
Tcline = 10
N = 40

# get header info
f = open('../data/plume_out_zr.txt', 'rb')
header_zr = f.readline().split()
f.close()
f = open('../data/plume_out_zw.txt', 'rb')
header_zw = f.readline().split()
f.close()

for i in range(len(header_zr)):
    header_zr[i] = header_zr[i].decode('utf-8')
for i in range(len(header_zw)):
    header_zw[i] = header_zw[i].decode('utf-8')

# ------------ vertical grid construction ------------------------------
h = np.ones((1, 1))*800
vgrd = pyroms.vgrid.s_coordinate_4(h, theta_b, theta_s, Tcline, N)
zr = vgrd.z_r[:]
zw = vgrd.z_w[:]

# ------------ vertical grid from Carroll 2017 -------------------------
zw2 = np.array([0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 90, 100,
                120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
                340, 360, 380, 400, 420, 440, 460, 480, 500, 550, 600,
                650, 700, 750, 800])
zw2 = -zw2[::-1]
zr2 = 0.5*(zw2[1:] + zw2[:-1])

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

temp2 = np.interp(zr2, depth_raw, temp_raw)
salt2 = np.interp(zr2, depth_raw, salt_raw)

# ------------ run executable ------------------------------------------
# first, with ROMS grid
np.savetxt('../data/iceplume_zw.txt', zw)
np.savetxt('../data/iceplume_t.txt', temp)
np.savetxt('../data/iceplume_s.txt', salt)

subprocess.call('../iceplume_test.exe', shell=True)

data_zr = np.loadtxt('../data/plume_out_zr.txt', skiprows=1)
data_zw = np.loadtxt('../data/plume_out_zw.txt', skiprows=1)

data = {}

for (i, header) in enumerate(header_zr):
    if str(header) != 'lev':
        data[header] = data_zr[:, i]
for (i, header) in enumerate(header_zw):
    data[header] = data_zw[:, i]

# second, with Carroll 2017 grid
np.savetxt('../data/iceplume_zw.txt', zw2)
np.savetxt('../data/iceplume_t.txt', temp2)
np.savetxt('../data/iceplume_s.txt', salt2)

subprocess.call('../iceplume_test.exe', shell=True)

data_zr = np.loadtxt('../data/plume_out_zr.txt', skiprows=1)
data_zw = np.loadtxt('../data/plume_out_zw.txt', skiprows=1)

data2 = {}

for (i, header) in enumerate(header_zr):
    if str(header) != 'lev':
        data2[header] = data_zr[:, i]
for (i, header) in enumerate(header_zw):
    data2[header] = data_zw[:, i]
