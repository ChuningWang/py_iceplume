import subprocess
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

# ------------ load profiles -------------------------------------------
data = io.loadmat('./carroll_2017_JGR_oceans_initialTS.mat')
depth_raw = np.array(data['depth']).squeeze()
temp_raw = np.array(data['temperature']).squeeze()
salt_raw = np.array(data['salinity']).squeeze()
depth_raw = depth_raw[::-1]
temp_raw = temp_raw[::-1]
salt_raw = salt_raw[::-1]

temp = np.interp(zr, depth_raw, temp_raw)
salt = np.interp(zr, depth_raw, salt_raw)

np.savetxt('./data/iceplume_zw.txt', zw)
np.savetxt('./data/iceplume_t.txt', temp)
np.savetxt('./data/iceplume_s.txt', salt)

# ------------ build and run the executable ----------------------------
subprocess.call('./build.bash', shell=True)
subprocess.call('./iceplume_test.exe', shell=True)

# ------------ load results from txt files -----------------------------
f = open('./data/plume_out_zr.txt', 'rb')
header_zr = f.readline().split()
f.close()
data_zr = np.loadtxt('./data/plume_out_zr.txt', skiprows=1)
f = open('./data/plume_out_zw.txt', 'rb')
header_zw = f.readline().split()
f.close()
data_zw = np.loadtxt('./data/plume_out_zw.txt', skiprows=1)

data = {}
for (i, header) in enumerate(header_zr):
    data[header] = data_zr[:, i]
for (i, header) in enumerate(header_zw):
    data[header] = data_zw[:, i]
