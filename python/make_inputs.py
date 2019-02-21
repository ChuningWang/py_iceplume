import numpy as np
from scipy import io

# ------------ vertical grid construction ------------------------------
zw = np.array([0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80, 90,
               100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300,
               320, 340, 360, 380, 400, 420, 440, 460, 480, 500, 550,
               600, 650, 700, 750, 800])

zw = -zw[::-1]
zr = 0.5*(zw[1:] + zw[:-1])

# ------------ scalar inputs -------------------------------------------
N = len(zr)
dx = 200.
dy = 200.
dt = 30.
sg_runoff = 200.
sg_temp = 0.
sg_salt = 0.
sg_dye01 = 1.
sg_typ = 3
sg_dep = 1
sg_len = 100.

# ------------ load profiles -------------------------------------------
data = io.loadmat('../data/carroll_2017_JGR_oceans_initialTS.mat')
depth_raw = np.array(data['depth']).squeeze()
temp_raw = np.array(data['temperature']).squeeze()
salt_raw = np.array(data['salinity']).squeeze()
depth_raw = depth_raw[::-1]
temp_raw = temp_raw[::-1]
salt_raw = salt_raw[::-1]

temp = np.interp(zr, depth_raw, temp_raw)
salt = np.interp(zr, depth_raw, salt_raw)
dye01 = np.zeros(N)
v = np.zeros(N)
w = np.zeros(N)

# temp = np.ones(N)*4
# salt = np.ones(N)*30

np.savetxt('../inputs/iceplume_zw.txt', np.array([zw]).T, fmt='%20.10e')
np.savetxt('../inputs/iceplume_zr.txt',
        np.array([temp, salt, v, w, dye01]).T, fmt='%20.10e')
np.savetxt('../inputs/iceplume_scalar.txt',
        np.array([[N, dx, dy, dt,
                   sg_runoff, sg_temp, sg_salt, sg_dye01,
                   sg_typ, sg_dep, sg_len]]),
        fmt='%03d ' + 9*'%20.10e ' + '%20.10e')
