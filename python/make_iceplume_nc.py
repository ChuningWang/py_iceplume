import subprocess
import csv
import numpy as np
from scipy import io
import netCDF4 as nc
import pyroms

# ------------ build the executable ------------------------------------
subprocess.call('../build.bash', shell=True)

# ------------ load data -----------------------------------------------
fin = nc.Dataset('../tidal_fjord/tidal_fjord_w02000_s100_r0250_m2amp0.00_s2amp0.00_dth050_dthp0.025_his_00061.nc', 'r')
time = fin.variables['ocean_time'][:]
zeta = fin.variables['zeta'][:]
salt = fin.variables['salt'][:, :, 136, 2]
temp = fin.variables['temp'][:, :, 136, 2]
dye01 = fin.variables['dye_01'][:, :, 136, 2]
v = 0.5*(fin.variables['v'][:, :, 137, 2] +
         fin.variables['v'][:, :, 136, 2])
w = 0.5*(fin.variables['w'][:, 1:, 136, 2] +
         fin.variables['w'][:, :-1, 136, 2])
# fin.close()

grd = pyroms.grid.get_ROMS_grid(None, zeta=zeta,
        hist_file='../tidal_fjord/tidal_fjord_w02000_s100_r0250_m2amp0.00_s2amp0.00_dth050_dthp0.025_his_00061.nc',
        grid_file='../tidal_fjord/tidal_fjord_grid_w02000_s100_r0250_m2amp0.00_s2amp0.00_dth050_dthp0.025.nc')
zw = grd.vgrid.z_w[:, :, 136, 2]
zr = grd.vgrid.z_r[:, :, 136, 2]

nt, nw = zw.shape
nr = nw - 1

# ------------ initiate variables --------------------------------------
f = open('../outputs/plume_out_zr.txt', 'rb')
header_zr = f.readline().split()
f.close()
f = open('../outputs/plume_out_zw.txt', 'rb')
header_zw = f.readline().split()
f.close()

for i in range(len(header_zr)):
    header_zr[i] = header_zr[i].decode('utf-8')
for i in range(len(header_zw)):
    header_zw[i] = header_zw[i].decode('utf-8')

data = {}
for (i, header) in enumerate(header_zr):
    data[header] = np.zeros((nt, nr))
for (i, header) in enumerate(header_zw):
    data[header] = np.zeros((nt, nw))


print('Run the module')

for ti in range(nt):
    np.savetxt('../data/iceplume_zw.txt', zw[ti])
    np.savetxt('../data/iceplume_t.txt', temp[ti])
    np.savetxt('../data/iceplume_s.txt', salt[ti])
    np.savetxt('../data/iceplume_dye01.txt', dye01[ti])
    np.savetxt('../data/iceplume_v.txt', v[ti])
    np.savetxt('../data/iceplume_w.txt', w[ti])

    # ------------ run the executable --------------------------------------
    subprocess.call('../iceplume_test.exe', shell=True)

    # ------------ load results from txt files -----------------------------
    data_zr = np.loadtxt('../outputs/plume_out_zr.txt', skiprows=1)
    data_zw = np.loadtxt('../outputs/plume_out_zw.txt', skiprows=1)

    for (i, header) in enumerate(header_zr):
        if str(header) != 'lev':
            data[header][ti, :] = data_zr[:, i]
    for (i, header) in enumerate(header_zw):
        data[header][ti, :] = data_zw[:, i]

det = data['det']
ent = data['ent']
det[det==0] = np.NaN
ent[ent==0] = np.NaN
