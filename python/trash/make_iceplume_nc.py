from datetime import datetime, timedelta
import subprocess
import csv
import numpy as np
from scipy import io
import netCDF4 as nc
import pyroms

# ------------ build the executable ------------------------------------
subprocess.call('cd ..; ./build.bash', shell=True)

# ------------ load data -----------------------------------------------
W = 2000.
Ds = 100.
m2amp = 0.00
s2amp = 0.00
dth = 50.
dthp = 0.025
runoff = 250.

islice0 = 131
islice1 = 141
jslice0 = 200
jslice1 = 300

dx0 = 200
dy0 = 200

fnum = 100
tnum = 24
tN = fnum*tnum

t0 = datetime(1900, 1, 1, 1) + timedelta(0)
dt = timedelta(1)

apd = '_w%05d_s%03d_r%04d_m2amp%4.2f_s2amp%4.2f_dth%03d_dthp%5.3f' \
        % (W, Ds, runoff, m2amp, s2amp, dth, dthp)
in_dir = '/glade/u/home/chuning/tidal_fjord/'
out_dir = '/glade/scratch/chuning/tmpdir_tidal_fjord/outputs' \
        + apd + '/outputs/1900/'

tlist = [(t0 + i*dt) for i in range(fnum)]
tlist[0] = tlist[0] - timedelta(hours=1)
fname_list = [tlist[i].strftime(out_dir + 'tidal_fjord' + apd \
        + '_his_%Y-%m-%dT%H:%M:%S.nc') for i in range(fnum)]

fname = fname_list[0]
gname = '/glade/u/home/chuning/tidal_fjord/tidal_fjord_grid' + apd + '.nc'

# load grid
grd = pyroms.grid.get_ROMS_grid(None, hist_file=fname, grid_file=gname)
N = grd.vgrid.N
eta, xi = grd.hgrid.mask_rho.shape

# initiate variables
time = np.zeros(tN)
zeta = np.zeros((tN, eta, xi))
salt = np.zeros((tN, N))
temp = np.zeros((tN, N))
dye01 = np.zeros((tN, N))
v = np.zeros((tN, N))
w = np.zeros((tN, N))

for i, fname in enumerate(fname_list):
    fin = nc.Dataset(fname, 'r')
    time[i*tnum:(i+1)*tnum] = fin.variables['ocean_time'][-tnum:]
    zeta[i*tnum:(i+1)*tnum] = fin.variables['zeta'][-tnum:]
    salt[i*tnum:(i+1)*tnum] = fin.variables['salt'][-tnum:, :, 136, 2]
    temp[i*tnum:(i+1)*tnum] = fin.variables['temp'][-tnum:, :, 136, 2]
    dye01[i*tnum:(i+1)*tnum] = fin.variables['dye_01'][-tnum:, :, 136, 2]
    v[i*tnum:(i+1)*tnum] = 0.5*(fin.variables['v'][-tnum:, :, 136, 2] +
                                fin.variables['v'][-tnum:, :, 135, 2])
    w[i*tnum:(i+1)*tnum] = 0.5*(fin.variables['w'][-tnum:, 1:, 136, 2] +
                                fin.variables['w'][-tnum:, :-1, 136, 2])
    fin.close()

grd = pyroms.grid.get_ROMS_grid(None, zeta=zeta,
        hist_file=fname, grid_file=gname)
zw = grd.vgrid.z_w[:, :, 136, 2]
zr = grd.vgrid.z_r[:, :, 136, 2]

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
    data[header] = np.zeros((tN, N))
for (i, header) in enumerate(header_zw):
    data[header] = np.zeros((tN, N+1))

print('Run the module')

for ti in range(tN):
    np.savetxt('../data/iceplume_zw.txt', zw[ti])
    np.savetxt('../data/iceplume_t.txt', temp[ti])
    np.savetxt('../data/iceplume_s.txt', salt[ti])
    np.savetxt('../data/iceplume_dye01.txt', dye01[ti])
    np.savetxt('../data/iceplume_v.txt', v[ti])
    np.savetxt('../data/iceplume_w.txt', w[ti])

    # ------------ run the executable --------------------------------------
    subprocess.call('cd ..; ./iceplume_test.exe', shell=True)

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
