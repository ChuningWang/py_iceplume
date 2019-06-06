import pyroms
import numpy as np
import pickle
import netCDF4 as nc
import matplotlib.pyplot as plt

hist_file = './fjord_his.nc'
grid_file = '/glade/work/chuning/roms_archive/fjord_ks_full/fjord_grid.nc'
river_file = '/glade/work/chuning/roms_archive/fjord_ks_full/fjord_river.nc'
pick_file = './py_iceplume.pickle'

fh = nc.Dataset(river_file)
xloc = fh.variables['subglacial_Erange'][:][:, 1]
yloc = fh.variables['subglacial_Xrange'][:][:, 1]
fh.close()

grd = pyroms.grid.get_ROMS_grid('grd1',
    hist_file=hist_file,
    grid_file=grid_file)
zr_raw = grd.vgrid.z_r[:][:, xloc[0]:xloc[1]+1, yloc[0]:yloc[1]+1]
zw_raw = grd.vgrid.z_w[:][:, xloc[0]:xloc[1]+1, yloc[0]:yloc[1]+1]
zr = zr_raw.mean(axis=(1, 2))
zw = zw_raw.mean(axis=(1, 2))
dz = np.diff(zw)
dy = 300
nt = 11

fh = nc.Dataset(hist_file)
u0 = fh.variables['u'][:, :, 31, 1]
fh.close()

data = pickle.load(open(pick_file, 'rb'))
det = data['det'][:, :, 0]
ent = data['ent'][:, :, 0]
u = (det+ent)/np.tile(dz, (nt, 1))/dy
