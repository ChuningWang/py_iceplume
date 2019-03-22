
# coding: utf-8

# ### Standalone verison of ICEPLUME using ROMS output files.

# In[1]:


# get_ipython().run_line_magic('matplotlib', 'inline')

# -------------- import modules --------------------------------
from datetime import datetime, timedelta
import subprocess
import pickle

import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

import pyroms
from cmocean import cm

plt.rcParams['figure.figsize'] = [18, 16]
plt.rcParams['font.size'] = 20


# In[2]:


# -------------- functionals --------------------------------
def get_zr(zeta, h, vgrid):
    """ get z at rho points from grid and zeta info. """

    ti = zeta.shape[0]
    zr = np.empty((ti, vgrid.N) + h.shape, 'd')
    if vgrid.Vtrans == 1:
        for k in range(vgrid.N):
            z0 = vgrid.hc * vgrid.s_rho[k] + (h - vgrid.hc) * vgrid.Cs_r[k]
            zr[:, k, :] = z0 + zeta * (1.0 + z0 / h)
    elif vgrid.Vtrans == 2 or vgrid.Vtrans == 4 or vgrid.Vtrans == 5:
        for k in range(vgrid.N):
            z0 = (vgrid.hc * vgrid.s_rho[k] + h * vgrid.Cs_r[k]) / (vgrid.hc + h)
            zr[:, k, :] = zeta + (zeta + h) * z0

    return zr


def get_zw(zeta, h, vgrid):
    """ get z at rho points from grid and zeta info. """

    ti = zeta.shape[0]
    zw = np.empty((ti, vgrid.Np) + h.shape, 'd')
    if vgrid.Vtrans == 1:
        for k in range(vgrid.Np):
            z0 = vgrid.hc * vgrid.s_w[k] + (h - vgrid.hc) * vgrid.Cs_w[k]
            zw[:, k, :] = z0 + zeta * (1.0 + z0 / h)
    elif vgrid.Vtrans == 2 or vgrid.Vtrans == 4 or vgrid.Vtrans == 5:
        for k in range(vgrid.Np):
            z0 = (vgrid.hc * vgrid.s_w[k] + h * vgrid.Cs_w[k]) / (vgrid.hc + h)
            zw[:, k, :] = zeta + (zeta + h) * z0

    return zw


def make_input(roms_his_file, roms_river_file, roms_grid_file, tracer1d=True):
    """ make input files for py_iceplume from roms outputs. """
    
    fh = nc.Dataset(roms_his_file, 'r')
    fh_river = nc.Dataset(roms_river_file, 'r')
    grd = pyroms.grid.get_ROMS_grid('grd_temp', hist_file=roms_his_file,
                                    grid_file=roms_grid_file)
    ntracer = len(fh.dimensions['tracer'])-2
    N = len(fh.dimensions['N'])

    epos = fh_river.variables['river_Eposition'][:]
    xpos = fh_river.variables['river_Xposition'][:]
    rdir = fh_river.variables['river_direction'][:]
    sgdep = fh_river.variables['subglacial_depth'][:]
    sgtyp = fh_river.variables['subglacial_type'][:]
    sglen = fh_river.variables['subglacial_length'][:]

    river_time = fh_river.variables['river_time'][:]
    sgtrs_raw = fh_river.variables['subglacial_transport'][:]
    sgtemp_raw = fh_river.variables['subglacial_temp'][:]
    sgsalt_raw = fh_river.variables['subglacial_salt'][:]
    sgdye_raw = []
    for j in range(ntracer):
        sgdye_raw.append(fh_river.variables['subglacial_dye_%02d' % (j+1)][:])
    fh_river.close()
    sgdye_raw = np.array(sgdye_raw)

    time = fh.variables['ocean_time'][:]
    dt = np.diff(time).mean()

    ntime = len(time)
    nriver = len(epos)

    sgtrs = np.ma.zeros((ntime, nriver))
    sgtemp = np.ma.zeros((ntime, nriver))
    sgsalt = np.ma.zeros((ntime, nriver))
    sgdye = np.ma.zeros((ntime, ntracer, nriver))

    mask_rho = grd.hgrid.mask_rho
    epos_rho = np.zeros(nriver).astype('int')
    xpos_rho = np.zeros(nriver).astype('int')
    rdir2 = np.zeros(nriver)
    
    if tracer1d:
        for i in range(nriver):
            sgtrs[:, i] = np.interp(time, river_time, sgtrs_raw[:, i])
            sgtemp[:, i] = np.interp(time, river_time, sgtemp_raw)
            sgsalt[:, i] = np.interp(time, river_time, sgsalt_raw)
            for j in range(ntracer):
                sgdye[:, j, i] = np.interp(time, river_time, sgdye_raw[j, :])
    else:
        for i in range(nriver):
            sgtrs[:, i] = np.interp(time, river_time, sgtrs_raw[:, i])
            sgtemp[:, i] = np.interp(time, river_time, sgtemp_raw[:, i])
            sgsalt[:, i] = np.interp(time, river_time, sgsalt_raw[:, i])
            for j in range(ntracer):
                sgdye[:, j, i] = np.interp(time, river_time, sgdye_raw[j, :, i])
    
    for i in range(len(epos)):
        if rdir[i] == 0:
            if ((mask_rho[epos[i], xpos[i]] == 0) & (mask_rho[epos[i]-1, xpos[i]] == 1)):
                epos_rho[i] = int(epos[i]-1)
                xpos_rho[i] = int(xpos[i])
                rdir2[i] = -1
            if ((mask_rho[epos[i], xpos[i]] == 1) & (mask_rho[epos[i]-1, xpos[i]] == 0)):
                epos_rho[i] = int(epos[i])
                xpos_rho[i] = int(xpos[i])
                rdir2[i] = 1
            else:
                epos_rho[i] = int(epos[i])
                xpos_rho[i] = int(xpos[i])
                rdir2[i] = 0
        if rdir[i] == 1:
            if ((mask_rho[epos[i], xpos[i]] == 0) & (mask_rho[epos[i], xpos[i]-1] == 1)):
                epos_rho[i] = int(epos[i])
                xpos_rho[i] = int(xpos[i]-1)
                rdir2[i] = -1
            if ((mask_rho[epos[i], xpos[i]] == 1) & (mask_rho[epos[i], xpos[i]-1] == 0)):
                epos_rho[i] = int(epos[i])
                xpos_rho[i] = int(xpos[i])
                rdir2[i] = 1
            else:
                epos_rho[i] = int(epos[i])
                xpos_rho[i] = int(xpos[i])
                rdir2[i] = 0

    h = np.ma.zeros(nriver)
    dx = np.ma.zeros(nriver)
    dy = np.ma.zeros(nriver)
    for i in range(nriver):
        h[i] = grd.vgrid.h[epos_rho[i], xpos_rho[i]]
        dx[i] = grd.hgrid.dx[epos_rho[i], xpos_rho[i]]
        dy[i] = grd.hgrid.dy[epos_rho[i], xpos_rho[i]]

    zeta = np.ma.zeros((ntime, nriver))
    salt = np.ma.zeros((ntime, N, nriver))
    temp = np.ma.zeros((ntime, N, nriver))
    v = np.ma.zeros((ntime, N, nriver))
    w = np.ma.zeros((ntime, N, nriver))
    if ntracer>0:
        dye = np.ma.zeros((ntime, ntracer, N, nriver))

    for i in range(nriver):
        zeta[:, i] = fh.variables['zeta'][:, epos_rho[i], xpos_rho[i]]
        salt[:, :, i] = fh.variables['salt'][:, :, epos_rho[i], xpos_rho[i]]
        temp[:, :, i] = fh.variables['temp'][:, :, epos_rho[i], xpos_rho[i]]
        for j in range(ntracer):
            dye[:, j, :, i] = fh.variables['dye_%02d' % (j+1)][:, :, epos_rho[i], xpos_rho[i]]
        w[:, :, i] = 0.5*(fh.variables['w'][:, 1:, epos_rho[i], xpos_rho[i]]+
                          fh.variables['w'][:, :-1, epos_rho[i], xpos_rho[i]])
        if rdir[i] == 0:
            v[:, :, i] = 0.5*(fh.variables['v'][:, :, epos_rho[i]-1, xpos_rho[i]]+
                              fh.variables['v'][:, :, epos_rho[i], xpos_rho[i]])
        elif rdir[i] == 1:
            v[:, :, i] = 0.5*(fh.variables['u'][:, :, epos_rho[i], xpos_rho[i]-1]+
                              fh.variables['u'][:, :, epos_rho[i], xpos_rho[i]])
    fh.close()

    zr = get_zr(zeta, h, grd.vgrid)
    zw = get_zw(zeta, h, grd.vgrid)
    
    roms_input = {}
    roms_input['N'] = N
    roms_input['ntime'] = ntime
    roms_input['nriver'] = nriver
    roms_input['ntracer'] = ntracer
    roms_input['dt'] = dt
    
    roms_input['dx'] = dx
    roms_input['dy'] = dy
    roms_input['zr'] = zr
    roms_input['zw'] = zw
    
    roms_input['epos'] = epos
    roms_input['xpos'] = xpos
    roms_input['rdir'] = rdir
    roms_input['sgdep'] = sgdep
    roms_input['sgtyp'] = sgtyp
    roms_input['sglen'] = sglen

    roms_input['sgtrs'] = sgtrs
    roms_input['sgtemp'] = sgtemp
    roms_input['sgsalt'] = sgsalt
    roms_input['sgdye'] = sgdye
    
    roms_input['time'] = time
    roms_input['salt'] = salt
    roms_input['temp'] = temp
    roms_input['dye'] = dye
    roms_input['w'] = w
    roms_input['v'] = v
    
    return roms_input


# In[7]:


# -------------- generate inputs --------------------------------
hist_file = '/Users/cw686/roms_archive/fjord_test2/outputs_archive/fjord_his_0.1s.nc'
grid_file = '/Users/cw686/roms_archive/fjord_test2/fjord_grid.nc'
river_file = '/Users/cw686/roms_archive/fjord_test2/fjord_river.nc'
roms_input = make_input(hist_file, river_file, grid_file)


# In[8]:


# ------------ build the executable ------------------------------------
subprocess.call('cd ..; ./build.bash', shell=True)


# In[9]:


# ------------ build input file and run --------------------------------------
trange = range(roms_input['ntime'])
# irange = range(roms_input['nriver'])
irange = range(1, 2)

N = roms_input['N']
ntime = len(trange)
nriver = len(irange)
ntracer = roms_input['ntracer']

iceplume_out = {}
iceplume_out['time'] = roms_input['time'][trange]
iceplume_out['epos'] = roms_input['epos'][irange]
iceplume_out['xpos'] = roms_input['xpos'][irange]
iceplume_out['zw'] = np.zeros((ntime, N+1, nriver))
iceplume_out['f'] = np.zeros((ntime, N+1, nriver))
iceplume_out['w'] = np.zeros((ntime, N+1, nriver))
iceplume_out['t'] = np.zeros((ntime, N+1, nriver))
iceplume_out['s'] = np.zeros((ntime, N+1, nriver))
iceplume_out['a'] = np.zeros((ntime, N+1, nriver))

iceplume_out['zr'] = np.zeros((ntime, N, nriver))
iceplume_out['ent'] = np.zeros((ntime, N, nriver))
iceplume_out['det'] = np.zeros((ntime, N, nriver))
iceplume_out['m'] = np.zeros((ntime, N, nriver))

iceplume_out['dye'] = np.zeros((ntime, ntracer, nriver))

for i, ti in enumerate(trange):
    for j, ri in enumerate(irange):
        np.savetxt('../inputs/iceplume_zw.txt',
                np.array([roms_input['zw'][ti, :, ri]]).T, fmt='%20.10e')
        np.savetxt('../inputs/iceplume_zr.txt',
                np.vstack((roms_input['temp'][ti, :, ri],
                           roms_input['salt'][ti, :, ri],
                           roms_input['v'][ti, :, ri],
                           roms_input['w'][ti, :, ri],
                           roms_input['dye'][ti, :, :, ri]
                          )).T, fmt='%20.10e')
        np.savetxt('../inputs/iceplume_scalar.txt',
                np.array([roms_input['N'], roms_input['dx'][ri],
                          roms_input['dy'][ri], roms_input['dt'],
                          roms_input['sgtrs'][ti, ri], roms_input['sgtemp'][ti, ri],
                          roms_input['sgsalt'][ti, ri], roms_input['sgdye'][ti, :, ri],
                          roms_input['sgtyp'][ri], roms_input['sgdep'][ri],
                          roms_input['sglen'][ri]
                         ]), fmt='%20.10e')

        # ------------ run the executable --------------------------------------
        subprocess.call('cd ..; ./iceplume_test.exe', shell=True)

        # ------------ load results from txt files -----------------------------
        data_zr = np.loadtxt('../outputs/iceplume_zr.txt', skiprows=1)
        data_zw = np.loadtxt('../outputs/iceplume_zw.txt', skiprows=1)
        data_dye = np.loadtxt('../outputs/iceplume_dye.txt')

        iceplume_out['zw'][i, :, j] = data_zw[:, 1]
        iceplume_out['f'][i, :, j] = data_zw[:, 2]
        iceplume_out['w'][i, :, j] = data_zw[:, 3]
        iceplume_out['t'][i, :, j] = data_zw[:, 4]
        iceplume_out['s'][i, :, j] = data_zw[:, 5]
        iceplume_out['a'][i, :, j] = data_zw[:, 6]
        iceplume_out['zr'][i, :, j] = data_zr[:, 1]
        iceplume_out['ent'][i, :, j] = data_zr[:, 2]
        iceplume_out['det'][i, :, j] = data_zr[:, 3]
        iceplume_out['m'][i, :, j] = data_zr[:, 9]
        iceplume_out['dye'][i, :, j] = data_dye[2:]


# In[10]:


# --------------------- save outputs ----------------------
pickle.dump(iceplume_out, open('./py_iceplume.pickle', 'wb'))

