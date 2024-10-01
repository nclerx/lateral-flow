import datetime
import os
import numpy as np
import pandas as pd
from scipy import interpolate


def import_grid(input_fn):
    """"Function reading x- and z-coordinates from a specified xlsx-file
    :param input_fn: complete path of input file
    :return:
    """
    inputdata = pd.read_excel(input_fn, header=0, index_col=0, engine='openpyxl')
    inputdata.columns = inputdata.columns.str.lower()
    inputdata = inputdata.set_index(np.round(inputdata.index))  # rounding off the x-coordinates to the m-scale
    x_array = np.array(list(inputdata.index))
    z_array = np.array(inputdata.iloc[:, 2])
    x_start = x_array[0]
    x_end = x_array[-1]
    dx_total = x_end - x_start
    if not all(element != np.diff(x_array)[0] for element in np.diff(x_array)):  # check if input data are evenly spaced
        delta_x = 100
        n_km = 10
        x_array2 = np.arange(0, dx_total - delta_x, delta_x)
        z_array2 = np.interp(x_array2, inputdata.index, inputdata.iloc[:, 2])
        # return [delta_x, len(x_array2), n_km, x_array2, z_array2]
        x_array3 = x_array2
        x_nrs = np.arange(0, len(x_array2))
        z_array3 = z_array2
    else:
        delta_x = np.diff(x_array)[0]
        n_km = 1000 / delta_x
        # return [delta_x, len(x_array), n_km, x_array, z_array]
        x_array3 = x_array
        x_nrs = np.arange(0, len(x_array3))
        z_array3 = z_array
    return [delta_x, len(x_array3), n_km, x_nrs, x_array3, z_array3]


# Set path to input directory
# input_dir_base = r'C:/Users/ClerxN/switchdrive/PhD/modelling/1DFM_input/'
input_dir_base = r'C:/Users/ClerxN/switchdrive/PhD/_pchapter3/input'


# ===============================================================================
### OUTPUT PARAMETERS
# ===============================================================================

# output_dir_base = r'C:/Users/ClerxN/switchdrive/PhD/modelling/1DFM_output'
output_dir_base = r'C:/Users/ClerxN/switchdrive/PhD/_paper2/sensitivity'
outputname = 'noSL_noSI_2020_BBH'


# ===============================================================================
### PHYSICAL PROCESSES
# ===============================================================================

surfacelowering = False
SIformation = False
surfacerunoff = False
gridsmoothing = True         # surface slope smoothing or not (5km-rolling average)
snowpackwarming = False      # no melt percolation before first runoff according to SEBM (bucket) model


# ===============================================================================
### MELT INPUT PARAMETERS
# ===============================================================================

# Choose one of the following options for melt input:
# melt_inputfile = None
# melt_inputfile = 'manual'
melt_inputfile = 'SEBM'

melt_inputfn_manual = os.path.join(input_dir_base, 'melt01.xlsx')
melt_inputfn_SEBM = os.path.join(input_dir_base, 'meltwater(mm)_IP_2023_bins.csv')

if melt_inputfile:
    if melt_inputfile == 'manual':
        melt_inputfn = melt_inputfn_manual
    elif melt_inputfile == 'SEBM':
        melt_inputfn = melt_inputfn_SEBM

if melt_inputfile == 'manual':
    meltdata = pd.read_excel(melt_inputfn, engine='openpyxl').dropna(axis=0, how='all')
elif melt_inputfile == 'SEBM':
    meltdata = pd.read_csv(melt_inputfn, engine='python').dropna(axis=0, how='all')
else:
    melt_input = 1e-3 / 3600  # melt rate [m/s]
    meltdays = 7  # melt duration [days]


# ===============================================================================
### SNOWPACK WARMING
# ===============================================================================

if snowpackwarming:
    runoff_inputfn = os.path.join(input_dir_base, 'runoff(mm)_IP_2023_bins.csv')
    runoffdata = pd.read_csv(runoff_inputfn, engine='python').dropna(axis=0, how='all')


# ===============================================================================
### GRID INPUT PARAMETERS
# ===============================================================================

lineargrid = True
manualgrid = False
gentlesteepgentlegrid = False
gridinputfile = False

gridcell_height = 2.0  # unit height of 1 grid cell [m w.e.]


# ===============================================================================
# ## INPUT FILE
# ===============================================================================

# Set path to input file describing slope of model transect
# grid_inputfn = r'C:/Users/ClerxN/switchdrive/PhD/modelling/1DFM_input/K_transect_main.xlsx'
# grid_inputfn = r'C:/Users/ClerxN/switchdrive/PhD/modelling/1DFM_input/K_transect_short.xlsx'
grid_inputfn = os.path.join(input_dir_base, 'Ktransectflowline.xlsx')
# grid_inputfn = os.path.join(input_dir_base, 'hydrological_channel_01.xlsx')     # only between 1900 and 1600 m a.s.l.!


# ===============================================================================
### MANUAL INPUT
# ===============================================================================
z0 = float(input('highest elevation? [m] '))
z_base = float(input('lowest elevation? [m] '))
if lineargrid or manualgrid or gentlesteepgentlegrid:
    dx = 100    # [m] width of each grid cell
    nx_km = 10
    dx_input = import_grid(grid_inputfn)[0]
    z_Ktransect = import_grid(grid_inputfn)[-1]
    # build in check
    z0_pos = np.abs(np.array(z_Ktransect)-z0).argmin()
    z_base_pos = np.abs(np.array(z_Ktransect)-z_base).argmin()
    if dx_input == dx:
        L = z_base_pos - z0_pos
    else:
        L = ((z_base_pos - z0_pos) * dx_input) / dx
    slope = -(z0 - z_base) / (L * dx)
    x = np.arange(0, L, 1)  # IDs of the grid cells e.g. [0, 1, 2, 3, 4, 5, 6, 7, 8]

if lineargrid:
    z = (z0 + x * slope * dx).tolist()  # [m a.s.l.] surface topography, e.g. [1, 0.985, 0.970, 0.955]
elif manualgrid:
    z = [10, 9.95, 9.9, 9.85, 9.75, 9.65, 9.55, 9.55, 9.55, 9.45, 9.4, 9.35, 9.3, 9.25, 9.2]
elif gentlesteepgentlegrid:
    # specify the irregularities in the slope
    # One can choose as many slopes as needed and over arbitrary distances
    # The slopes that differ from standard "slope" should not overlap. It is not possible to use
    # a slope of absolutely zero as np.arange cannot deal with zero steps
    x_ds_s = [60]   # [grid index] starting point of area of different slope (ds)
    L_ds = [40]     # [grid cells] length of area of different slope
    slope_ds = [-(20 / 1000)]  # [m/m] slope (cannot be zero)

    # now establish the array of the grid cell elevations
    for i in range(0, len(x_ds_s)):
        z = np.arange(z0, 0, slope * dx)[0:L]
        z[x_ds_s[i]:x_ds_s[i] + L_ds[i]] = np.arange(z[x_ds_s[i]], 0, slope_ds[i] * dx)[0:L_ds[i]]
        z[x_ds_s[i] + L_ds[i] - 1:] = np.arange(z[x_ds_s[i] + L_ds[i] - 1], 0, slope * dx)[
                                      0:L - x_ds_s[i] - L_ds[i] + 1]
    z = z.tolist()

else:
    [dx, L_full, nx_km, x_full, x_coords_full, z_full] = import_grid(grid_inputfn)
    z0_pos = np.abs(np.array(z_full)-z0).argmin()
    if gridsmoothing:
        z_pd = pd.DataFrame(z_full)
        z_pd = z_pd.rolling(window=50, center=True, min_periods=5).mean()
        z_full_smooth = z_pd[0].values
        z_base_pos = np.abs(np.array(z_full_smooth)-z_base).argmin()
    # complete to use only subsection of input file!
        z = z_full_smooth[z0_pos:z_base_pos]
    else:
        z_base_pos = np.abs(np.array(z_full)-z_base).argmin()
        z = z_full[z0_pos:z_base_pos]
    L = len(z)
    x = np.arange(0, L, 1)

dX = np.diff(x) * dx
dZ = np.diff(z)


# ===============================================================================
### TIMESTEP CALCULATIONS
# ===============================================================================

if melt_inputfile == 'manual':
    dates = [
        datetime.datetime(int(meltdata.year[j]), int(meltdata.month[j]), int(meltdata.day[j]), meltdata.time[j].hour)
        for j in np.arange(len(meltdata))]
    t_0 = dates[0]
    t_end = dates[-1]
    delta_t = dates[1] - dates[0]
    dt = 3600  # in seconds

    nt_meltinput = len(meltdata)
    if delta_t.seconds == dt:
        nt = nt_meltinput
    else:
        nt = (delta_t.seconds / dt) * nt_meltinput
    nt_day = (24 * 3600) / dt  # number of time steps per day [-]
    ndays = (t_end - t_0).days

    t = np.arange(0, nt * dt, dt)

    date_start = datetime.datetime.strftime(t_0, '%d/%m/%Y %H:%M')
    date_end = datetime.datetime.strftime((t_end + delta_t), '%d/%m/%Y %H:%M')

    t_index = [(t_0 + datetime.timedelta(seconds=int(dt * time))).strftime('%Y-%m-%d %H:%M:%S') for time in
               np.arange(nt)]

    meltrate = pd.DataFrame(index=t_index, columns=x)
    for i, j in enumerate(t_index):
        if delta_t.seconds != dt:
            d = int(i / (delta_t.seconds / dt))
            meltrate.loc[j, :] = meltdata['melt (mm)'][d] * 1e-3 / (dt * (delta_t.seconds / dt))
        else:
            meltrate.loc[j, :] = meltdata['melt (mm)'][i] * 1e-3 / dt

elif melt_inputfile == 'SEBM':
    meltdata['newDay'] = meltdata['Day'].astype(str).str.zfill(3)
    meltdata['newHour'] = meltdata['Hour'].astype(str).str.zfill(2)
    dates = list(pd.to_datetime(meltdata.Year.astype(str) + meltdata.newDay.astype(str) + meltdata.newHour.astype(str),
                                format='%Y%j%H'))
    meltdata.index = dates
    t_0 = dates[0]
    t_1 = dates[-1]

    delta_t = dates[1] - dates[0]
    dt = 3600

    meltrate_elev = (meltdata.iloc[:, 3:-2].astype(float) / 3600)
    colnames = []
    for i in meltrate_elev.columns:
        colnames.append(int(i[:-1]))
    meltrate_elev.columns = colnames

    t_start = input('start date (dd/mm/YYYY)? ') + ' 00:00:00'
    t_end = input('end date (dd/mm/YYYY)? ') + ' 00:00:00'
    date_start = datetime.datetime.strptime(t_start, '%d/%m/%Y %H:%M:%S')
    date_end = datetime.datetime.strptime(t_end, '%d/%m/%Y %H:%M:%S')

    meltrate = pd.DataFrame(columns=z[::-1], index=meltrate_elev.loc[date_start:date_end].index, dtype=float)

    for i in meltrate_elev.loc[date_start:date_end].index:
        meltrate.loc[i, :] = np.interp(z, meltrate_elev.columns.values, meltrate_elev.loc[i, :] * 1e-3)

    nt_day = (24 * 3600) / dt
    nt = ((date_end - date_start) * nt_day).days
    nt_meltinput = nt
    ndays = (date_end - date_start).days

    t = np.arange(meltrate.index.get_loc(date_start), nt * dt, dt)
    t_index = [(date_start + datetime.timedelta(seconds=int(dt * time))).strftime('%Y-%m-%d %H:%M:%S') for time in
               np.arange(nt)]

    if snowpackwarming:
        runoffdata = runoffdata.iloc[:, 3:]
        runoffdata.columns = colnames
        runoffdata.index = dates

        runoff_start = pd.DataFrame(index=colnames, columns=['start'], data=runoffdata.loc[date_start:date_end, :].ne(0).idxmax(0))

        runoff_start['seconds'] = [(i - date_start).total_seconds() for i in runoff_start.start]
        runoff_start['seconds'] = [None if j == 0 else j for j in runoff_start['seconds']]  # workaround to avoid runoff starting on April 1st in case there is no runoff at all

        f = interpolate.interp1d(runoff_start.seconds.dropna().index, runoff_start.seconds.dropna(), kind='linear', copy=True, fill_value='extrapolate')
        runoff_start['seconds'] = f(runoff_start.index)

        full_delay_seconds = np.interp(z[::-1], colnames, runoff_start.seconds)

        # create list of runoff/melt starting datetimes for all z-values
        melt_start = pd.DataFrame(index=z, columns=['start'], data=[date_start + datetime.timedelta(seconds=i) for i in full_delay_seconds])

        for i in z:
            meltrate.loc[date_start:melt_start.loc[i, 'start'], i] = 0

    meltrate = meltrate.loc[t_index].set_axis(x, axis=1)

else:
    ndays = 21  # number of days to simulate [-]

    date_start = '01/06/2022  00:00'
    t_0 = datetime.datetime.strptime(date_start, '%d/%m/%Y %H:%M')
    t_end = t_0 + pd.Timedelta(days=ndays)
    date_end = datetime.datetime.strftime(t_end, '%d/%m/%Y %H:%M')

    dt = 1 * 3600  # timestep duration [s] --> 6 hour
    nt_day = (24 * 3600) / dt  # number of time steps per day [-]
    nt = ndays * nt_day  # number of timesteps [-]
    nt_meltinput = meltdays * nt_day  # number of time steps during which melt input occurs [-]

    t = np.arange(0, ndays * nt_day)
    t_index = [(t_0 + datetime.timedelta(seconds=int(dt * time))).strftime('%Y-%m-%d %H:%M:%S') for time in
               np.arange(nt)]

    meltrate = pd.DataFrame(columns=x, index=t_index)
    for i in np.arange(int(nt)):
        if i <= nt_meltinput:
            meltrate.loc[t_index[i], :] = melt_input
        else:
            meltrate.loc[t_index[i], :] = 0

t_index_refreeze = ['first_refreeze', *t_index]
melt_ts = meltrate * dt
