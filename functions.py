import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import glob
from PIL import Image
import os
import matplotlib
import matplotlib.colors
import matplotlib.cm as cm
from input_parameters import Sw_IRR, HYD_COND, RHO_I, RHO_W
from input import dt, dx, dX, z


# ===============================================================================
### GRIDDING
# ===============================================================================

def linear_grid(length, nr_cells):
    """Returns
    - grid cell 'width' (dx)
    - number of grid cells per km (nx_km)
    - an array of x-coordinates (x)
    - an array of length (nr_cells - 1) of dx-values
    :param length: total transect length [m]
    :param nr_cells: total number of grid cells in along-transect (X-)direction
    :return:
    """
    dx = length / nr_cells
    nx_km = int(1000 / dx)
    x = np.arange(0, length, dx)
    dX = np.diff(x)
  
    return [dx, nx_km, x, dX]


def gridheight_linear(length, nr_cells, slope_avg):
    """Function calculating/assigning grid cell heights for a linear slope along the full transect.
    Returns
    - an array of grid cell heights (z)
    - an array of length (nr_cells - 1) of dz-values
    :param length: total transect length [m]
    :param nr_cells: total number of grid cells in along-transect (X-)direction
    :param slope_avg: average slope along transect (height difference [m] / transect length [m]
    :return:
    """
    z = []
    dz = slope_avg * length
    dz_x = dz / nr_cells
    for i in np.arange(nr_cells):
        z.append(dz - dz_x * i - 0.5 * dz_x)
    dZ = np.diff(z)
    dZ = np.append(dZ, dZ[-1])
    dz = np.average(dZ)

    return [dz, z, dZ]


# ===============================================================================
### CALCULATIONS
# ===============================================================================

def hydhead(z, psi):
    """Function calculating the hydraulic head in each grid cell
    :param z: height of grid cell [m]
    :param psi: water table height in grid cell [m]
    :return: hydraulic head [m]
    """
    h = z + psi

    return h


def hydgrad(h, dx):
    """Calculate the 'downward' hydraulic gradient in each grid cell (hydraulic gradient between grid cell i and i+1)
    :param h: hydraulic head [m]
    :param dx: grid cell width [m]
    :return: hydraulic gradient [m/m]
    """
    dh = np.diff(h)

    dhdx = dh / dx
    dhdx = np.append(dhdx, dhdx[-1])

    return dhdx


def matrix_mwe_to_height(porosity, height_mwe):
    """
    Converts matrix height [m w.e.] to height in m
    :param porosity: matrix porosity [-]
    :param height_mwe: height to be converted into m [m w.e.]
    :return:
    """
    height = height_mwe / ((RHO_I / RHO_W) * (1 - porosity))

    return height


def matrix_height_to_mwe(porosity, height):
    """
    Converts matrix height [m] to  height in m w.e.
    :param porosity: matrix porosity [-]
    :param height: height to be converted into m w.e. [m]
    :return:
    """
    height_mwe = height * ((RHO_I / RHO_W) * (1 - porosity))

    return height_mwe


def water_mwe_to_height(porosity, height_mwe):
    """
    Converts water height [m w.e.] to height in m.
    :param porosity: matrix porosity [-]
    :param height_mwe: height to be converted into m [m w.e.]
    :return:
    """
    height = height_mwe / porosity

    return height


def water_height_to_mwe(porosity, height):
    """
    Converts water height [m] to height in m w.e.
    :param porosity: matrix porosity
    :param height: height to be converted into m w.e. [m]
    :return:
    """
    height_mwe = height * porosity

    return height_mwe


def snowpack_flux(melt, storage, Sw_irr, height):
    """Calculate retention and outflux of meltwater in snowpack. The following condition is assumed to be met:
            travel time (percolation velocity * snowpack height) << timestep
    Using percolation velocities of 8-9 m hr^-1, snowpack heights in the meter-scale and timesteps of 1 hr it is assumed
    that all supplied excess meltwater (i.e. > storage capacity) is directly transferred into the lateral transport scheme.
    No surface lowering is included in this function.
    :param melt: supplied meltwater [m w.e.]
    :param storage: water equivalent stored in snowpack [m w.e.]
    :param Sw_irr: irreducible water content / snow retention capacity of complete snowpack (porosity-independent!) [-]
    :param height: total height of snowpack [m w.e.] in previous timestep, including snowpack storage
    :return:
    """
    output = pd.DataFrame(columns=['storage', 'outflow', 'height'], index=storage.index)

    for i, j in enumerate(melt):
        k = melt.index[i]
        if j + storage[k] > Sw_irr * height[k]:
            output.loc[k, 'storage'] = Sw_irr * height[k]
            output.loc[k, 'outflow'] = j + storage[k] - Sw_irr * height[k]
        else:
            output.loc[k, 'storage'] = j + storage[k]
            output.loc[k, 'outflow'] = 0

    return [output.storage, output.outflow]


def snowpack_flux_melt(melt, storage, Sw_irr, height):
    """
    As above, but including surface lowering due to melt.
    :param melt: supplied meltwater, can be either positive (melt) or zero (refreezing) [m w.e.]
    :param storage: water equivalent stored in snowpack [m w.e.]
    :param Sw_irr: irreducible water content / snow retention capacity of complete snowpack (porosity-independent!) [-]
    :param height: total height of snowpack at the start of the timestep [m w.e.]
    :return:
    """
    output = pd.DataFrame(columns=['storage', 'outflow', 'height_new'])

    for k in melt.index:
        j = melt[k]

        if j > 0:
            height_new = max(height[k] - j, 0)
            storage_max = height_new * Sw_irr
            if storage[k] > storage_max:
                output.loc[k, 'storage'] = storage_max
                output.loc[k, 'outflow'] = j + (storage[k] - storage_max)
                output.loc[k, 'height_new'] = height_new
            else:
                if j + storage[k] > storage_max:
                    output.loc[k, 'storage'] = storage_max
                    output.loc[k, 'outflow'] = j - (storage_max - storage[k])
                    output.loc[k, 'height_new'] = height_new
                else:
                    output.loc[k, 'storage'] = j + storage[k]
                    output.loc[k, 'outflow'] = 0
                    output.loc[k, 'height_new'] = height_new

        else:  # j = 0, i.e. refreezing/no melt
            height_new = height[k]
            output.loc[k, 'storage'] = storage[k]
            output.loc[k, 'outflow'] = 0
            output.loc[k, 'height_new'] = height_new

    return [output.storage, output.outflow, output.height_new]


def snowpack_flux_melt_runoff(melt, storage, Sw_irr, height):
    """
    As above, but including surface lowering due to melt and surface runoff. Introduces the term "excess melt" which
    captures the melt that cannot be accommodated by surface lowering/snowpack melt.
    :param melt: supplied meltwater, can be either positive (melt) or zero (refreezing) [m w.e.]
    :param storage: water equivalent stored in snowpack [m w.e.]
    :param Sw_irr: irreducible water content / snow retention capacity of complete snowpack (porosity-independent!) [-]
    :param height: total height of snowpack at the start of the timestep [m w.e.]
    :return:
    """
    output = pd.DataFrame(columns=['storage', 'outflow', 'height_new', 'runoff', 'excess_melt'])

    for k in melt.index:
        j = melt[k]

        if j > 0:
            height_new = max(height[k] - j, 0)
            storage_max = height_new * Sw_irr
            if height_new > 0:                          # no runoff if melt < snowpack height
                excess_melt = 0
                if storage[k] > storage_max:            # instantaneous percolation of all meltwater since snowpack is fully saturated
                    storage_updated = storage_max
                    outflow = j + (storage[k] - storage_max)
                else:
                    if j + storage[k] > storage_max:    # saturate snowpack before transferring water into lateral flow domain
                        storage_updated = storage_max
                        outflow = j - (storage_max - storage[k])
                    else:
                        storage_updated = j + storage[k]
                        outflow = 0
            else:                                       # all snowpack is melted away, remaining melt runs off at the surface
                storage_updated = 0
                outflow = storage[k]
                excess_melt = j - storage[k]
        else:                                           # no changes in snowpack height and -storage
            height_new = height[k]
            storage_updated = storage[k]
            outflow = 0
            excess_melt = 0

        output.loc[k, 'storage'] = storage_updated
        output.loc[k, 'outflow'] = outflow
        output.loc[k, 'height_new'] = height_new
        output.loc[k, 'excess_melt'] = excess_melt

    return [output.storage, output.outflow, output.height_new, output.excess_melt]


def snowpack_flux_SL(melt, storage, Sw_irr, height):
    """Calculate retention and outflux of meltwater in snowpack. The following condition is assumed to be met:
            travel time (percolation velocity * snowpack height) << timestep
    Using percolation velocities of 8-9 m hr^-1, snowpack heights in the meter-scale and timesteps of 1 hr it is assumed
    that all supplied excess meltwater (i.e. > storage capacity) is directly transferred into the lateral transport scheme.
    No surface lowering is included in this function.
    :param melt: supplied meltwater [m w.e.]
    :param storage: water equivalent stored in snowpack [m w.e.]
    :param Sw_irr: irreducible water content / snow retention capacity of complete snowpack (porosity-independent!) [-]
    :param height: total matrix height [m w.e.] in previous timestep, including snowpack storage
    :return:
    """
    output = pd.DataFrame(columns=['storage', 'outflow', 'height'], index=storage.index)

    for k in melt.index:
        j = melt[k]

        if j > 0:
            storage_max = max(height[k] - j, 0) * Sw_irr
            if storage[k] > storage_max:
                output.loc[k, 'storage'] = storage_max
                output.loc[k, 'outflow'] = j + (storage[k] - storage_max)
            else:
                if j + storage[k] > storage_max:
                    output.loc[k, 'storage'] = storage_max
                    output.loc[k, 'outflow'] = j - (storage_max - storage[k])
                else:
                    output.loc[k, 'storage'] = j + storage[k]
                    output.loc[k, 'outflow'] = 0

        else:  # j = 0, i.e. refreezing/no melt
            output.loc[k, 'storage'] = storage[k]
            output.loc[k, 'outflow'] = 0
          
    return [output.storage, output.outflow]


def lateral_flux(dt, dx, dX, z, K_lateral, inflow, storage):
    """ Function calculating the lateral meltwater outflow and the resulting new water table height (meltwater storage)
    in case surface runoff is not applied.
    :param dt: timestep length [s]
    :param dx: typical length of a grid cell (distance to travel laterally) [m]
    :param dX: array of delta X [m]
    :param z: elevation head of grid cell [m]
    :param K_lateral: hydraulic conductivity [m s^-1]
    :param inflow: inflow from snowpack (vertical percolation) [m]
    :param storage: water table height in previous timestep (water 'stored' in grid cell) [m]
    :return:
    outflow: total volume moving out of fast reservoir [m]
    watertable: new water table height for timestep t [m]
    dhdx_final: new hydraulic gradient at the end of timestep t [m m^-1]
    """
    A = storage + inflow                        # calculate new hydraulic head after inflow

    dhdx_new = -(A + z).diff(periods=-1) / np.append(dX, dX[-1])
    dhdx_new.iloc[-1] = hydgrad(z, dx)[-1]

    dhdx_new.loc[dhdx_new > 0] = 0              # ignore any 'backward'/upslope flow

    outflow = -K_lateral * dhdx_new * A * dt
    outflow /= dx

    outflow.loc[abs(outflow) > A] = A           # ensure only water present in a grid cell flows out

    watertable = A - outflow
    dhdx_final = (watertable + z).diff(periods=-1) / np.append(dX, dX[-1])

    return [outflow, watertable, dhdx_final]


def darcy_flow_SL(porosity, melt, matrixheight, snowstorage, watertable):
    """Function calculating storage in irreducible water saturation in the snow/matrix, updated matrix height and water
    table height as a function of melt input and values of snow storage, water table and matrix height in the previous
    timestep. Assumes any excess water (> snow height) runs of laterally and instantaneously. All calculations are done
    in m w.e. so necessary conversions are made upfront and afterward.
    :param porosity:
    :param melt: supplied meltwater [m w.e.]
    :param matrixheight: total height of matrix material in previous timestep [m w.e.]
    :param snowstorage: water volume stored in snowpack in previous timestep [m w.e.]
    :param watertable: height of water table in previous timestep [m w.e.]
    :return:
    """
    runoff = pd.Series(data=0, index=matrixheight.index)
    runoff_location = pd.Series(index=matrixheight.index)
    iceslabmelt = pd.Series(data=0, index=matrixheight.index)
    matrixheight_new = pd.Series(data=None, index=matrixheight.index)
    watertable_calc = pd.Series(data=None, index=matrixheight.index)
    inflow = pd.Series(data=None, index=matrixheight.index)

    # calculate new amount of storage in snowpack and amount of melt percolation into the water table
    [snowstorage_new, percolation] = snowpack_flux_SL(melt, snowstorage, Sw_IRR, matrixheight - watertable)

    # check if there are grid cells in which the complete slush matrix melts away (melt > slush matrix), all the water
    # previously in this matrix is added to surface runoff
    empty_mask = matrixheight - melt < 0
    if any(empty_mask):     # add melted matrix to runoff where no more matrix is left, calculate updated matrix height
        empty_index = empty_mask.index[empty_mask]
        matrixheight_new[empty_index] = 0
        matrixheight_new[~empty_mask] = matrixheight[~empty_mask] - melt[~empty_mask]
        iceslabmelt.loc[empty_index] = -(matrixheight - melt).loc[empty_index]
        runoff.loc[empty_index] += watertable.loc[empty_index] + percolation.loc[empty_index]
        watertable_calc.loc[empty_index] = 0
        watertable_calc.loc[~empty_mask] = (watertable + percolation).loc[~empty_mask]
    else:   # slush matrix > melt, i.e. all meltwater can be accommodated by melt of existing slush matrix
        # calculate updated matrix height
        matrixheight_new = matrixheight - melt
        watertable_calc = watertable + percolation

    # Calculate lateral flow 
    dhdx_new = -(watertable_calc + z).diff(periods=-1) / np.append(dX, dX[-1])
    dhdx_new.iloc[-1] = hydgrad(z, dx)[-1]
    dhdx_new.loc[dhdx_new > 0] = 0                      # ignore any 'backward'/upslope flow

    outflow = - HYD_COND * watertable_calc * dhdx_new * dt
    outflow /= dx
    outflow.loc[abs(outflow) > watertable_calc] = watertable_calc     # ensure only water present in a grid cell flows out

    inflow.iloc[0] = 0
    inflow.iloc[1:] = outflow[:-1]

    watertable_new = watertable_calc - outflow + inflow

    # make sure that also after lateral inflow there is no excess water (i.e. water table height <= matrix height)
    height_water_new = water_mwe_to_height(porosity, watertable_new)
    height_matrix = matrix_mwe_to_height(porosity, matrixheight_new)

    excess_mask = height_water_new > height_matrix
    if any(excess_mask):
        excess_index = excess_mask.index[excess_mask]
        runoff[excess_index] += water_height_to_mwe(porosity, height_water_new - height_matrix)[excess_index]
    watertable_new -= runoff
    watertable_new = watertable_new.clip(lower=0)
    if any(runoff):
        runoff_location = z[min((runoff[runoff>0]).index.values)]

    dhdx_final = -(watertable_new + z).diff(periods=-1) / np.append(dX, dX[-1])

    return [matrixheight_new, snowstorage_new, watertable_new, outflow, runoff, runoff_location, iceslabmelt, dhdx_final]


def detect_error_condition(start_error, label, ts, flag):
    """
    check whether error condition (= flow not everywhere to the right) occurs
    if yes, safe the index of the time step. This will be done only once, the first time
    error conditions are observed. The information is then used to plot water level for the
    concerned time step.
    """
    if start_error[0] == 0:
        if not all([x == 'down' for x in flag.flag]):
            start_error[0] = label
            start_error[1] = ts
            print('* Error condition started at time step: ' + label + '; index ' + str(int(ts)) + ' *')

    return start_error


# ===============================================================================
### PLOTTING
# ===============================================================================

def plot_overview(days, t_index, tsppday, tspsmelt, length, dx, var1, var2, var3, var1_label, var2_label, var3_label, xlimits):
    fig, ax = plt.subplots(nrows=3, figsize=(28, 20), gridspec_kw={'height_ratios': [3, 1, 1]})
    colors = plt.cm.brg(np.linspace(start=0, stop=1, num=days))

    for j, label in enumerate(t_index):
        if j % tsppday == 0:
            if j <= tspsmelt:
                ax[0].plot(np.arange(0, length) * dx, var1.loc[label], ',', linestyle='dashed', label='%.f days' % (j / tsppday), colors=colors[int(j / tsppday)])
                ax[1].plot(np.arange(0, length) * dx, var2.loc[label], ',', linestyle='dashed', label='%.f days' % (j / tsppday), colors=colors[int(j / tsppday)])
            if j > tspsmelt:
                ax[0].plot(np.arange(0, length) * dx, var1.loc[label], ',', linestyle='dotted', label='%.f days' % (j / tsppday), colors=colors[int(j / tsppday)])
                ax[1].plot(np.arange(0, length) * dx, var2.loc[label], ',', linestyle='dotted', label='%.f days' % (j / tsppday), colors=colors[int(j / tsppday)])

    ax[2].plot(np.arange(0, length) * dx, var3, label=var3_label,)

    ax[0].set_ylabel(var1_label)
    ax[0].set_ylim()
    ax[0].set_xlim(xlimits)

    ax[1].set_ylabel(var2_label)
    ax[1].set_ylim()
    ax[1].set_xlim(xlimits)

    ax[2].set_ylabel(var3_label)
    ax[2].set_xlim(xlimits)


def make_gif(frame_folder, output_folder, output_fn):
    """
    Function creating a GIF for all .png-files in a certain (output-)folder
    :param frame_folder: location of .png-files to be used as individual frames
    :param output_folder: location for output file
    :param output_fn: name of output file
    :return:
    """
    timing = 500  # duration of 1 frame [ms]

    frames = [Image.open(image) for image in glob.glob(f'{frame_folder}/*.png')]
    frame_one = frames[0]
    frame_one.save(os.path.join(output_folder, output_fn), format='GIF', append_images=frames, save_all=True,
                   duration=timing, loop=1)


def create_bins(T_RES, t_index, t, dx):
    """
    Function returning bins
    :param T_RES: variable
    :param t_index:
    :param t:
    :param dx:
    :return:
    """
    V_max = dx / min([min(T_RES.loc[t_index[ts + 1]]) for ts, time in enumerate(t[1:])]) * 3600
    V_min = dx / max([max(T_RES.loc[t_index[ts + 1]]) for ts, time in enumerate(t[1:])]) * 3600
    # ignores first (initialisation) timestep

    dT = (V_max - V_min) / 25
    bins = np.arange(V_min, V_max, dT)

    return [bins, V_min, V_max]


def define_color(val_min, val_max, binval):
    """
    Function returning color based on bin value of item.
    :param val_min: minimum parameter value
    :param val_max: maximum parameter value
    :param binval: bin number of item to be plotted/colored
    :return:
    color: RGB-code
    """
    norm = matplotlib.colors.Normalize(vmin=val_min, vmax=val_max)
    cmap = cm.bone
    s_m = cm.ScalarMappable(cmap=cmap, norm=norm)
    s_m.set_array([])
    color = cm.bone(norm(binval))

    return [color, s_m]
  
