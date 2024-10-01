from constants import *
import numpy as np
import os
import glob

from input import surfacelowering, SIformation, surfacerunoff, outputname, output_dir_base


# ===============================================================================
### PARAMETER FUNCTIONS
# ===============================================================================

def porosity(rho):
    """
    Calculates initial porosity (assuming a dry snowpack).
    :param rho: bulk snow density [-]
    :return: porosity [-]
    """
    por = 1 - rho / RHO_W

    return por


def carmankozeny(por, grain_diameter):
    """
    Calculates permeability based on Carman-Kozeny approximation, which relates a medium’s permeability to pressure
    drop and fluid viscosity for laminar flow through a packed bed of solids (Kozeny, 1927; Carman, 1937; Bear, 1972)
    :param por: [-]
    :param grain_diameter: average grain size [m]
    :return: permeability [m^2]
    """
    permeability = EPSILON ** 2 * (por ** 3 * grain_diameter ** 2) / (150 * (1 - por) ** 2)

    return permeability


def calonne(density, grain_diameter):
    """
    Parametrisation by Calonne et al. (2012) to calculate permeability based on specific surface area, density and
    microstructural anisotropy. The equivalent sphere radius relates the specific surface area of a snow particle to
    ice density (German, 1996).
    r_es: equivalent sphere radius [= 3 / (SSA * rho_ice)]      [m]
    rho_s: snow density [kg m^-3]
    k = (3.0 ± 0.3) r_es**2 exp(−0.0130 ± 0.0003)*rho_s
    For simplicity, here the variability is not taken into account.
    :param density [kg m^-3]
    :param grain_diameter: average grain size [m]
    :return: permeability [m^2]
    """
    permeability = 3.0 * (grain_diameter / 2) ** 2 * np.exp(0.0130 * density)

    return permeability


def hydcond(vel_lateral, grad_meas):
    """Function calculating the hydraulic conductivity based on the measured lateral meltwater flow velocity.
    :param vel_lateral: measured lateral flow velocity [m s^-1]
    :param grad_meas: hydraulic gradient (slope) of transect along which lateral flow velocity was measured [m/m]
    :return:
    """
    conductivity = vel_lateral / grad_meas

    return conductivity


# ===============================================================================
### MATRIX/SNOWPACK PARAMETERS
# ===============================================================================

PERM_INIT = carmankozeny(POR_INIT, D)  # matrix permeability following Carman-Kozeny approximation [m^2]
K_CK = (PERM_INIT * G * RHO_W) / MU  # theoretical hydraulic conductivity based on fluid & matrix properties according to Carman-Kozeny [m/s]

HYD_COND = hydcond(U_LATERAL, dhdx_meas)


# ===============================================================================
### OUTPUT FILES
# ===============================================================================

gif = False  # to toggle creation of gif on or off

output_dir = os.path.join(output_dir_base, outputname)
if not os.path.exists(output_dir):
    os.mkdir(output_dir)
if not os.path.exists(output_dir + '\gif') and gif:
    os.mkdir(output_dir + '\gif')
if os.path.exists(output_dir + '\gif'):
    files = glob.glob(output_dir + 'gif')
    for f in files:
        os.remove(f)

if surfacelowering and not SIformation and not surfacerunoff:
    outputnames = ['HYD_HEAD', 'melt_ts', 'SNOWPACK_HEIGHT', 'SNOWPACK_FLUX', 'SNOWPACK_STORAGE', 'LATERAL_FLUX_OUT',
                   'LATERAL_STORAGE', 'TOTAL_HEIGHT', 'ERROR']
    outputnames_old = ['HYD_HEAD', 'melt_ts', 'SNOWPACK_HEIGHT', 'SNOWPACK_FLUX', 'SNOWPACK_STORAGE', 'LATERAL_FLUX_OUT',
                   'LATERAL_STORAGE', 'TOTAL_HEIGHT', 'DISCHARGE']
elif surfacelowering and SIformation and not surfacerunoff:
    outputnames = ['HYD_HEAD', 'melt_ts', 'MATRIX_HEIGHT', 'SNOWPACK_FLUX', 'SNOW_STORAGE', 'WATERTABLE_HEIGHT', 'SI',
                   'DISCHARGE']
    outputnames_old =  ['HYD_HEAD', 'melt_ts', 'SNOWPACK_HEIGHT', 'SNOWPACK_FLUX', 'SNOWPACK_STORAGE', 'LATERAL_FLUX_OUT',
                   'LATERAL_STORAGE', 'TOTAL_HEIGHT', 'SI', 'SI_TOTAL', 'DISCHARGE']
elif surfacelowering and SIformation and surfacerunoff:
    outputnames = ['HYD_HEAD', 'melt_ts', 'MATRIX_HEIGHT', 'SNOW_STORAGE', 'WATERTABLE_HEIGHT', 'SI', 'DISCHARGE',
                   'SURFACE_RUNOFF', 'RUNOFF_LOCATION', 'ICESLABMELT']
else:
    outputnames = ['HYD_HEAD', 'melt_ts', 'SNOWPACK_HEIGHT', 'SNOWPACK_FLUX', 'SNOWPACK_STORAGE', 'LATERAL_FLUX_OUT',
                   'LATERAL_STORAGE', 'DISCHARGE', 'ERROR']
    outputnames_old = ['HYD_HEAD', 'melt_ts', 'SNOWPACK_HEIGHT', 'SNOWPACK_FLUX', 'SNOWPACK_STORAGE', 'LATERAL_FLUX_OUT',
                   'LATERAL_STORAGE', 'DISCHARGE']
outputnames_daily = ['SNOWPACK_FLUX_daily', 'SNOWPACK_STORAGE_daily', 'LATERAL_FLUX_OUT_daily', 'LATERAL_STORAGE_daily',
                     'DISCHARGE_daily']
filename = outputname + '_fulloutput.xlsx'
