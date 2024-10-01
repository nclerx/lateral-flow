import numpy as np

#===============================================================================
### CONSTANTS
#===============================================================================

G = 9.81                    # gravitational acceleration [m/s^2]

RHO_W = 1000.               # fluid density of water at 0°C [kg/m^3]
RHO_I = 918.9               # ice density at -10°C [kg/m^3]
RHO_A = 1.1                 # air density at 0°C [kg/m^3]
BETA = 4.4e-10              # fluid compressibility of water [m^2/N]
MU = 0.0017916              # water dynamic viscosity at 0°C [Pa*s]

D = 2e-3                    # average grain diameter for slush matrix particles [m]
EPSILON = 1                 # sphericity (= 1 for perfectly spherical grains) [-]
POR_INIT = 0.45             # initial porosity [-], as measured in summer 2020
Sw_IRR = 0.01               # irreducible water saturation (% of total volume/m w.e., i.e. not corrected for porosity)
# Sw_IRR = 0.07               # irreducible water saturation - high case (% of total volume/m w.e., i.e. not corrected for porosity)

Cp = 2.0e3                  # specific heat capacity of ice @ -10°C [J kg^-1 K^-1]
K_THERM = 2.30              # thermal conductivity of ice @ -10°C [W m^-1 K^-1]
T_ICESLAB = -10             # initial condition of ice slab temperature ('base' temperature of ice slab) [°C]
T_SLUSH = 0                 # boundary condition at ice slab surface - slush/water temperature [°C]
ALPHA = K_THERM / (Cp * RHO_I) # thermal diffusivity [J m^-3 K^-1]
LAT = 334e3                 # latent heat of water [J kg^-1]

Z_ICESLAB = 10              # depth at which ice slab is always at initial temperature [m]
DZ_ICESLAB = 0.1            # ice slab thickness intervals for refreezing calculations [m]
DEPTHS_ICESLAB = np.arange(0, Z_ICESLAB, DZ_ICESLAB)

U_VERTICAL = 0.25 / 3600    # unsaturated flow velocity (vertical flow) - average from ROSA, 0.25 m/hr
U_LATERAL = 192e-5          # lateral flow velocity [m/s] ~ 7.2 m/hr, avg. value of measurements on GrIS in summer 2020
# U_LATERAL = 36.1e-5         # lateral flow velocity [m/s] ~ 1.3 m/hr, low value of measurements on GrIS in summer 2020
# U_LATERAL = 394e-5            # lateral flow velocity [m/s] ~ 14.2 m/hr, high value of measurements on GrIS in summer 2020
dhdx_meas = 5 / 1000        # approx. slope on which summer 2020 lateral flow velocity measurements were performed [m/m]


sort = {'peak': 0, 'down': 1, 'up': 2, 'dip': 3}
