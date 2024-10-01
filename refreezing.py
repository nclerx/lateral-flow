from constants import *
from input import *

n = len(DEPTHS_ICESLAB)     # number of layers
store_output = False        # store output y/n

# Check if model timestep is small enough for numerical stability of heat flux- and refreezing, otherwise set timestep to 10 minutes
if dt > 600:
    dt_model = 600
else:
    dt_model = dt
t_final = datetime.datetime.strptime(t_index[-1], '%Y-%m-%d %H:%M:%S')
t_start = datetime.datetime.strptime(t_index[0], '%Y-%m-%d %H:%M:%S')
t_model = np.arange(0, (t_final-t_start).total_seconds() + dt_model, dt_model)


# ===============================================================================
### CALCULATIONS
# ===============================================================================

def refreeze(nz, model_time):
    """"""
    T = np.append(np.insert(np.ones(nz) * T_ICESLAB, 0, 0), T_ICESLAB)  # vector of temperatures for each layer (n+2 to
    # include bottom and top, given that temperatures are calculated at the nodes which are in the middle of the grid
    # cells). Top set to zero, will then be updated each model step.

    # T[1:11] = np.linspace(T_SLUSH, T_ICESLAB, 10) # linear gradient for temperature to evolve from 0 to -10°C in 1 m
    T[1:21] = np.linspace(T_SLUSH, T_ICESLAB, 20)   # linear gradient for temperature to evolve from 0 to -10°C in 2 m

    T_evol = np.ones([nz+2, len(model_time)]) * T_ICESLAB  # array of temperatures for each layer and time step
    for i in range(nz+2):
        T_evol[:, i] = T
    dTdt = np.empty(nz)  # time derivative of temperature at each node
    heatflux = np.zeros([nz+1, len(model_time)])  # array of heat flux for each layer and time step
    refrozen = np.zeros(len(model_time))  # total amount of refreezing to the top per timestep

    for j in range(0, len(model_time)-1):
        T[0] = T_SLUSH
        T[-1] = T_ICESLAB

        dTdt[:] = ALPHA * (-(T[1:-1] - T[:-2]) / DZ_ICESLAB**2 + (T[2:] - T[1:-1]) / DZ_ICESLAB**2)
        T[1:-1] = T[1:-1] + dTdt * dt_model
        T_evol[:, j] = T
        heatflux[:, j] = K_THERM * (T[:-1] - T[1:]) / DZ_ICESLAB
        refrozen[j] = heatflux[0, j] * dt_model / LAT

    return [T_evol, heatflux, refrozen]


[temp, hf, refr] = refreeze(n, t_model)
TEMP = pd.DataFrame(index=DEPTHS_ICESLAB, columns=t_model, data=temp[1:-1])
HEATFLUX = pd.DataFrame(index=DEPTHS_ICESLAB, columns=t_model, data=hf[:-1])
REFREEZE = pd.DataFrame(index=t_model, columns=['mm_cum'], data=np.cumsum(refr))

output_heatflux = [TEMP, HEATFLUX, REFREEZE]
outputnames_heatflux = ['TEMP', 'HEATFLUX', 'REFREEZE']


# ===============================================================================
### STORING OUTPUT
# ===============================================================================

if store_output:
    fn = os.path.join(output_dir_base, 'refreezing_2mgradient.xlsx') # adjust filename if wanted
    with pd.ExcelWriter(fn, engine='openpyxl') as writer:
        [dataframe.to_excel(writer, sheet_name=sheet, startrow=0, startcol=0) for dataframe, sheet in zip(output_heatflux, outputnames_heatflux)]
  
