from input_parameters import *
from input import *
from functions import *
import matplotlib.pyplot as plt
import pandas as pd
from refreezing import REFREEZE

plt.close('all')

time_start = datetime.datetime.now()
start_error = [0, 0]  # will record the first time step where error conditions (flow not everywhere to the right) occur
print(time_start, 'Calculations starting')

# ===============================================================================
### INITIALIZING VARIABLES
# ===============================================================================

HYD_HEAD = pd.DataFrame(data=None, index=t_index, columns=x, dtype='float')  # hydraulic head [m]
DHDX = pd.DataFrame(data=None, index=t_index, columns=x, dtype='float')  # hydraulic gradient [m/m]
SNOWPACK_FLUX = pd.DataFrame(data=0, index=t_index, columns=x,
                             dtype='float')  # flow out of snowpack (vertical percolation) [m w.e.]
MATRIX_HEIGHT = pd.DataFrame(data=None, index=t_index, columns=x,
                             dtype='float')  # height of slush matrix (replaces snowpack and lateral storage) [m]
SNOW_STORAGE = pd.DataFrame(data=0, index=t_index, columns=x,
                            dtype='float')  # water stored in snow (irreducible water) [m]
WATERTABLE_HEIGHT = pd.DataFrame(data=0, index=t_index, columns=x, dtype='float')  # water table height [m]

if surfacelowering:
    TOTAL_HEIGHT = pd.DataFrame(data=None, index=t_index, columns=x,
                                dtype='float')  # total height of snowpack (water table + snowpack)
    TOTAL_HEIGHT.iloc[0, :] = gridcell_height
if surfacerunoff:
    excess_melt = pd.DataFrame(data=0, index=t_index, columns=x,
                               dtype='float')  # excess melt in case no remaining snowpack
    SURFACE_RUNOFF = pd.DataFrame(data=0, index=t_index, columns=x,
                                  dtype='float')  # volume of meltwater in excess of the water table per grid cell
    SURFACE_RUNOFF_TOTAL = pd.DataFrame(data=0, index=t_index,
                                        columns=['runoff'])  # total volume of meltwater in excess of the water table
    ICESLABMELT = pd.DataFrame(data=0, index=t_index, columns=x,
                               dtype='float')  # excess melt not accommodated by the vertical or lateral flow domain
    RUNOFF_LOCATION = pd.DataFrame(data=None, index=t_index, columns=['altitude'],
                                   dtype='float')  # highest elevation of surface runoff occurrence
WATERTABLE_FLUX_IN = pd.DataFrame(data=0, index=t_index, columns=x,
                                  dtype='float')  # volume of laterally transported meltwater into grid cell
WATERTABLE_FLUX_OUT = pd.DataFrame(data=0, index=t_index, columns=x,
                                   dtype='float')  # volume of laterally transported meltwater out of grid cell
DISCHARGE = pd.DataFrame(data=None, index=t_index, columns=['outflow'])  # total discharge
SI = pd.DataFrame(data=0, index=t_index, columns=x, dtype='float')  # superimposed ice formed each timestep [m]
SI_TOTAL = pd.DataFrame(data=None, index=t_index, columns=x, dtype='float')  # cumulative superimposed ice [m]
ERROR = pd.DataFrame(data=None, index=t_index, columns=['error'])

MATRIX_HEIGHT.iloc[0, :] = matrix_mwe_to_height(POR_INIT,
                                                gridcell_height)  # calculate initial height of snow/slush matrix

# ===============================================================================
### CALCULATIONS
# ===============================================================================

flag_non_zero = False
label_zero = False

for ts, j in enumerate(t):
    label = t_index[ts]
    label_prev = t_index[ts - 1] if ts > 0 else None

    if all(meltrate.loc[label] == 0) and (label_prev is None or (
            all(SNOW_STORAGE.loc[label_prev] == 0) and all(WATERTABLE_HEIGHT.loc[label_prev] == 0))):
        if not flag_non_zero:
            SNOW_STORAGE.loc[label] = 0
            MATRIX_HEIGHT.loc[label] = gridcell_height
            WATERTABLE_HEIGHT.loc[label] = 0
            HYD_HEAD.loc[label] = np.array(z).T
            DHDX.loc[label] = hydgrad(HYD_HEAD.loc[label], dX)
            DISCHARGE.loc[label] = 0
            ERROR.loc[label] = 0
            continue
    else:
        flag_non_zero = True

    if (ts + 1) % nt_day == 1 and not label_zero:
        print('day', (ts / nt_day) + 1, label[:-9], 'time elapsed', datetime.datetime.now() - time_start)

    if flag_non_zero and all(MATRIX_HEIGHT.loc[label_prev] == 0):
        if not label_zero:
            label_zero = label
            print('no more matrix at', label_zero)
        SNOW_STORAGE.loc[label] = 0
        MATRIX_HEIGHT.loc[label] = 0
        WATERTABLE_HEIGHT.loc[label] = 0
        HYD_HEAD.loc[label] = np.array(z).T
        DHDX.loc[label] = hydgrad(HYD_HEAD.loc[label], dX)
        DISCHARGE.loc[label] = DISCHARGE.loc[label_prev]
        ERROR.loc[label] = ERROR.loc[label_prev]
        if surfacerunoff:
            SURFACE_RUNOFF_TOTAL.loc[label] = SURFACE_RUNOFF_TOTAL.loc[label_prev]
            ICESLABMELT.loc[label] = meltrate.loc[label] * dt
        continue

    # Update HYD_HEAD and DHDX
    HYD_HEAD.loc[label] = np.array(z) + WATERTABLE_HEIGHT.loc[label_prev]
    DHDX.loc[label] = hydgrad(HYD_HEAD.loc[label], dX)

    # Calculate new matrix height, snow storage, water table, lateral outflow and hydraulic gradient
    if surfacelowering and SIformation and surfacerunoff:
        [MATRIX_HEIGHT.loc[label], SNOW_STORAGE.loc[label], WATERTABLE_HEIGHT.loc[label],
         WATERTABLE_FLUX_OUT.loc[label],
         SURFACE_RUNOFF.loc[label], RUNOFF_LOCATION.loc[label], ICESLABMELT.loc[label], DHDX.loc[label]] = (
            darcy_flow_SL(POR_INIT, meltrate.loc[label] * dt, MATRIX_HEIGHT.loc[label_prev],
                          SNOW_STORAGE.loc[label_prev], WATERTABLE_HEIGHT.loc[label_prev]))
    elif surfacelowering and SIformation and not surfacerunoff:
        [SNOW_STORAGE.loc[label], SNOWPACK_FLUX.loc[label], MATRIX_HEIGHT.loc[label]] = \
            snowpack_flux_melt(meltrate.loc[label] * dt, SNOW_STORAGE.loc[label_prev], Sw_IRR,
                               MATRIX_HEIGHT.loc[label_prev])
        [WATERTABLE_FLUX_OUT.loc[label], WATERTABLE_HEIGHT.loc[label], DHDX.loc[label]] = \
            lateral_flux(dt, dx, dX, z, HYD_COND, SNOWPACK_FLUX.loc[label], WATERTABLE_HEIGHT.loc[label_prev])
        WATERTABLE_FLUX_IN.loc[label, x[0]] = 0
        WATERTABLE_FLUX_IN.loc[label][1:] = (WATERTABLE_FLUX_OUT.loc[label, x[:-1]]).values
        WATERTABLE_HEIGHT.loc[label].iloc[1:] = WATERTABLE_HEIGHT.loc[label].iloc[1:].values + \
                                                WATERTABLE_FLUX_OUT.loc[label].iloc[:-1].values
        if any(WATERTABLE_HEIGHT.loc[label] > MATRIX_HEIGHT.loc[label]):
            maskb = WATERTABLE_HEIGHT.loc[label] > MATRIX_HEIGHT.loc[label]
            print('water table reaching snow surface at', label, 'in grid cells', maskb.index[maskb])
            break
    ## check if below still works
    elif surfacerunoff and not surfacelowering and not surfacerunoff:
        [MATRIX_HEIGHT.loc[label], SNOW_STORAGE.loc[label], WATERTABLE_HEIGHT.loc[label],
         WATERTABLE_FLUX_OUT.loc[label],
         DHDX.loc[label]] = darcy_flow(POR_INIT, meltrate.loc[label] * dt, MATRIX_HEIGHT.loc[label_prev],
                                       SNOW_STORAGE.loc[label_prev], WATERTABLE_HEIGHT.loc[label_prev])
    else:  # check if this is still all OK/working
        [SNOW_STORAGE.loc[label], SNOWPACK_FLUX.loc[label]] = \
            snowpack_flux(meltrate.loc[label] * dt, SNOW_STORAGE.loc[label_prev], Sw_IRR,
                          MATRIX_HEIGHT.loc[label_prev])
        [WATERTABLE_FLUX_OUT.loc[label], WATERTABLE_HEIGHT.loc[label], DHDX.loc[label]] = \
            lateral_flux(dt, dx, dX, z, HYD_COND, SNOWPACK_FLUX.loc[label], WATERTABLE_HEIGHT.loc[label_prev])
        WATERTABLE_FLUX_IN.loc[label, x[0]] = 0
        WATERTABLE_FLUX_IN.loc[label][1:] = (WATERTABLE_FLUX_OUT.loc[label, x[:-1]]).values
        WATERTABLE_HEIGHT.loc[label].iloc[1:] = WATERTABLE_HEIGHT.loc[label].iloc[1:].values + \
                                                WATERTABLE_FLUX_OUT.loc[label].iloc[:-1].values

    # Calculate amount of SI formation
    if SIformation:
        # Find grid cells with positive WATERTABLE_HEIGHT
        positive_storage_mask = WATERTABLE_HEIGHT.loc[label] > 0

        if any(positive_storage_mask):
            # Find the indices where WATERTABLE_HEIGHT is positive
            positive_indices = positive_storage_mask.index[positive_storage_mask]

            # Calculate amount of refreezing for all positive indices
            time_wet = dt
            refrozen = (REFREEZE.loc[dt * ts].mm_cum - REFREEZE.loc[dt * ts - time_wet].mm_cum) / 1000

            # Subtract refrozen from WATERTABLE_HEIGHT, clipping at zero
            updated_watertable = WATERTABLE_HEIGHT.loc[label, positive_indices] - refrozen
            updated_watertable = updated_watertable.clip(lower=0)

            # Calculate the amount of refreezing
            SI.loc[label] = WATERTABLE_HEIGHT.loc[label, positive_indices] - updated_watertable
            SI.loc[label] = SI.loc[label].fillna(0)

            # Update LATERAL_STORAGE with the clipped values
            WATERTABLE_HEIGHT.loc[label, positive_indices] = updated_watertable

        # Calculate SI_TOTAL
        SI_TOTAL.loc[label] = SI.loc[label] + SI_TOTAL.loc[label_prev]

    if any(MATRIX_HEIGHT.loc[label] < 0):
        mask2 = MATRIX_HEIGHT.loc[label] < 0
        if surfacerunoff:
            print('error, matrix thickness less than 0 at', label)
            break
        else:
            print('water table reaching snow surface at', label, 'in grid cells', mask2.index)
            break

    DISCHARGE.loc[label] = WATERTABLE_FLUX_OUT.loc[label, x[-1]] + DISCHARGE.loc[label_prev]
    if surfacerunoff:
        SURFACE_RUNOFF_TOTAL.loc[label] = sum(SURFACE_RUNOFF.loc[label]) + SURFACE_RUNOFF_TOTAL.loc[label_prev]

    ## DETERMINE SIMULATION ERROR - TO CHECK IF THESE CALCULATIONS ARE CORRECT
    # if not surfacerunoff:
    #     if surfacelowering:
    #         gap = (sum(sum(meltrate.loc[:label].values)) * dt - sum(MATRIX_HEIGHT.loc[label]) - sum(WATERTABLE_HEIGHT.loc[label]) - sum(sum(SI.loc[:label, :].values)) - DISCHARGE.loc[label, 'outflow']) - sum(sum(ICESLABMELT.loc[:label, :].values))
    #         if gap != 0:
    #             ERROR.loc[label] = gap
    #         else:
    #             ERROR.loc[label] = 0
    #     else:
    #         gap = (sum(sum(meltrate.loc[:label].values)) * dt - sum(MATRIX_HEIGHT.loc[label]) - sum(WATERTABLE_HEIGHT.loc[label]) - sum(sum(SI.loc[:label, :].values)) - DISCHARGE.loc[label, 'outflow'])
    # else:
    #     gap = sum(sum(meltrate.loc[:label].values)) * dt - (sum(MATRIX_HEIGHT.loc[label]) + sum(WATERTABLE_HEIGHT.loc[label]) + sum(sum(SI.loc[:label].values)) + DISCHARGE.loc[label, 'outflow'] + SURFACE_RUNOFF_TOTAL.loc[label, 'runoff'])
    #     if gap != 0:
    #         ERROR.loc[label] = gap
    #     else:
    #         ERROR.loc[label] = 0

time_end_calc = datetime.datetime.now()
print('runtime', time_end_calc - time_start)

SNOW_STORAGE_m = water_mwe_to_height(POR_INIT, SNOW_STORAGE)
WATERTABLE_HEIGHT_m = water_mwe_to_height(POR_INIT, WATERTABLE_HEIGHT)
MATRIX_HEIGHT_m = matrix_mwe_to_height(POR_INIT, MATRIX_HEIGHT)

# update (and beautify) plotting routines
# ===============================================================================
### PLOTTING
# ===============================================================================

time_start_plot = datetime.datetime.now()
# print('start plotting at', time_start_plot)

# overview plot
plt.rcParams.update({'font.size': 28})
fig, ax = plt.subplots(3, figsize=(28, 20), gridspec_kw={'height_ratios': [3, 1, 1]})
colors = plt.cm.brg(np.linspace(0, 1, ndays))

for i, label in enumerate(t_index):
    if ndays <= 30:
        if i % nt_day == 0:
            nday = int(i / nt_day)
            if melt_inputfile == 'SEBM':
                lab = label
            else:
                lab = '%.f days' % nday
            if not WATERTABLE_HEIGHT.loc[label].isnull().all():
                if i <= nt_meltinput:
                    ax[0].plot(x, WATERTABLE_HEIGHT.loc[label], ',', linestyle='dotted', label=lab,
                               color=colors[nday - 1])
                    ax[1].plot(x, MATRIX_HEIGHT.loc[label], ',', linestyle='dotted', color=colors[nday - 1])
                    # if surfacelowering:
                    #     ax[1].plot(x, TOTAL_HEIGHT.loc[label], 'x', linestyle='dotted', color=colors[nday-1])
                else:
                    ax[0].plot(x, WATERTABLE_HEIGHT.loc[label], ',', linestyle='dashed', label=lab,
                               color=colors[nday - 1])
                    ax[1].plot(x, MATRIX_HEIGHT.loc[label], ',', linestyle='dashed', color=colors[nday - 1])
                    # if surfacelowering:
                    #     ax[1].plot(x, TOTAL_HEIGHT.loc[label], 'x', linestyle='dotted', color=colors[nday-1])
    else:
        if i % (7 * nt_day) == 0:
            nday = int(i / nt_day)
            if melt_inputfile == 'SEBM':
                lab = label
            else:
                lab = '%.f days' % nday
            if not WATERTABLE_HEIGHT.loc[label].isnull().all():
                if i <= nt_meltinput:
                    ax[0].plot(x, WATERTABLE_HEIGHT.loc[label], ',', linestyle='dotted', label=lab,
                               color=colors[nday - 1])
                    ax[1].plot(x, MATRIX_HEIGHT.loc[label], ',', linestyle='dotted', color=colors[nday - 1])
                    # if surfacelowering:
                    #     ax[1].plot(x, TOTAL_HEIGHT.loc[label], 'x', linestyle='dotted', color=colors[nday - 1])
                else:
                    ax[0].plot(x, WATERTABLE_HEIGHT.loc[label], ',', linestyle='dashed', label=lab,
                               color=colors[nday - 1])
                    ax[1].plot(x, MATRIX_HEIGHT.loc[label], ',', linestyle='dashed', color=colors[nday - 1])
                    # if surfacelowering:
                    #     ax[1].plot(x, TOTAL_HEIGHT.loc[label], 'x', linestyle='dotted', color=colors[nday - 1])

# # also plot the base value of k Q_fast, meaning  k when there is no water
# ax[1].scatter(np.arange(0, L, dx), T_RES_FAST.iloc[0], marker ='x', s = 200,
#                    color='gray', label='$k_{Qfast base}$')

ax[0].set_ylabel('water table height [m]')
if WATERTABLE_HEIGHT.max().max() > 0.5:
    ax[0].set_ylim([0, WATERTABLE_HEIGHT.max().max()])
else:
    ax[0].set_ylim([0, 0.5])
# if gentlesteepgentlegrid is True:
#     ax[0].broken_barh([(x_steep_start, L_steep)], (0, 1), facecolors='grey', alpha=0.2)

ax[0].set_xlim([0, L])
ax[0].set_ylim([0, max(TOTAL_HEIGHT.iloc[0])])
ax[0].set_xticks(np.arange(0, L + 1, nx_km * 5), np.arange(0, ((L + 1) * dx) / 1000, 5))
# ax[1].set_ylabel('$Q_{fast}$ [m$^{3}$ hr$^{-1}$]')
ax[1].set_xticks(np.arange(0, L + 1, nx_km * 5), np.arange(0, ((L + 1) * dx) / 1000, 5))
ax[1].set_ylabel('total & snow-\npack height [m]')
ax[1].set_xlim([0, L])

# ax[2].plot(np.arange(0, L, dx), z, label='height [m]')
ax[2].set_xticks(np.arange(0, L + 1, nx_km * 5), np.arange(0, ((L + 1) * dx) / 1000, 5))
ax[2].step(x, z, where='post', color='gray')
ax[2].set_xlabel('horizontal distance along slope [km]')
ax[2].set_ylabel('height [m a.s.l.]')
ax[2].set_xlim([0, L])
# ax[2].set_xlim([(-0.5) * dx, L - dx * 0.5])

# Plot water level at the following time steps:
# - if error conditions occurred: first time step with error condition (water flow not unanimously to the left)
# - last time step with melt
# - the very last time step

# ax[2].plot(np.arange(0, L, dx), HYD_HEAD.loc[label], '.', linestyle='--')
if start_error[0] == 0:  # plot final result only if no error condition observed
    ax[2].step(x, HYD_HEAD.loc[label].to_list(), where='post', linestyle='--', color='tab:orange',
               label='day {}'.format(int(len(HYD_HEAD) / int(86400 / dt) - 1)))
# plot the water depth at the last time step where there was melt. Do this only in case the model run
# consists of an initial period with melt, followed by a dry-out period
if 'meltdays' in locals():
    ax[2].step(x, HYD_HEAD.iloc[meltdays * int(86400 / dt)].to_list(),
               where='post', linestyle='--', color='tab:cyan',
               label='day {}'.format(meltdays))
if start_error[0] != 0:
    dt_err = (datetime.datetime.strptime(start_error[0], '%Y-%m-%d %H:%M:%S') - t_0) / datetime.timedelta(days=1)
    ax[2].step(np.arange(-0.5, L) * dx,
               HYD_HEAD.iloc[start_error[1]].to_list() + [HYD_HEAD.iloc[start_error[1], -1]],
               where='post', linestyle='-', color='tab:red',
               label='E@d{:.2f}'.format(dt_err), zorder=-10)
    ax[2].step(x, HYD_HEAD.iloc[start_error[1] - 1].to_list(),
               where='post', linestyle='-', color='tab:green',
               label='pre-E', zorder=-5)
# ax[2].set_ylim([0.98, 1.01])
ax[2].set_ylim([np.min(z), np.max(z) + 0.05])
# if gentlesteepgentlegrid is True:
#     ax[1].broken_barh([(x_steep_start, L_steep)], (0, 1), facecolors='grey', alpha=0.2)
fig.legend(loc='upper right', fontsize=22)
fig.suptitle('Water table height throughout modelled transect')
plt.savefig(os.path.join(output_dir, 'overview.png'))

fig, ax = plt.subplots(figsize=(28, 20))
for i, label in enumerate(t_index):
    if i % (7 * nt_day) == 0:
        nday = int(i / nt_day)
        lab = label
        ax.plot(x, WATERTABLE_HEIGHT_m.loc[label], ',', linestyle='dotted', label=lab, color=colors[nday - 1])
        ax.plot(x, MATRIX_HEIGHT_m.loc[label], ',', linestyle='dashed', color=colors[nday - 1])
ax.set_xlim([0, L])
if WATERTABLE_HEIGHT.max().max() > water_mwe_to_height(POR_INIT, gridcell_height):
    ax.set_ylim([0, WATERTABLE_HEIGHT.max().max()])
else:
    ax.set_ylim([0, water_mwe_to_height(POR_INIT, gridcell_height)])
ax.set_xticks(np.arange(0, L + 1, nx_km * 25), np.arange(0, ((L + 1) * dx) / 1000, 25))
fig.legend(loc='upper right', fontsize=22)
fig.suptitle('Water table height throughout modelled transect')
plt.savefig(os.path.join(output_dir, 'actual_heights.png'))

# if not gridinputfile:
#     fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(28, 10), gridspec_kw={'width_ratios': [1, 1]})
#     t_plot = [t_index[int(i*nt_day-1)] for i in np.arange(ndays)]
#     t_plot_labels = ['day'+i[-11:-9] for i in t_plot]
#     if not LATERAL_STORAGE.loc[label].isnull().all():
#         ax[0].plot(LATERAL_STORAGE.loc[:, x[int(len(x) / 2)]], color='b', lw=2, label='water table height')
#         ax[1].plot(LATERAL_STORAGE.loc[:, x[-1]], color='b', lw=2, label='water table height')
#         if SIformation:
#             ax[0].plot(-SI_TOTAL.loc[:, x[int(len(x) / 2)]], lw=2, color='k', ls='--', label='superimposed ice')
#             ax[1].plot(-SI_TOTAL.loc[:, x[-1]], color='k', lw=2, ls='--', label='superimposed ice')
#     ax[0].set_xticks(t_plot)
#     ax[1].set_xticks(t_plot)
#     ax[0].axhline(y=0, color='k', lw=1)
#     ax[1].axhline(y=0, color='k', lw=1)
#     ax[1].plot(LATERAL_STORAGE.loc[:, x[-1]])
#     ax[0].set_xticklabels(t_plot_labels, rotation=65)
#     ax[1].set_xticklabels(t_plot_labels, rotation=65)
#     ax[0].set_ylabel('water table height [m]')
#     ax[0].set_xlabel('middle grid cell of transect', labelpad=20)
#     ax[1].set_xlabel('final grid cell of transect', labelpad=20)
#     ax[0].set_xlim([LATERAL_STORAGE.index[0], LATERAL_STORAGE.index[-1]])
#     ax[1].set_xlim([LATERAL_STORAGE.index[0], LATERAL_STORAGE.index[-1]])
#     ax[0].set_ylim([-0.05, 0.2])
#     ax[1].set_ylim([-0.05, 0.2])
#     ax[0].axvline(x=nt_meltinput, color='k', lw=1, ls=':')
#     ax[1].axvline(x=nt_meltinput, color='k', lw=1, ls=':')
#     ax[0].grid(True)
#     ax[1].grid(True)
#     for tick in ax[0].get_xticklabels():
#         tick.set_horizontalalignment('right')
#     for tick in ax[1].get_xticklabels():
#         tick.set_horizontalalignment('right')
#     ax[0].legend(loc='upper right')
#     ax[1].legend(loc='upper right')
#     fig.suptitle('Water table height over time')
#     plt.subplots_adjust(bottom=0.25)
#     plt.savefig(os.path.join(output_dir, filename[:-5] + '_gridcells.png'))


# single timestep-plots for GIF
if gif:
    [bins, V_min, V_max] = create_bins(WATERTABLE_HEIGHT, t_index, t, dx)
    if V_min > V_max:
        le = V_max
        V_max = V_min
        V_min = le
    for j, ts in enumerate(t_index):
        label = t_index[j]
        if j % nt_day == 0:
            fig, ax = plt.subplots(2, figsize=(16, 9), gridspec_kw={'height_ratios': [3, 1]})
            if not WATERTABLE_HEIGHT.loc[label].isnull().all():
                ax[0].scatter(np.arange(0, L, dx), WATERTABLE_HEIGHT.loc[label], c=dx / WATERTABLE_HEIGHT.loc[label],
                              label='%.f days' % (j / nt_day), cmap='viridis', marker='.')
                # for i in range(len(x)-1):
            #     left = x[i]-0.5*dx
            #     right = left+dx
            #     vel_color = define_color(V_min, V_max, dx/T_RES_FAST.loc[label, x[i]]*3600)[0]
            #     ax[0].fill_between([left, right], Q_FAST_STORAGE.loc[label, x[i]], color=vel_color, alpha=0.7)
            ax[0].set_xlim(0, L - dx)
            ax[0].set_ylim(0, WATERTABLE_HEIGHT.max().max())
            # ax[0].set_ylim(0, 1)
            ax[0].set_ylabel('water table height [m]')
            ax[0].annotate('date: %s' % label[:10], xy=(0.5, 0.9 * WATERTABLE_HEIGHT.max().max()), xycoords='data',
                           fontsize='large', fontweight='bold')
            melt = sum(meltrate[int(j - nt_day):j]) * 3600
            ax[0].annotate('total melt input the last 24 hrs: %.3f mm' % melt,
                           xy=(0.5, 0.8 * WATERTABLE_HEIGHT.max().max()),
                           xycoords='data', fontsize='large')
            ax[0].annotate('total discharge: %.3f m' % DISCHARGE.loc[label][0],
                           xy=(0.5, 0.7 * WATERTABLE_HEIGHT.max().max()),
                           xycoords='data', fontsize='large')
            # ax[1].plot(np.arange(0, L, dx), z )
            ax[1].plot(np.arange(0, L, dx), z, label='height [m]')
            ax[1].set_xlabel('horizontal distance along slope [m]')
            ax[1].set_ylabel('height [m]')
            ax[1].set_xlim([0, L - dx])
            ax[1].set_ylim([0.8, 1])
            # if ts <= nt_meltinput:
            #     ax.annotate('melt rate: %.2f mm/hr' % (meltrate[j]*3600), xy=(0.5, 0.8), xycoords='data', fontsize='large', color='red')
            # else:
            #     ax.annotate('no melt input', xy=(0.5, 0.8), xycoords='data', fontsize='large')
            # if gentlesteepgentlegrid is True:
            #     ax[0].broken_barh([(x_steep_start, L_steep)], (0, 1), facecolors='grey', alpha=0.2)
            #     ax[1].broken_barh([(x_steep_start, L_steep)], (0, 1), facecolors='grey', alpha=0.2)
            plt.subplots_adjust(bottom=0.1, left=0.07, right=0.92, top=0.95)
            # plt.colorbar(define_color(V_min, V_max, dx / T_RES_FAST.loc[label, x[i]] * 3600)[-1],
            #              label='apparent lateral flow velocity [m/hr]', ax=ax, shrink=0.8, location='right')
            plt.savefig(os.path.join(output_dir, 'gif', str(label[:10] + '_' + label[11:13] + 'h00.png')))
            plt.close()

time_end_plot = datetime.datetime.now()
# print('finished plotting at', time_end_plot, 'start writing output')
print('elapsed time for plotting', time_end_plot - time_start_plot)

# ===============================================================================
### WRITING OUTPUT FILES
# ===============================================================================

# Writing output dataframes to one Excel-file (slow):
# with pd.ExcelWriter(os.path.join(output_dir, filename)) as writer:
#     [globals()[s].to_excel(writer, sheet_name=s) for s in outputnames]

# Writing output dataframes to individual .csv-files
[globals()[s].to_csv(os.path.join(output_dir, f'{s}.csv'), index=True) for s in outputnames]

time_end_output = datetime.datetime.now()
print('elapsed time for writing output', time_end_output - time_end_plot, '\ntime finished', datetime.datetime.now())

# # ===============================================================================
# ### CREATING GIF
# # ===============================================================================

if gif:
    print('start making gif at', time_end_plot)

    make_gif((os.path.join(output_dir, 'gif')), output_dir, filename[0:-5] + '.gif')

    time_end_gif = datetime.datetime.now()
    print('elapsed time for making gif', time_end_gif - time_end_plot)
