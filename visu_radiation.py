import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys

fig, axes = plt.subplots(1,2,subplot_kw={'projection': 'polar'})

filnm = 'ck_radpattern.txt'
names = ['STmode','deg','power']
df = pd.read_csv(filnm, header=None, names=names, delim_whitespace=True)

st_modes = df['STmode'].unique()

df_mode_S  = df[df['STmode'] == 1]
list_deg_S = df_mode_S['deg'].tolist()
list_rad_S = np.deg2rad(list_deg_S)
list_power_S = df_mode_S['power'].tolist()

df_mode_T  = df[df['STmode'] == 2]
list_deg_T = df_mode_T['deg'].tolist()
list_rad_T = np.deg2rad(list_deg_T)
list_power_T = df_mode_T['power'].tolist()

max_S = np.max(list_power_S[0:360])
max_T = np.max(list_power_T)
max_radia = max(max_S,max_T)

axes[0].plot(list_rad_S[0:360],list_power_S[0:360])
axes[0].set_title(f'Spheroidal mode')
axes[0].set_theta_offset(np.pi/2)
axes[0].set_theta_direction(-1)
axes[0].set_ylim([0,max_radia])

axes[1].plot(list_rad_T,list_power_T)
axes[1].set_title(f'Toroidal mode')
axes[1].set_theta_offset(np.pi/2)
axes[1].set_theta_direction(-1)
axes[1].set_ylim([0,max_radia])

plt.show()

#STnum = 0
#for idx, st_mode in enumerate(st_modes):
#    df_mode = df[df['STmode'] == st_mode]
#    list_deg   = df_mode['deg'].tolist()
#    list_rad = np.deg2rad(list_deg)
#    list_power = df_mode['power'].tolist()
 

#    axes[STnum].plot(list_rad,list_power)
#    axes[STnum].set_title(f'STmode {st_mode}')
#    axes[STnum].set_theta_offset(np.pi/2)
#    axes[STnum].set_theta_direction(-1)

#    STnum = STnum + 1

#plt.show()
