import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
from obspy.imaging.beachball import beach
module_dir='/home/higa/Store_House/visu_CMTresult'
sys.path.append(module_dir)
from plot_srcrcv_relation import plot_srcrcv_rela
from Useful_funcs_cmt import min2deci, Mxx2Mrr, calc_M0Mw, calc_NDC

## READ FILES ##
# obatain initial source location
ini_names = ['lat','lon','dep']
ini_df = pd.read_csv('../../../src_coordinates',names=ini_names,delim_whitespace=True)
ini_lat = float(ini_df['lat'])
ini_lon = float(ini_df['lon'])
ini_dep = float(ini_df['dep'])

# obtain ensemble
names = ['dt','dlat','dlon','ddep','Mrr','Mtt','Mpp','Mrt','Mrp','Mtp','scalar','halfd','mis','misWAV']
df = pd.read_csv('all_params.dat',names=names,delim_whitespace=True)
# obtain optimal model
best_row = df['mis'].idxmin()
r = df.loc[best_row]
dt_best = r['dt']
dlat_best = r['dlat']
dlon_best = r['dlon']
ddep_best = r['ddep']
Mrr_best = r['Mrr']
Mtt_best = r['Mtt']
Mpp_best = r['Mpp']
Mrt_best = r['Mrt']
Mrp_best = r['Mrp']
Mtp_best = r['Mtp']
scalar_best = r['scalar']
halfd_best = r['halfd']
MT_best = [Mrr_best,Mtt_best,Mpp_best,Mrt_best,Mrp_best,Mtp_best]
print('----- Best result:')
print('dt=',dt_best,'dlat=',dlat_best,'dlon=',dlon_best,'ddep=',ddep_best)
print('MT=',MT_best)
print('halfd=',halfd_best)
print('misfit=',r['mis'])
print('----- ')


# convert
best_lat = dlat_best + ini_lat
best_lon = dlon_best + ini_lon
best_dep = ddep_best + ini_dep


dts = df['dt']
dlats = df['dlat']
dlons = df['dlon']
ddeps = df['ddep']
Mrrs = df['Mrr']
Mtts = df['Mtt']
Mpps = df['Mpp']
Mrts = df['Mrt']
Mrps = df['Mrp']
Mtps = df['Mtp']
scalars = df['scalar']
halfds = df['halfd']
misfits = df['mis']
array_len = len(dts)
count = np.arange(1,array_len+1)

M0s, Mws = np.zeros(array_len), np.zeros(array_len)
for i in range(array_len):
    MT = Mrrs[i],Mtts[i],Mpps[i],Mrts[i],Mrps[i],Mtps[i]
    exp = scalars[i] - 7 # Since calc_M0Mw's argument is Nm, dyncm -> Nm
    M0s[i],Mws[i] = calc_M0Mw(MT,exp)

CMT_best = Mrr_best, Mtt_best, Mpp_best, Mrt_best, Mrp_best, Mtp_best
exp_best = scalar_best - 7
M0_best, Mw_best = calc_M0Mw(CMT_best,exp_best)

ISO_best,CLVD_best,DC_best,diff_best,stress_best = calc_NDC([Mrr_best,Mtt_best,Mpp_best,Mrt_best,Mrp_best,Mtp_best])

### PLOT ##
fig, axes = plt.subplots(3,3,figsize=(8,8))
beach1 = beach(MT_best,xy=(0,0.1),width=1.8,linewidth=0.8,facecolor='r')
axes[0][0].add_collection(beach1)
axes[0][0].set_aspect('equal')
axes[0][0].set_xlim([-1,1])
axes[0][0].set_ylim([-1,1])
axes[0][0].text(-1,-1.1,f'Mw: {Mw_best}',size=12)
axes[0][0].text(-1,-1.3,f'ISO: {ISO_best}  CLVD: {CLVD_best} DC: {DC_best}',size=12)
#axes[0][0].text(-0.6,-1.2,f'NDC: {NDC_best}',size=12)
axes[0][0].text(-1.3,-1.5,f'Lat: {np.round(best_lat,2)}(deg)  Lon: {np.round(best_lon,2)}(deg)',size=12)
#axes[0][0].text(-0.6,-1.6,f'Lon: {np.round(best_lat,2)}',size=12)
axes[0][0].text(-1.3,-1.8,f'Dep: {np.round(best_dep,2)}(km)  Halfd: {np.round(halfd_best,2)}(s)',size=12)
axes[0][0].text(-1.3,-2.1,f'Time shift: {np.round(dt_best,2)}(s)',size=12)
#axes[0][0].text(-0.6,-2.0,f'Halfd: {np.round(halfd_best,2)}',size=12)
axes[0][0].axis('off')

axes[0][1].plot(count,Mrrs,marker='o',color='c',linestyle='None',ms=0.5)
axes[0][1].plot(array_len+3,Mrr_best,marker='v',color='c',label='Mrr',linestyle='None',ms=2.5)
axes[0][1].plot(count,Mtts,marker='o',color='b',linestyle='None',ms=0.5)
axes[0][1].plot(array_len+3,Mtt_best,marker='v',color='b',label='M\u03b8\u03b8',linestyle='None',ms=2.5)
axes[0][1].plot(count,Mpps,marker='o',color='y',linestyle='None',ms=0.5)
axes[0][1].plot(array_len+3,Mpp_best,marker='v',color='y',label='M\u03c6\u03c6',linestyle='None',ms=2.5)
axes[0][1].plot(count,Mrts,marker='o',color='k',linestyle='None',ms=0.5)
axes[0][1].plot(array_len+3,Mrt_best,marker='v',color='k',label='Mr\u03b8',linestyle='None',ms=2.5)
axes[0][1].plot(count,Mrps,marker='o',color='g',linestyle='None',ms=0.5)
axes[0][1].plot(array_len+3,Mrp_best,marker='v',color='g',label='Mr\u03c6',linestyle='None',ms=2.5)
axes[0][1].plot(count,Mtps,marker='o',color='r',linestyle='None',ms=0.5)
axes[0][1].plot(array_len+3,Mtp_best,marker='v',color='r',label='M\u03b8\u03c6',linestyle='None',ms=2.5)
axes[0][1].set_title('Transition of moment tensor components',fontsize=10)
axes[0][1].set_xlabel('the number of models generated',fontsize=10)
axes[0][1].set_ylabel('moment tensor (dyn cm)',fontsize=10)
axes[0][1].legend(loc='lower left')

import matplotlib.colors  as mcolors
from matplotlib.colors import PowerNorm
norm = PowerNorm(gamma=0.5, vmin=misfits.min(), vmax=misfits.mean())
#norm = PowerNorm(gamma=0.5, vmin=misfits.min(), vmax=np.percentile(misfits,60))
#norm = mcolors.LogNorm(vmin=misfits.min(), vmax=misfits.mean()*0.02)
axes[0][2].scatter(count,np.log10(misfits),c=misfits,s=0.5,cmap='viridis',norm=norm)
axes[0][2].set_title('Transition of misfit',fontsize=10)
axes[0][2].set_xlabel('the number of models generated',fontsize=10)
axes[0][2].set_ylabel('log10(misfit)',fontsize=10)

axes[1][0].scatter(ini_lon+dlons,ini_lat+dlats,c=misfits,s=0.5,cmap='viridis',norm=norm)
axes[1][0].plot(best_lon,best_lat,ms=6,c='r',marker='v',mew=0.5,mec='w')
#axes[1][0].plot(jma_lon,jma_lat,ms=6,c='b',marker='v',mew=0.5,mec='k')
#axes[1][0].plot(fnet_lon,fnet_lat,ms=6,c='y',marker='v',mew=0.5,mec='k')
#axes[1][0].plot(gcmt_lon,gcmt_lat,ms=6,c='m',marker='v',mew=0.5,mec='w')
axes[1][0].plot(ini_lon,ini_lat,ms=6,c='w',marker='v',mew=0.5,mec='k')
axes[1][0].set_xlabel('longitude (deg)')
axes[1][0].set_ylabel('latitude (deg)')
axes[1][0].set_title('longitude-latitude')
ylim = np.round(axes[1][0].get_ylim(),1)
yticks = np.arange(ylim[0],ylim[1],0.25)
xlim = np.round(axes[1][0].get_xlim(),1)
xticks = np.arange(xlim[0],xlim[1],0.25)
axes[1][0].tick_params(axis='x',rotation=45)
axes[1][0].set_xticks(xticks)
axes[1][0].set_yticks(yticks)
axes[1][0].grid(True)
axes[1][0].set_aspect('equal')

#axes[1][1].scatter(dlats,ddeps,c=misfits,s=0.5,cmap='seismic',norm=norm)
#axes[1][1].scatter(dlats,ddeps,c=misfits,s=0.5,cmap='inferno',norm=norm)
axes[1][1].scatter(ini_lat+dlats,ini_dep+ddeps,c=misfits,s=0.5,cmap='viridis',norm=norm)
axes[1][1].plot(best_lat,best_dep,ms=6,c='r',marker='v',mew=0.5,mec='w')
#axes[1][1].plot(jma_lat,jma_dep,ms=6,c='b',marker='v',mew=0.5,mec='k')
#axes[1][1].plot(fnet_lat,fnet_dep,ms=6,c='y',marker='v',mew=0.5,mec='k')
#axes[1][1].plot(gcmt_lat,gcmt_dep,ms=6,c='m',marker='v',mew=0.5,mec='w')
axes[1][1].plot(ini_lat,ini_dep,ms=6,c='w',marker='v',mew=0.5,mec='k')
#axes[1][1].plot(true_lat,true_dep,ms=5,c='b',marker='v')
xlim = np.round(axes[1][1].get_xlim(),1)
xticks = np.arange(xlim[0],xlim[1],0.25)
axes[1][1].set_xticks(xticks)
axes[1][1].invert_yaxis()
axes[1][1].set_xlabel('latitude (deg)')
axes[1][1].set_ylabel('depth (km)')
axes[1][1].set_title('latitude-depth')
axes[1][1].grid(True)
axes[1][1].set_box_aspect(aspect=0.6)

#axes[1][2].scatter(dlons,ddeps,c=misfits,s=0.5,cmap='gist_rainbow',norm=norm)
#axes[1][2].scatter(dlons,ddeps,c=misfits,s=0.5,cmap='magma',norm=norm)
axes[1][2].scatter(ini_lon+dlons,ini_dep+ddeps,c=misfits,s=0.5,cmap='viridis',norm=norm)
axes[1][2].plot(best_lon,best_dep,ms=6,c='r',marker='v',mew=0.5,mec='w')
#axes[1][2].plot(jma_lon,jma_dep,ms=6,c='b',marker='v',mew=0.5,mec='k')
#axes[1][2].plot(fnet_lon,fnet_dep,ms=6,c='y',marker='v',mew=0.5,mec='k')
#axes[1][2].plot(gcmt_lon,gcmt_dep,ms=6,c='m',marker='v',mew=0.5,mec='w')
axes[1][2].plot(ini_lon,ini_dep,ms=6,c='w',marker='v',mew=0.5,mec='k')
#axes[1][2].plot(true_lon,true_dep,ms=5,c='b',marker='v')
xlim = np.round(axes[1][2].get_xlim(),1)
xticks = np.arange(xlim[0],xlim[1],0.25)
axes[1][2].set_xticks(xticks)
axes[1][2].invert_yaxis()
axes[1][2].set_xlabel('longitude (deg)')
axes[1][2].set_ylabel('depth (km)')
axes[1][2].set_title('longitude-depth')
axes[1][2].grid(True)
axes[1][2].set_box_aspect(aspect=0.6)


#sc = axes[2][0].scatter(dts+halfds,count,c=misfits,s=0.5,cmap='viridis',norm=norm)
#axes[2][0].plot(dt_best+halfd_best,array_len+1,c='r',ms=6,marker='v',mew=0.5,mec='w')
#axes[2][0].plot(true_dt,array_len+1,c='b',marker='v')
#axes[2][0].plot(jma_dt,array_len+1,c='b',ms=6,marker='v',mew=0.5,mec='k')
#axes[2][0].plot(gcmt_dt,array_len+1,c='m',ms=6,marker='v',mew=0.5,mec='w')
#axes[2][0].set_title('centroid time (s)')
#axes[2][0].set_ylabel('the number of models generated',fontsize=8)
#axes[2][0].set_xlabel('sec',fontsize=8)
#axes[2][0].set_box_aspect(aspect=0.8)

sc= axes[2][0].scatter(dts,count,c=misfits,s=0.5,cmap='viridis',norm=norm)
axes[2][0].plot(dt_best,array_len+1,c='r',ms=6,marker='v',mew=0.5,mec='w')
axes[2][0].set_title('time shift (s)')
axes[2][0].set_ylabel('the number of models generated',fontsize=8)
axes[2][0].set_xlabel('sec',fontsize=8)


sc = axes[2][1].scatter(np.log10(M0s),count,c=misfits,s=0.5,cmap='viridis',norm=norm)
axes[2][1].plot(np.log10(M0_best),array_len+1,ms=6,c='r',marker='v',mew=0.5,mec='w')
#axes[2][1].plot(np.log10(true_M0),array_len+1,ms=5,c='b',marker='v')
#axes[2][1].plot(np.log10(jma_M0),array_len+1,ms=6,c='b',marker='v',mew=0.5,mec='k')
#axes[2][1].plot(np.log10(fnet_M0),array_len+1,ms=6,c='y',marker='v',mew=0.5,mec='k')
#axes[2][1].plot(np.log10(gcmt_M0),array_len+1,ms=6,c='m',marker='v',mew=0.5,mec='k')
axes[2][1].set_title('seismic moment (dyn cm)')
axes[2][1].set_ylabel('the number of models generated',fontsize=8)
axes[2][1].set_xlabel('(dyn cm)',fontsize=8)
#axes[2][1].set_aspect('equal')

#sc = axes[2][2].scatter(halfds,count,c=misfits,s=0.5,cmap='gist_rainbow',norm=norm)
sc = axes[2][2].scatter(halfds,count,c=misfits,s=0.5,cmap='viridis',norm=norm)
axes[2][2].plot(halfd_best,array_len+1,ms=6,c='r',marker='v',mew=0.5,mec='w')
#axes[2][2].plot(true_halfd,array_len+1,ms=5,c='b',marker='v')
#axes[2][2].plot(gcmt_halfd,array_len+1,ms=6,c='m',marker='v',mew=0.5,mec='k')
axes[2][2].set_title('half duration (s)')
axes[2][2].set_ylabel('the number of models generated',fontsize=8)
axes[2][2].set_xlabel('(sec)',fontsize=8)
#axes[2][2].set_aspect(aspect=1)

#plt.colorbar(sc,ax=axes[2][2], label='misfit')
plt.tight_layout()
plt.show()


plot_srcrcv_rela(MT_best)

