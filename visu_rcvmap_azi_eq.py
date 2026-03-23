import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from obspy.imaging.mopad_wrapper import beach
import sys

# =============================================
# visu_rcvmap_azi_eq.py just draw map of didtribution of stations usig in inversion 
# Input file  : rcv_coordinates
# A prioiri info: src coordinates
# Outputs file: map_rcvs_distri.pdf
# ====== Higa (June 16 2024) ===============

def draw_stamap2():
 
    scl='110m'
    land = cfeature.NaturalEarthFeature('physical','land',scl,edgecolor='face',facecolor='w')
    ocean = cfeature.NaturalEarthFeature('physical','ocean',scl,edgecolor='face',facecolor='0.95')
    countries = cfeature.NaturalEarthFeature('cultural','admin_0_countries',scl,edgecolor='gray',facecolor='none')


    
    names_src = ["srclat","srclon","srcdep"]
    df_src = pd.read_csv("src_coordinates", header=None,names=names_src,delim_whitespace=True)
    src_lat = float(df_src['srclat'])
    src_lon = float(df_src['srclon'])

    names = ['rcv','lat','lon']
    df_rcv = pd.read_csv('rcv_coordinates',header=None,names=names,delim_whitespace=True)
    rcvs = []

    for i in range(len(df_rcv)):
        rcvs.append([df_rcv.at[i,"rcv"],df_rcv.at[i,"lat"],df_rcv.at[i,"lon"]])

    x_list = []
    y_list = []
    rcv_list = []
    #cate_list = []
    for rcvnm,rlat, rlon in rcvs:
        x_list.append(rlon)
        y_list.append(rlat)
        rcv_list.append(rcvnm)

    #print(x_list,y_list)
    x_max, x_min = max(x_list),min(x_list)
    y_max, y_min = max(y_list),min(y_list)
    dx, dy = 5, 5
    x_max, x_min = int(x_max + dx), int(x_min - dx)
    y_max, y_min = int(y_max + dy), int(y_min - dy)
    print(x_min,x_max,y_min,y_max)

    # Draw map
    fsz = 10
    plt.rcParams['font.size'] = fsz
    gs = gridspec.GridSpec(1,1)
    ax = [plt.subplot(gs[0,0],projection=ccrs.PlateCarree())]
    #ax = [plt.subplot(gs[0,0],projection=ccrs.AzimuthalEquidistant(central_longitude=src_lon, central_latitude=src_lat))]

    gl = ax[0].gridlines(draw_labels=True, linestyle='-',color='0.8',linewidth=0.1)
    gl.top_labels = False
    gl.right_labels = False

    ax[0].set_extent([x_min,x_max,y_min,y_max],ccrs.PlateCarree())
    #ax[0].set_xticks(list(range(x_min,x_max+dx,dx)),crs=ccrs.PlateCarree())
    #ax[0].set_yticks(list(range(y_min,y_max+dy,dy)),crs=ccrs.PlateCarree())
    ax[0].add_feature(land)
    ax[0].add_feature(ocean)
    ax[0].add_feature(countries)
    ax[0].set_title('The Distribution of the Stations',fontsize=15)
    ax[0].set_xlabel('Longitude (degree)',fontsize=10)
    ax[0].set_ylabel('Latitude (degree)',fontsize=10)
    ax[0].gridlines(linestyle='-',color='0.8',linewidth=0.1)

    # Draw statioms
#    print(sta)
    for rcvnm,rlat,rlon in rcvs:
        ax[0].plot(float(rlon),float(rlat),marker='v',color='b',mec='k',ms=10,transform=ccrs.PlateCarree())
        ax[0].text(float(rlon),float(rlat+2),f"{rcvnm}",horizontalalignment='center',verticalalignment='top',color='b',fontsize=13,transform=ccrs.PlateCarree())


    # Define mt_orig  & cent_orig but it should be input----------------
    mt_orig = [1 , 0 , -1 , 0, 0 , 0]
    cent_orig = [src_lat, src_lon, 50]
    # Get beach orig
    beachorig1 = beach(mt_orig,xy=(cent_orig[1],cent_orig[0]),width=1.5,\
                        linewidth=0.8 , facecolor='b') 

    #beachorig1 = beach(mt_orig,xy=(0,0),width=1.5,\
    #                    linewidth=0.8 , facecolor='b')  

    # Draw beachball on the map
    ax[0].add_collection(beachorig1)
    #ax[0].set_aspect('equal')
    #ax[0].set_xlim(X_lim)
    #ax[0].set_ylim(Y_lim)

    plt.show()
    #plt.savefig("map_rcvs_distri.pdf",dpi=600)


draw_stamap2()
