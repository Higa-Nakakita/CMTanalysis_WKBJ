import numpy as np
from obspy.geodetics import gps2dist_azimuth
from scipy.optimize import brentq
import sys
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.patches as patches
import os
import pandas as pd

m2km = 1e-3
km2m = 1e3
rad2deg = 180/np.pi 
deg2rad = np.pi/180

def sph2cart(latlon_list):
    # 3D Spherical coordi (theta,phi) -> 3D Euclid coordi (x,y,z)
    theta,phi = latlon_list[0],latlon_list[1]
    lat = theta * deg2rad
    lon = phi * deg2rad
    x = np.cos(lat) * np.cos(lon)
    y = np.cos(lat) * np.sin(lon)
    z = np.sin(lat)
    return np.array([x, y, z])

def cart2sph(cart_vec):
    # 3D Euclid coordi (x,y,z) -> 3D Spherical coordi (theta,phi) 
    x,y,z = cart_vec[0],cart_vec[1],cart_vec[2]
    phi = np.arctan2(y,x)
    r_xy = np.sqrt(x**2+y**2)
    theta = np.arctan2(z,r_xy)
    theta,phi = theta*rad2deg, phi*rad2deg
    return np.array([theta,phi])

def Slerp(A_cart,B_cart,t):
    # Spherial linear interpolation
    # https://swkagami.hatenablog.com/entry/lie_08slerp
    # Represent a point on the gcp by internal division
    if (t<0 or t>1 ):
        print('t must be 0<=t<=1')
        sys.exit()
    omega = np.arccos(np.dot(A_cart,B_cart))
    q1_vec = np.sin((1-t)*omega)/np.sin(omega) * A_cart
    q2_vec = np.sin(t*omega)/np.sin(omega) * B_cart
    q_vec  = q1_vec + q2_vec
    abs_q  = np.sqrt(q_vec[0]**2+q_vec[1]**2+q_vec[2]**2)
    point  = q_vec/abs_q
    return point

def rd_rcvs_info(filepath):
    receivers = []
    with open(filepath, 'r') as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) == 3:
                station, latitude, longitude = parts
                receivers.append({
                    'station': station,
                    'latitude': float(latitude),
                    'longitude': float(longitude)
                })
    return receivers


# Make Grids
d_grid = 0.5
lat_min, lat_max = 19.75, 49.75
lon_min, lon_max = 119.75, 149.75
lat_grids = np.arange(lat_min,lat_max+d_grid,d_grid)
lon_grids = np.arange(lon_min,lon_max+d_grid,d_grid)

## Set params for drawing figure
#fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()}) #,figsize=(15,15))
#ax.coastlines(resolution='10m')
#ax.set_extent([115,155,15,55])
#gl = ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=True,alpha=0.5)
#gl.top_labels = False
#gl.right_labels = False
# draw grids
#crossed_grids = []
#for lat in lat_grids:
#    for lon in lon_grids:
#        rect = patches.Rectangle((lon,lat),d_grid,d_grid,edgecolor='black', facecolor='none', linewidth=0.5, transform=ccrs.PlateCarree())
#        ax.add_patch(rect)

# source croodinates (tranform spherical -> cartesian) 
names = ['lat','lon','dep']
df = pd.read_csv('src_coordinates',names=names,delim_whitespace=True)
slat,slon,sdep = df['lat'][0],df['lon'][0],df['dep'][0]

#Prevent rounding 2025 0721
slat_str = f"{df['lat'][0]:.2f}"
slon_str = f"{df['lon'][0]:.2f}"

src_latlon = np.array([slat,slon])
src_cart = sph2cart(src_latlon)
#src_latlon_conf = cart2sph(src_cart)
#print(src_latlon_conf)
#ax.plot(slon,slat,marker='o',ms=2,c='k',linewidth=0,label='source')

fil = 'output_'+f'{slat_str}_{slon_str}'+'.txt'
if (os.path.isfile(fil)):
    os.remove(fil)
#with open(fil,mode='w') as f:
#    f.write(f'{slat:.3f} {slon:.3f}\n')

# receiver coordinates (transform spherical -> cartersian)
Fnet_stas = '/work94A/higa/Fnet_data/Fnet_stalist_ref'
#Fnet_stas = 'rcv_coordinates'
rcvs_info = rd_rcvs_info(Fnet_stas)
for rcv in rcvs_info:
    rlat,rlon,sta = rcv['latitude'],rcv['longitude'],rcv['station']
    rcv_latlon = np.array([rlat,rlon])
    #print(rcv_latlon)
    rcv_cart = sph2cart(rcv_latlon)
    #print(rcv_cart)
    #rcv_latlon_conf = cart2sph(rcv_cart)
    #print(rcv_latlon_conf)

    cross_points = []
    # Find intersect ponint of latitude
    for target_lat in lat_grids:
        #print(target_lat)
        def diff_lat(t):
            gcp_cart = Slerp(src_cart,rcv_cart,t)
            lat,lon = cart2sph(gcp_cart)
            diff = lat - target_lat
            return diff
        try:
            t_cross = brentq(diff_lat,0,1) # Find t with satisifie diff=0 within 0<=t<=1
            gcp_cart = Slerp(src_cart,rcv_cart,t_cross)
            lat,lon = cart2sph(gcp_cart)
            cross_points.append([lat,lon])
        except ValueError :
            # When diff_lat(0)*diff_lat(1) > 0
            pass
    
    for target_lon in lon_grids:
        #print(target_lon)
        def diff_lon(t):
            gcp_cart = Slerp(src_cart,rcv_cart,t)
            lat,lon = cart2sph(gcp_cart)
            diff = lon - target_lon
            return diff
        try:
            t_cross = brentq(diff_lon,0,1) # Find t with satisifie diff=0 within 0<=t<=1
            gcp_cart = Slerp(src_cart,rcv_cart,t_cross)
            lat,lon = cart2sph(gcp_cart)
            cross_points.append([lat,lon])
        except ValueError :
            # When diff_lat(0)*diff_lat(1) > 0
            pass
    
    # add src and rcv coordinates
    cross_points.append(src_latlon)
    cross_points.append(rcv_latlon)
    
    #print(cross_points)
    #print(len(cross_points))
    num_cross = len(cross_points)
    
    # Avoid overlap point
    # I think overlap cannot occur unless grid line and ray path do completely overlap
    unique_points = []
    for i in range(num_cross):
        lat_ref = cross_points[i][0]
        lon_ref = cross_points[i][1]
        k = 0
        for j in range(num_cross):
            lat = cross_points[j][0]
            lon = cross_points[j][1]
            if lat_ref==lat and lon_ref==lon :
                 k = k + 1
                 if k== 2 :
                    print('Detect overlap at',lat,lon)
                    #print(k)
                    #print('stop??')
                    #sys.exit()
    
    # Sort cross points by epicentral dist
    dist_list = []
    for i in range(num_cross):
        lat = cross_points[i][0]
        lon = cross_points[i][1]
        dist,_,_ = gps2dist_azimuth(slat,slon,lat,lon)
        dist_list.append(dist)
    
    # Make list
    idxsort_list = np.argsort(dist_list)
    cross_points_sort = []
    cross_points_sort = np.zeros((4,num_cross))
    j = 0
    for i in idxsort_list:
        lat = cross_points[i][0]
        lon = cross_points[i][1]
        dist = dist_list[i]
        dist = dist * m2km
        cross_points_sort[0][j] = lat
        cross_points_sort[1][j] = lon
        cross_points_sort[3][j] = dist
        j = j + 1
    
    # add d_delta info
    for i in range(num_cross-1):
        lat1, lon1 = cross_points_sort[0][i],cross_points_sort[1][i]
        lat2, lon2 = cross_points_sort[0][i+1],cross_points_sort[1][i+1]
        dist,_,_ = gps2dist_azimuth(lat1,lon1,lat2,lon2)
        dist = dist * m2km
        cross_points_sort[2][i] = dist
    #print(cross_points_sort) 
    
    
    # Draw passed grids
    for i in range(num_cross-1):
        lat1, lon1 = cross_points_sort[0][i],cross_points_sort[1][i]
        lat2, lon2 = cross_points_sort[0][i+1],cross_points_sort[1][i+1]
        lat_mid = (lat1 + lat2)*0.5
        lon_mid = (lon1 + lon2)*0.5
        lon_g = int((lon_mid - lon_min)/d_grid) * d_grid + lon_min
        lat_g = int((lat_mid - lat_min)/d_grid) * d_grid + lat_min
        #print(lon_g,lat_g)
        passed_grids = patches.Rectangle((lon_g,lat_g),d_grid,d_grid,edgecolor='red',facecolor='none',linewidth=2,transform=ccrs.PlateCarree())
        #ax.add_patch(passed_grids)
    
    
    # draw great circle path by divided points
    lats, lons = cross_points_sort[0][:], cross_points_sort[1][:]
    #ax.plot(lons,lats,marker='o',ms=1,c='b',linewidth=0.5)
    
    #ax.plot(rlon,rlat,marker='o',ms=2,c='k',linewidth=0)#,label=f'{sta}')
    #ax.legend()
    
    total_dist = 0
    for i in range(num_cross-1):
        lat1, lon1 = cross_points_sort[0][i],cross_points_sort[1][i]
        lat2, lon2 = cross_points_sort[0][i+1],cross_points_sort[1][i+1]
        dist = cross_points_sort[2][i]
        #print(f"{lat1:.3f}, {lon1:.3f} -> {lat2:.3f}, {lon2:.3f} : {dist :.3f} km")
        total_dist += dist
    
    #print('sum of div: ',total_dist,' (km)')
    dist,_,_ = gps2dist_azimuth(slat,slon,rlat,rlon)
    conf_dist = dist * m2km
    #print('gps2dist_azimuth: ',conf_dist,' (km)')
    
    
    
    with open(fil,mode='a') as f:
       for i in range(num_cross-1):
            lat1, lon1 = cross_points_sort[0][i],cross_points_sort[1][i]
            lat2, lon2 = cross_points_sort[0][i+1],cross_points_sort[1][i+1]
            lat_mid = (lat1 + lat2)*0.5
            lon_mid = (lon1 + lon2)*0.5
            lon_g = int((lon_mid - lon_min)/d_grid) * d_grid + lon_min
            lat_g = int((lat_mid - lat_min)/d_grid) * d_grid + lat_min
            #print(lon_g,lon_g+d_grid,lat_g,lat_g+d_grid)
            lon_w,lon_e,lat_s,lat_n = lon_g,lon_g+d_grid,lat_g,lat_g+d_grid # west,east,south,north
            dist = cross_points_sort[2][i]
            f.write(f"{lat1:.3f} {lon1:.3f} {lat2:.3f} {lon2:.3f} {lat_s:.3f} {lat_n:.3f} {lon_w:.3f} {lon_e:.3f} {dist :.3f}\n")

#plt.show()
