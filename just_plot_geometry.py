import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from obspy.imaging.beachball import beach

MT = [1,0,-1,0,0,0]
b = beach(MT,xy=(slon,slat))

df = pd.read_csv('src_coordinates',names=['slat','slon','sdep'])
slat,slon,ddep = df['slat'][0],df['slon'][0],df['sdep'][0]

df = pd.read_csv('rcv_coordinates',names=['rcvnm','slat','slon'])
rcv_list = []
for i in range(len(df)):
    rcv_list.append([df['rcvnm'][i],df['slat'][i],df['slon'][i]])


fig, axes = pl
