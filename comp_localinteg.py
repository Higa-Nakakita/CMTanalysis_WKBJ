import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

d2r = np.pi/180 # degree -> radian
d2km = d2r * 6371 # degree -> km
rad2km = 6371 #radian -> km

filnm = '../ck_local_integ.txt'
names = ['sta','idx', 'del', 'STmode', 'n', 'w', 'k', 'Q', 'gv', 'k_prem', 'Q_prem', 'gv_prem']
df = pd.read_csv(filnm, names=names, delim_whitespace=True)

df_eigp = df.astype({
    'sta': str,
    'idx': float,
    'del': float,
    'STmode': int,
    'n': int,
    'w': float,
    'k': float,
    'Q': float,
    'gv': float,
    'k_prem': float,
    'Q_prem': float,
    'gv_prem': float,
})

#omg = float(input("input frq(mHz):"))
omg = 10
omg = 2*np.pi*omg*1e-3
print(omg)

#nmode = 0
nmode = 1


sta_idx = input('input station idx: ')
df_sta = df_eigp[df_eigp['sta'] == sta_idx]

min_idx = (df_sta['w']-omg).abs().idxmin()
min_dif_w = df_sta['w'].loc[min_idx]
df_w = df_sta[df_sta['w']==min_dif_w]
#print(df_w)


fig, axes = plt.subplots(2, 2, figsize=(12, 8))

unique_ST = df_w['STmode'].unique()
for ST in unique_ST:
    df_ST = df_w[df_w['STmode']==ST]
    #print(df_ST['idx'])

    unique_n = df_ST['n'].unique()
    for n_val in unique_n:
        #if n_val != nmode:
        if n_val > nmode:
            break 
        df_n = df_ST[df_w['n'] == n_val]
        #print(df_n)

        if (ST==1):
            list_del = list(df_n['del'])
            list_gv = list(df_n['gv'])
            #print(list_del)
            #print(list_gv)
            axes[0,0].plot(df_n['del']*rad2km,df_n['k'],ms=3,c='0.6',marker='o',linewidth=0,label=f'k at n={n_val} omg={min_dif_w}')
            axes[0,0].plot(df_n['del']*rad2km,df_n['k_prem'],ms=1,c='r',label=f'k at n={n_val} omg={min_dif_w}')
            #axes[1,0].plot(df_n['del'],df_n['gv'],ms=5,c='0.6',marker='o',label=f'gv at n={n_val} omg={min_dif_w}')
            axes[1,0].plot(df_n['del']*rad2km,list_gv,ms=3,c='0.6',marker='o',linewidth=0,label=f'gv at n={n_val} omg={min_dif_w}')
            axes[1,0].plot(df_n['del']*rad2km,df_n['gv_prem'],ms=1,c='r',label=f'gv at n={n_val} omg={min_dif_w}')
            #print(df_n['idx'])
            #print(df_n['gv'])
            #print(df_n['k'])
            #print(df_n['idx'][0:len(df_n['idx'])])
            #sys.exit()
        if (ST==2):
            axes[0,1].plot(df_n['del']*rad2km,df_n['k'],ms=1,c='0.6',label=f'k at n={n_val} omg={min_dif_w}')
            axes[0,1].plot(df_n['del']*rad2km,df_n['k_prem'],ms=1,c='r',label=f'k at n={n_val} omg={min_dif_w}')
            axes[1,1].plot(df_n['del']*rad2km,df_n['gv'],ms=1,c='0.6',label=f'gv at n={n_val} omg={min_dif_w}')
            axes[1,1].plot(df_n['del']*rad2km,df_n['gv_prem'],ms=1,c='r',label=f'gv at n={n_val} omg={min_dif_w}')


plt.tight_layout()
plt.legend()
plt.show()
