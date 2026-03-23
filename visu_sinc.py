import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys

fig, axes = plt.subplots(2,2)
names = ['n','omg','Re','Im']

df = pd.read_csv('ck_sinc.txt',names=names,delim_whitespace=True)

n_list = df['n'].tolist()
omg_list = df['omg'].tolist()
Re_list = df['Re'].tolist()
Im_list = df['Im'].tolist()

modes = df['n'].unique()
print(modes)
#sys.exit()

for idx, mode in enumerate(modes):
    df_fil = df[df['n'] == mode]
    axes[0][0].plot(df_fil['omg']/(2*np.pi),df_fil['Re'],label=f'Re (n={mode})',marker='o',ms=1.5,c='k')
    axes[0][1].plot(df_fil['omg']/(2*np.pi),df_fil['Im'],label=f'Im (n={mode})',marker='o',ms=1.5,c='k')
    axes[1][0].plot(df_fil['omg']/(2*np.pi),df_fil['Re'],label=f'Re (n={mode})',marker='o',ms=1.5,c='k')
    axes[1][0].plot(df_fil['omg']/(2*np.pi),df_fil['Im'],label=f'Im (n={mode})',marker='o',ms=1.5,c='c')
plt.show()
