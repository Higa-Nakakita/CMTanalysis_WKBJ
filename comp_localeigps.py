import pandas as pd
import matplotlib.pyplot as plt


### Attension kk of this code is not wavenumber, but station idx.
filnm = '../ck_local_eigp.txt'
names = ['idx', 'STmode', 'n', 'l', 'k', 'w', 'Q', 'gv']
df = pd.read_csv(filnm, names=names, delim_whitespace=True)

df_eigp = df.astype({
    'idx': str,
    'STmode': int,
    'n': int,
    'l': int,
    'k': int,
    'w': float,
    'Q': float,
    'gv': float,
})


maxN = 1

# ユニークなkごとに分ける
unique_k = df_eigp['k'].unique()

fig, axes = plt.subplots(3, 4, figsize=(12, 8))
for k_val in unique_k:
    df_k = df_eigp[df_eigp['k'] == k_val]

    # ユニークなidxごとに分ける
    unique_idx = df_k['idx'].unique()
    for idx_val in unique_idx:
        df_idx = df_k[df_k['idx'] == idx_val]

        unique_ST = df_idx['STmode'].unique()
        for ST in unique_ST:
           df_ST = df_idx[df_idx['STmode']==ST]


           # ユニークなnごとに分ける
           unique_n = df_ST['n'].unique()
           for n_val in unique_n:
               if n_val > maxN:
                   break 
               df_n = df_ST[df_idx['n'] == n_val]


               if (ST==1):
                   axes[0,0].plot(df_n['l'],df_n['w'],ms=1,c='0.6',label=f'{k_val},  idx={idx_val}, n={n_val}')
                   axes[1,0].plot(df_n['l'],df_n['Q'],ms=1,c='0.6',label=f'{k_val},  idx={idx_val}, n={n_val}')
                   axes[2,0].plot(df_n['l'],df_n['gv'],ms=1,c='0.6',label=f'{k_val},  idx={idx_val}, n={n_val}')
                   axes[0,2].plot(df_n['w'],df_n['l'],ms=1,c='0.6',label=f'{k_val},  idx={idx_val}, n={n_val}')
                   axes[1,2].plot(df_n['w'],df_n['Q'],ms=1,c='0.6',label=f'{k_val},  idx={idx_val}, n={n_val}')
                   axes[2,2].plot(df_n['w'],df_n['gv'],ms=1,c='0.6',label=f'{k_val},  idx={idx_val}, n={n_val}')
               if (ST==2):
                   axes[0,1].plot(df_n['l'],df_n['w'],ms=1,c='0.7',label=f'{k_val},  idx={idx_val}, n={n_val}')
                   axes[1,1].plot(df_n['l'],df_n['Q'],ms=1,c='0.6',label=f'{k_val},  idx={idx_val}, n={n_val}')
                   axes[2,1].plot(df_n['l'],df_n['gv'],ms=1,c='0.6',label=f'{k_val},  idx={idx_val}, n={n_val}')
                   axes[0,3].plot(df_n['w'],df_n['l'],ms=1,c='0.6',label=f'{k_val},  idx={idx_val}, n={n_val}')
                   axes[1,3].plot(df_n['w'],df_n['Q'],ms=1,c='0.6',label=f'{k_val},  idx={idx_val}, n={n_val}')
                   axes[2,3].plot(df_n['w'],df_n['gv'],ms=1,c='0.6',label=f'{k_val},  idx={idx_val}, n={n_val}')



for k_val in unique_k:
    df_k = df_eigp[df_eigp['k'] == k_val]

    df_prem = df_k[df_k['idx'] == 'PREM']
    unique_ST = df_prem['STmode'].unique()

    for ST in unique_ST:
        df_ST = df_prem[df_prem['STmode']==ST]

        # ユニークなnごとに分ける
        unique_n = df_ST['n'].unique()
        for n_val in unique_n:
            if n_val > maxN:
                break 
            df_n = df_ST[df_prem['n'] == n_val]


            if (ST==1):
                axes[0,0].plot(df_n['l'],df_n['w'],ms=1,c='r',label=f'{k_val},  PREM, n={n_val}')
                axes[1,0].plot(df_n['l'],df_n['Q'],ms=1,c='r',label=f'{k_val},  idx={idx_val}, n={n_val}')
                axes[2,0].plot(df_n['l'],df_n['gv'],ms=1,c='r',label=f'{k_val},  idx={idx_val}, n={n_val}')
                axes[0,2].plot(df_n['w'],df_n['l'],ms=1,c='r',label=f'{k_val},  idx={idx_val}, n={n_val}')
                axes[1,2].plot(df_n['w'],df_n['Q'],ms=1,c='r',label=f'{k_val},  idx={idx_val}, n={n_val}')
                axes[2,2].plot(df_n['w'],df_n['gv'],ms=1,c='r',label=f'{k_val},  idx={idx_val}, n={n_val}')
            if (ST==2):
                axes[0,1].plot(df_n['l'],df_n['w'],ms=1,c='r',label=f'{k_val},  PREM, n={n_val}')
                axes[1,1].plot(df_n['l'],df_n['Q'],ms=1,c='r',label=f'{k_val},  idx={idx_val}, n={n_val}')
                axes[2,1].plot(df_n['l'],df_n['gv'],ms=1,c='r',label=f'{k_val},  idx={idx_val}, n={n_val}')
                axes[0,3].plot(df_n['w'],df_n['l'],ms=1,c='r',label=f'{k_val},  idx={idx_val}, n={n_val}')
                axes[1,3].plot(df_n['w'],df_n['Q'],ms=1,c='r',label=f'{k_val},  idx={idx_val}, n={n_val}')
                axes[2,3].plot(df_n['w'],df_n['gv'],ms=1,c='r',label=f'{k_val},  idx={idx_val}, n={n_val}')

axes[0,2].set_xlim([0.005,0.4])
axes[1,2].set_xlim([0.005,0.4])
axes[2,2].set_xlim([0.005,0.4])
axes[0,3].set_xlim([0.005,0.4])
axes[1,3].set_xlim([0.005,0.4])
axes[2,3].set_xlim([0.005,0.4])

plt.show()
