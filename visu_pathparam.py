import pandas as pd
import matplotlib.pyplot as plt
import sys
import numpy as np

pi2 = np.pi * 2

def visu_pathparam():
    filnm = 'ck_pathparam.txt'
    names = ['n','jjtmp','omgtmp','STmode','rcv','Qp','gv','atten','geosp','phstmp','zphs']
    df = pd.read_csv(filnm, header=None, names=names, delim_whitespace=True)

    list_n =  df['n'].tolist()
    list_jj=  df['jjtmp'].tolist()
    list_omg= df['omgtmp'].tolist()
    list_ST = df['STmode'].tolist()
    list_Qp = df['Qp'].tolist()
    list_gv = df['gv'].tolist()
    list_atte=df['atten'].tolist()
    list_geosp=df['geosp'].tolist()
    list_phs = df['phstmp'].tolist()
    list_zphs =df['zphs'].tolist()

    rcvnum = int(input('Enter rcvnum:'))
    
    df_filtered = df[df['rcv'] == rcvnum]

    st_modes = df_filtered['STmode'].unique()
    fig, axes = plt.subplots(2,2,figsize=(10,7))

    STnum = 0
    for idx, st_mode in enumerate(st_modes):
        df_mode = df_filtered[df_filtered['STmode'] == st_mode]
        n_values = df_mode['n'].unique()

        for n in n_values:
            df_n = df_mode[df_mode['n'] == n]
            axes[0, STnum].plot((df_n['omgtmp'])/pi2, df_n['Qp'], label=f' Qp(n={n})', marker='o',ms=1.0,linewidth=0)
            axes[1, STnum].plot((df_n['omgtmp'])/pi2, df_n['gv'], label=f' gv(n={n})', marker='o',ms=1.0,linewidth=0)
            #axes[2, STnum].plot((df_n['omgtmp'])/pi2, df_n['atten'], label=f' atten(n={n})', marker='o',ms=1.0,linewidth=0)
            #axes[3, STnum].plot((df_n['omgtmp'])/pi2, df_n['geosp'], label=f' geosp(n={n})', marker='o',ms=1.0,linewidth=0)
            #axes[4, STnum].plot((df_n['omgtmp'])/pi2, df_n['phstmp'], label=f' phs(n={n})', marker='o',ms=1.0,linewidth=0)

        # Set titles and labels
        axes[0, STnum].set_title(f'mode={st_mode}')
        axes[0, STnum].set_xlabel('Frequency (frqprt)')
        axes[0, STnum].set_ylabel('Qp')
        axes[0, STnum].legend(fontsize='small')
        #axes[1, STnum].set_title(f'mode={st_mode}')
        axes[1, STnum].set_xlabel('Frequency (frqprt)')
        axes[1, STnum].set_ylabel('gv')
        axes[1, STnum].legend(fontsize='small')
        STnum = STnum + 1

    # Display the plots
    plt.show()

# Call the function to execute

visu_pathparam()
