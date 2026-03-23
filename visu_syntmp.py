import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

def visu_syntemp():
    filnm = 'ck_syntmp.txt'
    #names = ['nt','STmode','comp','rcv','Re','Im','syn']
    names = ['nt','STmode','comp','rcv','Re','syn']
    df = pd.read_csv(filnm,header=None,names=names,delim_whitespace=True)

    list_nt = df['nt'].tolist()
    list_ST = df['STmode'].tolist()
    #list_Re = df['Re'].tolist()
    #list_Im = df['Im'].tolist()
    list_syn= df['syn'].tolist()

    # Get the rcvnum from user input
    rcvnum = int(input('Enter rcvnum: '))

    # Filter the dataframe by the selected rcvnum
    df_filtered = df[df['rcv'] == rcvnum]

    # Get the unique STmode values
    st_comps = df_filtered['comp'].unique()

    fig, axes = plt.subplots(3,3,figsize=(10,6))

    compnum = 0
    for idx, st_comp in enumerate(st_comps):
        df_comp = df_filtered[df_filtered['comp'] == st_comp]
        print(df_comp['nt'])
        #print(df_comp['Re'])
        print(df_comp['syn'])
        axes[0,compnum].plot(df_comp['nt'],df_comp['Re'],label=f'Re',marker='o',ms=1.5, linewidth=0.5)
        axes[1,compnum].plot(df_comp['nt'],df_comp['syn'],label=f'Im',marker='o',ms=1.5, linewidth=0.5)
        axes[2,compnum].plot(df_comp['nt'],df_comp['syn'],label=f'syn',marker='o',ms=1.5, linewidth=0.5)

        axes[0, compnum].set_title(f'comp={st_comp}')
        axes[0, compnum].set_xlabel('Frequency (frqprt)')
        axes[0, compnum].set_ylabel('Real(zwavefrq_conjg)')
        axes[0, compnum].legend(fontsize='small')
        axes[1, compnum].set_title(f'comp={st_comp}')
        axes[1, compnum].set_xlabel('Frequency (frqprt)')
        axes[1, compnum].set_xlim([0,600])
        axes[1, compnum].set_ylabel('Imag(zwavefrq_conjg)')
        axes[1, compnum].legend(fontsize='small')
        axes[2,compnum].set_ylabel('syntmp')
        axes[2,compnum].legend(fontsize='small')
        compnum = compnum + 1


    # Display the plots
    plt.show()


visu_syntemp()
