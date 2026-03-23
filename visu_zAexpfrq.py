import pandas as pd
import matplotlib.pyplot as plt
import sys

def visu_zAexpfrq():
    filnm = 'ck_zAexpfrq.txt'
    names = ['n','j','frqprt','STmode','icomp','rcv','Re','Im']   
    df = pd.read_csv(filnm, header=None, names=names, delim_whitespace=True)

    list_n = df['n'].tolist()
    list_j = df['j'].tolist()
    list_frq=df['frqprt'].tolist()
    list_ST = df['STmode'].tolist()
    list_Re = df['Re'].tolist()
    list_Im = df['Im'].tolist()

    # Get the rcvnum from user input
    rcvnum = int(input('Enter rcvnum: '))

    # Filter the dataframe by the selected rcvnum
    df_filtered = df[df['rcv'] == rcvnum]

    # Get the unique STmode values
    st_comps = df_filtered['icomp'].unique()

    fig, axes = plt.subplots(2,3,figsize=(10,6))
    
    compnum = 0
    for idx, st_comp in enumerate(st_comps):
        df_comp = df_filtered[df_filtered['icomp'] == st_comp]
        n_values = df_comp['n'].unique()

        for n in n_values:
            if n==0 or n==1 :
                continue
            df_n = df_comp[df_comp['n'] == n]
            axes[0,compnum].plot(df_n['frqprt'],df_n['Re'],label=f'Re (n={n})',marker='o',ms=1.5, linewidth=0.5)
            axes[1,compnum].plot(df_n['frqprt'],df_n['Im'],label=f'Im (n={n})',marker='o',ms=1.5, linewidth=0.5)

        axes[0, compnum].set_title(f'comp={st_comp}')
        axes[0, compnum].set_xlabel('Frequency (frqprt)')
        axes[0, compnum].set_ylabel('Real(zAexpfrq)')
        axes[0, compnum].legend(fontsize='small')
        axes[1, compnum].set_title(f'comp={st_comp}')
        axes[1, compnum].set_xlabel('Frequency (frqprt)')
        axes[1, compnum].set_ylabel('Imag(zAexpfrq)')
        axes[1, compnum].legend(fontsize='small')
        compnum = compnum + 1

    # Display the plots
    plt.show()

def visu_zwavefrq():
    filnm = 'ck_zwavefrq.txt'
    names = ['j','frqprt','STmode','icomp','rcv','Re','Im']
    df = pd.read_csv(filnm, header=None, names=names, delim_whitespace=True)

    list_j = df['j'].tolist()
    list_frq=df['frqprt'].tolist()
    list_ST = df['STmode'].tolist()
    list_Re = df['Re'].tolist()
    list_Im = df['Im'].tolist()

    # Get the rcvnum from user input
    rcvnum = int(input('Enter rcvnum: '))

    # Filter the dataframe by the selected rcvnum
    df_filtered = df[df['rcv'] == rcvnum]

    # Get the unique STmode values
    st_comps = df_filtered['icomp'].unique()

    fig, axes = plt.subplots(2,3,figsize=(10,6))

    compnum = 0
    for idx, st_comp in enumerate(st_comps):
        df_comp = df_filtered[df_filtered['icomp'] == st_comp]

        axes[0,compnum].plot(df_comp['frqprt'],df_comp['Re'],label=f'Re',marker='o',ms=1.5, linewidth=0.5)
        axes[1,compnum].plot(df_comp['frqprt'],df_comp['Im'],label=f'Im',marker='o',ms=1.5, linewidth=0.5)

        axes[0, compnum].set_title(f'comp={st_comp}')
        axes[0, compnum].set_xlabel('Frequency (frqprt)')
        axes[0, compnum].set_ylabel('Real(zwavefrq)')
        axes[0, compnum].legend(fontsize='small')
        axes[1, compnum].set_title(f'comp={st_comp}')
        axes[1, compnum].set_xlabel('Frequency (frqprt)')
        axes[1, compnum].set_ylabel('Imag(zwavefrq)')
        axes[1, compnum].legend(fontsize='small')
        compnum = compnum + 1

    # Display the plots
    plt.show()

def visu_zwavefrq_conjg():
    filnm = 'ck_zwavefrq_conjg.txt'
    names = ['nt','STmode','icomp','rcv','Re','Im']
    df = pd.read_csv(filnm, header=None, names=names, delim_whitespace=True)

    list_nt = df['nt'].tolist()
    list_ST = df['STmode'].tolist()
    list_Re = df['Re'].tolist()
    list_Im = df['Im'].tolist()

    # Get the rcvnum from user input
    rcvnum = int(input('Enter rcvnum: '))

    # Filter the dataframe by the selected rcvnum
    df_filtered = df[df['rcv'] == rcvnum]

    # Get the unique STmode values
    st_comps = df_filtered['icomp'].unique()

    fig, axes = plt.subplots(2,3,figsize=(10,6))

    compnum = 0
    for idx, st_comp in enumerate(st_comps):
        df_comp = df_filtered[df_filtered['icomp'] == st_comp]

        axes[0,compnum].plot(df_comp['nt'],df_comp['Re'],label=f'Re',marker='o',ms=1.5, linewidth=0.5)
        axes[1,compnum].plot(df_comp['nt'],df_comp['Im'],label=f'Im',marker='o',ms=1.5, linewidth=0.5)

        axes[0, compnum].set_title(f'comp={st_comp}')
        #axes[0, compnum].set_xlabel('Frequency (frqprt)')
        axes[0, compnum].set_ylabel('Real(zwavefrq_conjg)')
        axes[0, compnum].legend(fontsize='small')
        axes[1, compnum].set_title(f'comp={st_comp}')
        #axes[1, compnum].set_xlabel('Frequency (frqprt)')
        axes[1, compnum].set_ylabel('Imag(zwavefrq_conjg)')
        axes[1, compnum].legend(fontsize='small')
        compnum = compnum + 1

    # Display the plots
    plt.show()

# Call the function to execute
visu_zAexpfrq()
#visu_zwavefrq()
#visu_zwavefrq_conjg()
