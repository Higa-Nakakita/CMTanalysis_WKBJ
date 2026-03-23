import pandas as pd
import matplotlib.pyplot as plt
import sys

def visu_srcterm():
    filnm = '../ck_zAexpsrc.txt'
    names = ['n','j','frqprt','STmode','rcv','Re','Im']   
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
    st_modes = df_filtered['STmode'].unique()

    fig, axes = plt.subplots(2,2,figsize=(10,7))
    
    STnum = 0
    for idx, st_mode in enumerate(st_modes):
        df_mode = df_filtered[df_filtered['STmode'] == st_mode]
        n_values = df_mode['n'].unique()

        for n in n_values:
            df_n = df_mode[df_mode['n'] == n]
            axes[0,STnum].plot(df_n['frqprt'],df_n['Re'],label=f'Re (n={n})',marker='o',ms=1.5, linewidth=0)
            axes[1,STnum].plot(df_n['frqprt'],df_n['Im'],label=f'Im (n={n})',marker='o',ms=1.5, linewidth=0)

        axes[0, STnum].set_title(f'mode={st_mode}')
        axes[0, STnum].set_xlabel('Frequency (frqprt)')
        axes[0, STnum].set_ylabel('Re zAexpsrc')
        axes[0, STnum].legend(fontsize='small',loc='upper left', bbox_to_anchor=(1, 1))
        axes[1, STnum].set_title(f'mode={st_mode}')
        axes[1, STnum].set_xlabel('Frequency (frqprt)')
        axes[1, STnum].set_ylabel('Im zAexpsrc')
        axes[1, STnum].legend(fontsize='small',loc='upper left', bbox_to_anchor=(1, 1))
        STnum = STnum + 1

    # Display the plots
    plt.show()

# Call the function to execute
visu_srcterm()
