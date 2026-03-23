import pandas as pd
import matplotlib.pyplot as plt
import sys

def visu_fixterm():
    # Read the data from the file
    filnm = '../ck_zAexpfix.txt'
    names = ['n', 'j', 'frqprt', 'STmode','comp','rcv', 'Re', 'Im']
    df = pd.read_csv(filnm, header=None, names=names, delim_whitespace=True)

    # Get lists from the dataframe
    list_n = df['n'].tolist()
    list_j = df['j'].tolist()
    list_frq = df['frqprt'].tolist()
    list_ST = df['STmode'].tolist()
    list_comp= df['comp'].tolist()
    list_Re = df['Re'].tolist()
    list_Im = df['Im'].tolist()

    # Get the rcvnum from user input
    rcvnum = int(input('Enter rcvnum: '))
    print(f'You entered rcvnum: {rcvnum}')

    # Filter the dataframe by the selected rcvnum
    df_filtered = df[df['rcv'] == rcvnum]

    # Get the unique STmode values
    st_modes = df_filtered['STmode'].unique()
    st_comps = df_filtered['comp'].unique()

    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(10, 7))
    fig.tight_layout(pad=4.0)
    
    compnum = 0
    for idx, st_comp in enumerate(st_comps):
        # Filter by STmode
        df_comp = df_filtered[df_filtered['comp'] == st_comp]

        # Get unique n values
        n_values = df_comp['n'].unique()

        # Plot Re and Im on separate axes
        for n in n_values:
            # Filter by n
            df_n = df_comp[df_comp['n'] == n]

            # Plot Re and Im
            axes[0, compnum].plot(df_n['frqprt'], df_n['Re'], label=f'Re (n={n})', marker='o',ms=1.5,linewidth=0.5)
            axes[1, compnum].plot(df_n['frqprt'], df_n['Im'], label=f'Im (n={n})', marker='o',ms=1.5,linewidth=0.5)
        
        # Set titles and labels
        axes[0, compnum].set_title(f'comp={st_comp}')
        axes[0, compnum].set_xlabel('Frequency (frqprt)')
        axes[0, compnum].set_ylabel('Re zAexpfix')
        axes[0, compnum].legend(fontsize='small',loc='upper left', bbox_to_anchor=(1, 1))
        axes[1, compnum].set_title(f'comp={st_comp}')
        axes[1, compnum].set_xlabel('Frequency (frqprt)')
        axes[1, compnum].set_ylabel('Im zAexpfix')
        axes[1, compnum].legend(fontsize='small',loc='upper left', bbox_to_anchor=(1, 1))
        compnum = compnum + 1

    # Display the plots
    plt.show()

# Call the function to execute
visu_fixterm()

