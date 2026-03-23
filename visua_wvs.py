import matplotlib.pyplot as plt
import pandas as pd
from obspy import read
from scipy.fft import fft, fftfreq
from obspy.signal.filter import bandpass
import glob
import os
import sys
from matplotlib.backends.backend_pdf import PdfPages

def visu_waves():

    OptorFor = 0 # OptorFor = 0: Optimal, 1:Forward
    OptorFor = 1
    if OptorFor == 0 :
        pdf_fig = PdfPages(f'waves_optimal.pdf')
    if OptorFor == 1 :
        pdf_fig = PdfPages(f'waves_forward.pdf')

    # obtain stanm
    rcv_coordi = 'rcv_coordinates'
    names = ['sta','lat','lon']
    df_rcv = pd.read_csv(rcv_coordi,names=names,delim_whitespace=True)
    sta_list = df_rcv['sta']
    stanm = len(sta_list)

    numtw = int(input('Enter num of tw:'))
    obs_dirpath = './Waveforms/'
    syn_dirpath = './Waveforms/'
    comp_list = ['LHZ','LHR','LHT']
    plots_perPage = 2
    #plots_perPage = 6
    pages = stanm // plots_perPage
    All_pages = pages * numtw

    print(All_pages,pages)
    cmpnm = len(comp_list)
    for tw in range(numtw):
        twnum = tw + 1
        for page in range(pages):
            fig, axes = plt.subplots(plots_perPage,cmpnm,figsize=(cmpnm*4,plots_perPage*2))
            for row in range(plots_perPage):
                sta_idx = page * plots_perPage + row
                sta = sta_list[sta_idx]
                print(sta_idx)
                for col in range(cmpnm):
                    comp = comp_list[col]
                    print(comp)

                    obs_path = obs_dirpath + f'obs_{sta}_tw{twnum}.{comp}.sac' 
                    if OptorFor == 0:
                        syn_path = syn_dirpath + f'*.syn_{sta}_tw{twnum}.{comp}.sac'
                    elif OptorFor == 1:
                        syn_path = syn_dirpath + f'syn_{sta}_tw{twnum}.{comp}.sac'
                    syn_path = glob.glob(syn_path)
                    print(syn_path)
                    syn_path = syn_path[0]

                    obs_TorF = os.path.exists(obs_path)
                    syn_TorF = os.path.exists(syn_path)

                    if (obs_TorF == False or syn_TorF== False) :
                         print('No file:',sta)
                         continue
                     
                    obs_st = read(obs_path)
                    syn_st = read(syn_path)

                    obs_time = obs_st[0].times()
                    obs_data = obs_st[0].data
                    syn_time = syn_st[0].times()
                    syn_data = syn_st[0].data

                    axes[row][col].plot(obs_time,obs_data,color='k')
                    axes[row][col].plot(syn_time,syn_data,color='r')
                    axes[row][col].set_xlabel('Time (s)')
                    axes[row][col].set_ylabel('(m)')
                    axes[row][col].set_title(f'{sta} {comp} tw{twnum}')
                    axes[row][col].ticklabel_format(style="sci",axis="y",scilimits=(0,0))
            plt.tight_layout()
            pdf_fig.savefig(fig)
    pdf_fig.close()

visu_waves()

