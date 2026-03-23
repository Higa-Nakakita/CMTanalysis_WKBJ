import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

rate = 20

rcvfil = '../rcv_coordinates'
names = ['sta','lat','lon']
df_rcv = pd.read_csv(rcvfil,names=names,delim_whitespace=True)
sta_list = df_rcv['sta']
print(sta_list)

filnm = '../VRs_info.txt' 
names = ['ita','jj','sta','comp','iw','VRs','VR']
df = pd.read_csv(filnm,names=names,delim_whitespace=True)

plots_perpage = 3
numtw = 3
unique_sta = df['sta'].unique()
stanm = len(unique_sta)
pages = stanm // plots_perpage

comps = ['LHZ','LHR','LHT']
cols = ['r','b','g']
tws = ['tw1','tw2','tw3']


VRs_record = []
for page in range(pages):
    fig, axes = plt.subplots(plots_perpage,numtw,figsize=(numtw*4,plots_perpage*2))
    for row in range(plots_perpage):
        sta_idx = page*plots_perpage + row
        #print('staidx',sta_idx)
        sta = sta_list[sta_idx]

        df_sta = df[df['sta']==sta_idx+1]

        for col in range(numtw):
            tw = tws[col]
            df_iw = df_sta[df_sta['iw']==col+1]
            unique_comp = df_iw['comp'].unique()
            for cmp_idx in  unique_comp:
                df_comp = df_iw[df_iw['comp']==cmp_idx]
                #print(row,col,cmp_idx)
                #print(df_comp)
                comp = comps[cmp_idx-1]
                axes[row][col].plot(df_comp['jj'],df_comp['VRs'],label=f'{comp}',marker='o',linewidth=0,ms=0.5,c=cols[cmp_idx-1])

                N = int((1-rate*1e-2)*len(df_comp['VRs']))
                mean = df_comp['VRs'][N:].mean()
                std = df_comp['VRs'][N:].std()
                #print(sta,col,cmp_idx,mean,std)
                VRs_record.append([sta,col+1,cmp_idx,mean,std])


            axes[row][col].set_xlabel('The number of models') 
            axes[row][col].set_ylabel('Variance Reduction (%)') 
            axes[row][col].set_ylim(0,100)
            axes[row][col].set_title(f'{sta} {tw}')
        plt.tight_layout()
        #plt.legend()


df_vr = pd.DataFrame(VRs_record, columns=['sta', 'tw', 'comp', 'mean', 'std'])

# tw と comp ごとに平均を算出
mean_by_tw_comp = df_vr.groupby(['tw', 'comp'])['mean'].mean()
print(mean_by_tw_comp)



with open('vrs_record.txt','w') as f:
    for i in range(len(VRs_record)):
        sta,tw,cmp_idx,mean,std = VRs_record[i]
        f.write(f'{sta} {tw} {cmp_idx} {mean} {std}\n')
    f.write(f'{mean_by_tw_comp}')

plt.show()



