import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys


def rm_based_list(rm_list):

    rm_stanmlist = rm_list
    rcvcoordi = 'rcv_coordinates'
    obsdatapa = 'obsdata_path'
    names_rcvco = ['name','lat','lon']
    df_rcvco = pd.read_csv(rcvcoordi,names=names_rcvco,delim_whitespace=True)
    names_datapa = ['datapath']
    df_datapa = pd.read_csv(obsdatapa,names=names_datapa,delim_whitespace=True)

    extract_rcvs = []
    for k in range(len(df_rcvco['name'])):
        stanm = df_rcvco['name'][k]
        if stanm in rm_stanmlist:
            print('Remove ',stanm)
            continue
        else :
            extract_nm,extract_la,extract_lo = df_rcvco['name'][k],df_rcvco['lat'][k],df_rcvco['lon'][k]
            extract_rcvs.append([extract_nm,extract_la,extract_lo])

    extract_datas = []
    for k in range(len(df_datapa['datapath'])):
        datapath = df_datapa['datapath'][k]
        print(datapath)
        if any(rmsta+'.LH' in datapath for rmsta in rm_stanmlist):
            print('Remove ',datapath)
        else:
            extract_datapath = datapath
            extract_datas.append(datapath)

    with open(rcvcoordi,'w') as f:
        for k in range(len(extract_rcvs)):
            name,lat,lon = extract_rcvs[k][0],extract_rcvs[k][1],extract_rcvs[k][2]
            f.write(f'{name} {lat} {lon}\n')
    with open(obsdatapa ,'w') as f:
        for k in range(len(extract_datas)):
            datapath = extract_datas[k]
            f.write(f'{datapath}\n')

def rm_based_rcvcoordi():
    rcvcoordi = 'rcv_coordinates'
    obsdatapa = 'obs_datas.txt'
    #obsdatapa = 'obsdata_path'
    new_obsdatapa = 'new_obsdata_path.txt'
    names_rcvco = ['name','lat','lon']
    df_rcvco = pd.read_csv(rcvcoordi,names=names_rcvco,delim_whitespace=True)
    names_datapa = ['datapath']
    df_datapa = pd.read_csv(obsdatapa,names=names_datapa,delim_whitespace=True)

    extract_datas = []
    save_stanmlist = df_rcvco['name']
    for k in range(len(df_datapa['datapath'])):
        datapath = df_datapa['datapath'][k]
        print(datapath)
        if any(save_sta+'_tw1_full.LH' in datapath for save_sta in save_stanmlist):
            extract_datapath = datapath
            extract_datas.append(datapath)
        else:
            print('Remove ',datapath)


    with open(new_obsdatapa ,'w') as f:
        for k in range(len(extract_datas)):
            datapath = extract_datas[k]
            f.write(f'{datapath}\n')

def rm_chat():
    # ファイルパス
    file_A = "rcv_coordinates"  # ステーションリスト（ADM 37.9 138.43 ...）
    file_B = "obsdata_path"  # SACファイルリスト（フルパス）
    
    # ステーション名リストを作成
    with open(file_A, "r") as fa:
        sta_names = set(line.strip().split()[0] for line in fa)
    
    # Bから該当ステーション名の行だけ抽出
    filtered_lines = []
    with open(file_B, "r") as fb:
        for line in fb:
            for sta in sta_names:
                if f"{sta}" in line:
                    filtered_lines.append(line)
                    break  # 一致したら次の行へ
    
    # 新たなファイルに書き出し
    with open("obsdata_path", "w") as fout:
        fout.writelines(filtered_lines)
    


#rm_stanmlist = ['TMC','TKD','INN','YSI','NRW','YZK','SRN','SBT','NOK','OOW','NOP','KNP','NKG','HJO','WTR','TOW','TAS','KZS','CHS','TMC','SIB','SGN','FUJ','KNY','GJM','NSK','TOW','KZK','YHZ','SBT','NRW','UMJ','KSN','ONS','TGA','YTY','HID','ZMM','HRO','KSN','KMM','YMZ','ORH','KSR','NAA','KMT']
#rm_based_list(rm_stanmlist)
#rm_based_rcvcoordi()
rm_chat()
