import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pandas as pd
from obspy.imaging.beachball import beach

#basenm = "20200912024409_38.72_142.27_42.93_6.1"
#basenm = "20190509234841_31.8_131.97_25.46_6.2"
basenm = '.'
basedir = f'./{basenm}/'

#20200912024409_38.72_142.27_42.93_6.1
GCMT_loca = [38.70, 142.18, 43.3]
Mrr,Mtt,Mpp,Mrt,Mrp,Mtp = 1.350, 0.021, -1.370, 0.543, 1.300, -0.372

#20190509234841_31.8_131.97_25.46_6.2
#GCMT_loca = [31.78,131.92,30.4]
#Mrr,Mtt,Mpp,Mrt,Mrp,Mtp = 2.29,-0.608,-1.69,1.16,1.78,-1.01

GCMT = [Mrr,Mtt,Mpp,Mrt,Mrp,Mtp]


supfrq_list = ["10-20","10-30","10-35","10-40"] #mHz

strucnm = 'premiso'
#strucnm = 'premani'
#strucnm = 'kazuiso'
#strucnm = 'kazuani'
#strucnm = 'jmakazuiso'
#strucnm = 'jmakazuani'

dir_path = f'./Waveforms_{strucnm}'
dir_path = basedir + dir_path
###########
df_rcv = pd.read_csv(basedir+'rcv_coordinates',names=['sta','lat','lon'],delim_whitespace=True)

stanms,lats,lons = [],[],[]
for i in range(len(df_rcv)):
    stanms.append(df_rcv['sta'][i])
    lats.append(df_rcv['lat'][i])
    lons.append(df_rcv['lon'][i])

# VRデータ読み込み
cols = ['k', 'icomp', 'iw', 'vr_tmp1', 'vr_tmp2','mis']
df_vr = pd.read_csv(dir_path+'/VRs_forward.txt', delim_whitespace=True, names=cols)
print(df_vr)

# インデックス調整
df_vr['k'] -= 1
df_vr['icomp'] -= 1
df_vr['iw'] -= 1

# サイズ取得
rcvnum = df_vr['k'].max() + 1
ncomp = 3
nwave = df_vr['iw'].max() + 1

# VR配列構築
vr_array = np.full((rcvnum, ncomp, nwave), np.nan)
mis_array = np.full((rcvnum, ncomp, nwave), np.nan)
for _, row in df_vr.iterrows():
    k = int(row['k'])
    icomp = int(row['icomp'])
    iw = int(row['iw'])
    vr_array[k, icomp, iw] = row['vr_tmp1'] #Actually, this VR is invalid
    mis_array[k, icomp, iw]= row['mis']

##############


components = ['Vertical', 'Radial', 'Transverse']  # icomp = 0,1,2 に対応
ncomp = len(components)

# VR のレンジ（共通に 0–100 % とする想定なら固定）
vmin, vmax = 0.0, 100.0
cmap = 'viridis'  # 好きなカラーマップに変えてOK

# 図全体のサイズを時間窓数に応じて調整
fig_width = 12
fig_height = 3 * nwave  # 1行あたり高さ3インチくらい
fig, axes = plt.subplots(
    nwave, ncomp,
    figsize=(fig_width, fig_height),
    subplot_kw={"projection": ccrs.PlateCarree()}
)

# axes を 2次元配列として扱いやすくする
# （nwave=1 や ncomp=1 の場合でも形を揃えるため）
if nwave == 1:
    axes = np.expand_dims(axes, axis=0)
if ncomp == 1:
    axes = np.expand_dims(axes, axis=1)

for iw in range(nwave):         # 行: 時間窓
    for icomp, comp in enumerate(components):  # 列: 成分
        ax = axes[iw, icomp]

        # この時間窓・この成分の VR（観測点ごと）
        VR_val = vr_array[:, icomp, iw]  # shape = (rcvnum,)

        # NaN を除いた平均値（タイトルに表示）
        valid_mask = ~np.isnan(VR_val)
        if np.any(valid_mask):
            VR_mn_value = np.nanmean(VR_val)
        else:
            VR_mn_value = np.nan

        # ベースマップ
        ax.set_extent([120, 150, 23, 47])   # 日本周辺の例。必要に応じて変更
        ax.add_feature(cfeature.COASTLINE)
        ax.add_feature(cfeature.BORDERS)
        ax.add_feature(cfeature.LAND, facecolor='lightgray')
        ax.add_feature(cfeature.OCEAN, facecolor='lightblue')

        # VR を色として scatter（NaN は落とす）
        lons = df_rcv['lon'].values[valid_mask]
        lats = df_rcv['lat'].values[valid_mask]
        vals = VR_val[valid_mask]

        sc = ax.scatter(
            lons, lats,
            c=vals,
            cmap=cmap,
            s=50,
            vmin=vmin, vmax=vmax,
            edgecolor='k',
            transform=ccrs.PlateCarree()
        )

        # タイトル：成分名 + 平均VR + 時間窓インデックス
        print(iw+1,supfrq_list[iw])
        supfrq = supfrq_list[iw]
        #ax.set_title(f"{comp} {supfrq}mHz VR={VR_mn_value:.0f}%", fontsize=12)
        ax.set_title(f"{comp} {supfrq}mHz", fontsize=14)


        b = beach(GCMT,xy=(GCMT_loca[1],GCMT_loca[0]),width=2, linewidth=0.8, facecolor='r', edgecolor='k')
        ax.add_collection(b)


# 共通カラーバーを右側に配置
#cbar = fig.colorbar(sc, ax=axes.ravel().tolist(), orientation='vertical', fraction=0.01, pad=0.02)
#cbar.set_label("VR")

# 共通カラーバーを右端の図の右横に配置
# 3列目の各行にカラーバーをつける
for row in range(axes.shape[0]):  # 各行をループ
    ax = axes[row, -1]  # 3列目（最後の列）の軸を取得
    cbar = fig.colorbar(sc, ax=ax, orientation='vertical', fraction=0.05, pad=0.02)
    #cbar.set_label("VR")

plt.tight_layout()

#fig.suptitle(f"VR maps for {strucnm}", fontsize=16)
#plt.tight_layout(rect=[0, 0, 0.9, 0.9])  # 右にカラーバー、上にタイトル

#savedir = f'./Figs'
savedir =f'./'
plt.savefig(f"{savedir}/VR_maps_{strucnm}_{basenm}.pdf", dpi=300)
#plt.show()

