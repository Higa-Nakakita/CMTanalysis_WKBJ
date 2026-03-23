import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from obspy.imaging.beachball import beach

# ---------------------------------------------------------
# 1. データ読み込み (ユーザーのファイル指定に従う)
# ---------------------------------------------------------

try:
    # 震源 (Source) 座標の読み込み
    df_src = pd.read_csv('src_coordinates', names=['slat', 'slon', 'sdep'],delim_whitespace=True)
    # 震源緯度、経度を取得
    slat, slon = df_src['slat'][0], df_src['slon'][0]
    sdep = df_src['sdep'][0] # 深度はここではプロットに使わないが取得
    print(slat)
    print(slon)

    # 観測点 (Receiver) 座標の読み込み
    df_rcv = pd.read_csv('rcv_coordinates', names=['rcvnm', 'slat', 'slon'],delim_whitespace=True)
    # 観測点リストをデータフレームから直接利用
    
except FileNotFoundError as e:
    print(f"エラー: {e}")
    print("指定されたファイル名 'src_coordinates' または 'rcv_coordinates' が見つかりません。ファイル名を確認してください。")
    # 処理を中断
    exit()
except Exception as e:
    print(f"データの読み込み中にエラーが発生しました: {e}")
    exit()

# ---------------------------------------------------------
# 2. 震源メカニズム解の設定
# ---------------------------------------------------------

# モーメントテンソル (MT): [Mrr, Mtt, Mpp, Mrt, Mrp, Mtp]
# ユーザーコードの例を採用: [1, 0, -1, 0, 0, 0] (Double Coupleではないが、例として使用)
MT = [1, 0, -1, 0, 0, 0]

# ---------------------------------------------------------
# 3. Cartopyによるプロット
# ---------------------------------------------------------

# FigureとAxesの準備 (PlateCarree投影を使用)
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

# 地図要素の追加
ax.coastlines(resolution='50m', color='black', linewidth=0.5)
ax.add_feature(cfeature.LAND, facecolor='lightgray', zorder=0)
ax.add_feature(cfeature.OCEAN, facecolor='azure', zorder=0)
ax.add_feature(cfeature.BORDERS, linestyle=':', alpha=0.5)
ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False, color='gray', alpha=0.5)

# --- A. 観測点のプロット ---
ax.scatter(df_rcv['slon'], df_rcv['slat'], 
           color='blue', marker='^', s=100, label='観測点',
           transform=ccrs.PlateCarree(), zorder=5)

# 観測点名のラベル表示
for index, row in df_rcv.iterrows():
    ax.text(row['slon'] + 0.1, row['slat'], row['rcvnm'],
            transform=ccrs.PlateCarree(), fontsize=15, fontweight='bold', zorder=6)

# --- B. Beachballのプロット ---
# width: ビーチボールのサイズ (度単位)。地図のスケールに合わせて調整
beach_width = 2.0 
print(slon,slat)
b = beach(MT, xy=(slon, slat), width=beach_width, linewidth=1, facecolor='red', zorder=10)

# CartopyのAxesにビーチボール(Collection)を追加
ax.add_collection(b)

# --- C. 表示範囲の調整 ---
# 全データが含まれるように範囲計算（観測点、震源、ビーチボールを含む）
all_lons = list(df_rcv['slon']) + [slon]
all_lats = list(df_rcv['slat']) + [slat]
## ビーチボールのサイズも考慮
#margin = beach_width + 1.0 

buf = 1
ax.set_extent([min(all_lons)-buf, max(all_lons)+buf, 
               min(all_lats)-buf, max(all_lats)+buf], crs=ccrs.PlateCarree())

#plt.legend(loc='upper left')
plt.show()
