import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

def visu_eif_nl(Flag):
    if Flag ==0: #rcv
        print('rcv')
        names = ['n', 'l', 'U', 'V', 'W']
        filnm = 'ck_eifrcv_nl'
        nm = 'rcv'
    elif Flag ==1: #src
        print('src')
        names = ['n', 'l', 'U', 'dU','V','dV','W','dW']
        filnm = 'ck_eifsrc_nl'
        nm = 'src'
    else:
        print("Wrong in Flag?")

    # ファイル読み込み
    df = pd.read_csv(f'../{filnm}.txt', names=names, delim_whitespace=True)
    n_values = sorted(df['n'].unique())

    # PDF出力の準備
    with PdfPages(f'{filnm}.pdf') as pdf:
        for page_start in range(0, len(n_values), 5):  # 5つのnごとに1ページ
            fig, axes = plt.subplots(3, 5, figsize=(9, 9))  # 5行×3列（U, V, W）×3モード
            fig.subplots_adjust(hspace=0.5)

            for i in range(5):
                if page_start + i >= len(n_values):
                    break
                n = n_values[page_start + i]
                df_n = df[df['n'] == n]

                # 各成分を描画
                axes[0, i].plot(df_n['l'], df_n['U'], marker='o', ms=1, linewidth=0)
                axes[0, i].set_title(f'n={n} (U)')
                axes[0, i].set_xlabel('l')
                axes[0, i].set_ylabel(f'U{nm}')

                axes[1, i].plot(df_n['l'], df_n['V'], marker='o', ms=1, linewidth=0)
                axes[1, i].set_title(f'n={n} (V)')
                axes[1, i].set_xlabel('l')
                axes[1, i].set_ylabel(f'V{nm}')

                axes[2, i].plot(df_n['l'], df_n['W'], marker='o', ms=1, linewidth=0)
                axes[2, i].set_title(f'n={n} (W)')
                axes[2, i].set_xlabel('l')
                axes[2, i].set_ylabel(f'W{nm}')

            # 空のサブプロットを非表示にする
            for row in range(3):
                for col in range(3):
                    if page_start + col >= len(n_values):
                        axes[row, col].axis('off')

            pdf.savefig(fig)
            plt.close(fig)


Flag = 0
Flag = 1
visu_eif_nl(Flag)
