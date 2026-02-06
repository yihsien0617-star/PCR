# @title 🧬 PCR 戰情室 (v9.0 動態熱模擬版)
# 加入「溫度滑桿」，可即時模擬不同 Annealing 溫度下的電泳結果

%matplotlib inline
import ipywidgets as widgets
from IPython.display import display, clear_output
import matplotlib.pyplot as plt
import numpy as np
import threading
import time

# 確保 Biopython 已載入
try:
    import Bio
    from Bio.Blast import NCBIWWW
    from Bio.Blast import NCBIXML
    from Bio.SeqUtils import MeltingTemp as mt
    from Bio.Seq import Seq
except ImportError:
    print("⚠️ 尚未偵測到 Biopython，請在 Anaconda Prompt 輸入 'pip install biopython' 安裝。")

# ================= 全域變數 (儲存暫存資料用) =================
raw_blast_data = {"f": [], "r": []} # 儲存原始 BLAST 結果
primer_info = {"tm_f": 0, "tm_r": 0, "seq_f": "", "seq_r": ""}

# ================= 核心計算與繪圖 =================

def calculate_adjusted_tm(perfect_tm, identity_percent):
    """
    估算不完美結合時的 Tm 值
    經驗法則：每 1% 的錯配 (Mismatch) 會導致 Tm 下降約 1-1.5°C
    """
    mismatch_percent = 100 - identity_percent
    drop = mismatch_percent * 1.2 # 設定係數
    return perfect_tm - drop

def draw_dynamic_gel(products, current_temp, organism):
    """繪製電泳圖 (互動版)"""
    # 清除舊圖 (避免重疊)
    plt.close('all')
    
    fig = plt.figure(figsize=(8, 5), dpi=100)
    ax = fig.add_subplot(111)
    ax.set_facecolor('black')
    
    # 繪製 Marker
    ladder = [100, 200, 300, 400, 500, 600, 800, 1000, 1500, 2000, 3000]
    for b in ladder:
        y = -np.log10(b)
        ax.hlines(y, 0.5, 1.5, colors='white', alpha=0.5)
        ax.text(0.1, y, f"{b}", color='white', fontsize=7, va='center')

    # 繪製通過溫度篩選的產物
    visible_count = 0
    for i, p in enumerate(products):
        size = p['size']
        y = -np.log10(size)
        
        # 判斷是否為主要產物 (通常 Identity 也是最高的)
        is_target = (p['f_ident'] == 100 and p['r_ident'] == 100)
        
        if is_target:
            color, alpha, lw = '#00FF00', 1.0, 3.5 # 亮綠色
            label = f"{size}bp (Target)"
        else:
            color, alpha, lw = 'yellow', 0.6, 2 # 黃色雜訊
            label = f"{size}bp (Non-specific)"
            
        ax.hlines(y, 2.5, 3.5, colors=color, alpha=alpha, linewidth=lw)
        # 為了避免文字重疊，只標示前幾個
        if visible_count < 5:
            ax.text(3.7, y, label, color=color, fontsize=9, va='center')
        visible_count += 1

    ax.set_xlim(0, 5)
    ax.set_ylim(-np.log10(3500), -np.log10(80))
    ax.set_xticks([1, 3])
    ax.set_xticklabels(['Marker', f'PCR @ {current_temp}°C'])
    ax.set_yticks([])
    ax.set_title(f"In-Silico PCR Simulation ({organism})", color='white')
    
    plt.tight_layout()
    display(fig) # 顯示圖表

# ================= 互動模擬邏輯 =================

def update_simulation(temp_val):
    """當滑桿移動時，重新過濾並畫圖"""
    
    with sim_output:
        clear_output(wait=True)
        
        # 取得儲存的原始資料
        f_hits = raw_blast_data['f']
        r_hits = raw_blast_data['r']
        tm_perfect_f = primer_info['tm_f']
        tm_perfect_r = primer_info['tm_r']
        
        products = []
        
        # 交叉比對並加入溫度過濾
        for f in f_hits:
            # 計算這條引子在這個位置的實際 Tm
            hit_tm_f = calculate_adjusted_tm(tm_perfect_f, f['ident'])
            
            # 如果結合力低於設定溫度，引子會脫落 (不發生反應)
            if hit_tm_f < temp_val: continue

            for r in r_hits:
                # 判斷是否在同一條染色體且方向正確
                if f['id'] == r['id'] and f['str'] != r['str']:
                    
                    # 計算 Reverse 引子的實際 Tm
                    hit_tm_r = calculate_adjusted_tm(tm_perfect_r, r['ident'])
                    if hit_tm_r < temp_val: continue
                    
                    # 計算產物大小
                    dist = abs(r['s'] - f['s']) + 1
                    if 50 < dist < 5000:
                        # 避免重複
                        if not any(p['size'] == dist for p in products):
                            products.append({
                                'size': dist, 
                                'gene': f['t'],
                                'f_ident': f['ident'],
                                'r_ident': r['ident']
                            })
        
        products.sort(key=lambda x: x['size'])
        
        # 顯示文字報告
        print(f"🌡️ 目前設定溫度: {temp_val}°C")
        print(f"--------------------------------")
        if not products:
            print("❄️ 溫度過高！引子無法結合，沒有 PCR 產物。")
            print("   (試著調低溫度滑桿)")
        else:
            print(f"📊 模擬結果: 共發現 {len(products)} 條產物")
            for p in products:
                status = "✅ 完美結合" if (p['f_ident']==100 and p['r_ident']==100) else "⚠️ 非專一性"
                print(f"   - {p['size']} bp [{status}] (F:{p['f_ident']}%/R:{p['r_ident']}%)")
            
            # 畫圖
            draw_dynamic_gel(products, temp_val, w_org.value)

# ================= 介面宣告 =================

style = {'description_width': 'initial'}
layout = widgets.Layout(width='98%')

# 輸入區
w_org = widgets.Dropdown(options=[('人類', 'Homo sapiens'), ('小鼠', 'Mus musculus')], value='Homo sapiens', description='物種:', style=style)
w_db = widgets.Dropdown(options=[('mRNA', 'refseq_rna'), ('Genomic DNA', 'nt')], value='refseq_rna', description='資料庫:', style=style)
w_fwd = widgets.Text(placeholder='Forward Primer', description='Fwd:', style=style, layout=layout)
w_rev = widgets.Text(placeholder='Reverse Primer', description='Rev:', style=style, layout=layout)

# 控制區
btn_run = widgets.Button(description="1. 取得引子資料 (連線 NCBI)", button_style='primary', icon='cloud-download', layout=widgets.Layout(width='100%'))
w_progress = widgets.IntProgress(value=0, max=100, description='進度:', bar_style='info', layout=widgets.Layout(width='100%'))
log_output = widgets.Output() # 用來顯示連線進度

# 模擬區 (一開始隱藏，跑完才出來)
w_slider = widgets.FloatSlider(value=55.0, min=40.0, max=75.0, step=0.5, description='Annealing Temp (°C):', style=style, layout=layout)
sim_output = widgets.Output() # 用來顯示即時圖表

# ================= 執行邏輯 =================

def run_blast_thread(f_seq, r_seq, organism, database):
    global raw_blast_data, primer_info
    
    with log_output:
        try:
            w_progress.value = 10
            print(f"🚀 開始連線 NCBI... 目標: {organism}")
            
            # 計算完美 Tm
            tm_f = mt.Tm_NN(Seq(f_seq.upper()))
            tm_r = mt.Tm_NN(Seq(r_seq.upper()))
            primer_info = {"tm_f": tm_f, "tm_r": tm_r, "seq_f": f_seq, "seq_r": r_seq}
            print(f"   引子 Tm: F={tm_f:.1f}°C, R={tm_r:.1f}°C")

            # BLAST (使用較寬鬆的 expect 值以抓取非專一性 hits)
            w_progress.value = 30
            print("   ⏳ 正在下載 Forward 引子結合位點...")
            rh_f = NCBIWWW.qblast("blastn", database, f_seq, entrez_query=f'"{organism}"[Organism]', expect=50000, word_size=7)
            rec_f = NCBIXML.read(rh_f)
            
            w_progress.value = 60
            print("   ⏳ 正在下載 Reverse 引子結合位點...")
            rh_r = NCBIWWW.qblast("blastn", database, r_seq, entrez_query=f'"{organism}"[Organism]', expect=50000, word_size=7)
            rec_r = NCBIXML.read(rh_r)
            
            w_progress.value = 90
            print("   ✅ 資料下載完成！正在啟動熱模擬引擎...")
            
            # 解析並儲存 Raw Data
            raw_blast_data['f'] = []
            raw_blast_data['r'] = []
            
            # 提取詳細資訊 (包含 identity)
            for a in rec_f.alignments:
                for h in a.hsps:
                    ident = (h.identities / h.align_length) * 100
                    raw_blast_data['f'].append({'id': a.accession, 's': h.sbjct_start, 'str': h.frame[1], 't': a.title, 'ident': ident})

            for a in rec_r.alignments:
                for h in a.hsps:
                    ident = (h.identities / h.align_length) * 100
                    raw_blast_data['r'].append({'id': a.accession, 's': h.sbjct_start, 'str': h.frame[1], 't': a.title, 'ident': ident})

            w_progress.value = 100
            w_progress.bar_style = 'success'
            print("🎉 準備就緒！請使用下方滑桿調整溫度。")
            
            # 啟動模擬介面
            w_slider.value = min(tm_f, tm_r) - 5 # 預設溫度
            widgets.interactive(update_simulation, temp_val=w_slider)
            update_simulation(w_slider.value) # 觸發第一次繪圖

        except Exception as e:
            print(f"❌ 錯誤: {e}")
        finally:
            btn_run.disabled = False

def on_button_click(b):
    if not w_fwd.value.strip() or not w_rev.value.strip(): return
    btn_run.disabled = True
    log_output.clear_output()
    sim_output.clear_output()
    w_progress.value = 0
    w_progress.bar_style = 'info'
    
    # 啟動背景執行緒
    threading.Thread(target=run_blast_thread, args=(w_fwd.value, w_rev.value, w_org.value, w_db.value)).start()

# 綁定事件
btn_run.on_click(on_button_click)
w_slider.observe(lambda change: update_simulation(change['new']), names='value')

# 組合介面
ui = widgets.VBox([
    widgets.HTML("<h3>🧬 PCR 戰情室 (v9.0 動態熱模擬版)</h3>"),
    widgets.HBox([w_org, w_db]),
    widgets.HBox([w_fwd, w_rev]),
    btn_run,
    w_progress,
    log_output,
    widgets.HTML("<hr><h4>🎛️ 溫度模擬實驗室</h4>"),
    w_slider,
    sim_output
])

display(ui)
