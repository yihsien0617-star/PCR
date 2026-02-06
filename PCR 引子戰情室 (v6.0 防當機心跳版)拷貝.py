Python 3.14.0 (tags/v3.14.0:ebf955d, Oct  7 2025, 10:15:03) [MSC v.1944 64 bit (AMD64)] on win32
Enter "help" below or click "Help" above for more information.
>>> # @title 🧬 PCR 戰情室 (v8.0 文字藝術版)
... # 放棄多執行緒與圖形庫，改用純文字輸出，保證看得到結果！
... 
... import ipywidgets as widgets
... from IPython.display import display, clear_output
... import sys
... import time
... 
... # 確保 Biopython 已載入
... try:
...     import Bio
...     from Bio.Blast import NCBIWWW
...     from Bio.Blast import NCBIXML
...     from Bio.SeqUtils import MeltingTemp as mt
...     from Bio.Seq import Seq
... except ImportError:
...     print("⚠️ 尚未偵測到 Biopython，請在 Anaconda Prompt 輸入 'pip install biopython' 安裝。")
... 
... # ================= 核心邏輯 =================
... 
... def calculate_pcr_conditions(product_size, tm_f, tm_r, enzyme):
...     tm_min = min(tm_f, tm_r)
...     anneal_temp = tm_min - 3 
...     if enzyme == "Taq":
...         rate, name = 60, "General Taq"
...     else:
...         rate, name = 30, "High-Fidelity"
...     ext_time = max(30, int((product_size / 1000) * rate))
...     return f"95°C 30s -> {anneal_temp:.1f}°C 30s -> 72°C {ext_time}s"
... 
... def draw_ascii_gel(products):
...     """
...     使用純文字畫出電泳圖 (ASCII Art)
...     保證所有電腦都能顯示
...     """
...     print("\n========= [ 虛擬電泳圖 (文字版) ] =========")
    print("Marker (bp)        您的 PCR 產物")
    print("|                  |")
    
    # 定義 Marker 刻度
    markers = [3000, 2000, 1500, 1000, 800, 500, 300, 200, 100]
    
    # 建立畫布
    lines = []
    
    for m in markers:
        # 檢查在這個 Marker 附近有沒有產物
        found_p = []
        for p in products:
            # 如果產物大小在 Marker 的 +/- 15% 範圍內，視為同一高度
            if m * 0.85 <= p['size'] <= m * 1.15:
                found_p.append(p)
        
        # 繪製 Marker 線
        line_str = f" - {m:<4} ----------"
        
        if found_p:
            # 畫出產物
            p_str = ""
            for p in found_p:
                if p == products[0]: # 主要產物
                    p_str += f" [█ {p['size']}bp 預期產物] "
                else:
                    p_str += f" [x {p['size']}bp 雜訊] "
            line_str += p_str
        else:
             line_str += " |"
             
        print(line_str)
        print(" |                  |") # 間隔空行
        
    print("===========================================\n")

# ================= 介面與執行 =================

style = {'description_width': 'initial'}
layout = widgets.Layout(width='98%')

w_org = widgets.Dropdown(options=[('人類 (Homo sapiens)', 'Homo sapiens'), ('小鼠 (Mus musculus)', 'Mus musculus')], value='Homo sapiens', description='1. 物種:', style=style)
w_db = widgets.Dropdown(options=[('mRNA (RefSeq RNA)', 'refseq_rna'), ('Genomic DNA (nt)', 'nt')], value='refseq_rna', description='2. 資料庫:', style=style)
w_fwd = widgets.Text(placeholder='輸入 Forward 引子', description='Fwd:', style=style, layout=layout)
w_rev = widgets.Text(placeholder='輸入 Reverse 引子', description='Rev:', style=style, layout=layout)
btn_run = widgets.Button(description="開始分析 (請耐心等待)", button_style='success', icon='play', layout=widgets.Layout(width='100%'))
output_area = widgets.Output()

def run_analysis_sync(b):
    """同步執行模式 (會暫時卡住畫面，但保證輸出)"""
    
    # 1. 鎖定按鈕
    btn_run.disabled = True
    btn_run.description = "分析中... (畫面凍結是正常的)"
    
    f_seq = w_fwd.value.strip()
    r_seq = w_rev.value.strip()
    organism = w_org.value
    database = w_db.value
    
    output_area.clear_output()
    
    with output_area:
        if not f_seq or not r_seq:
            print("❌ 錯誤：請輸入完整的引子序列。")
            btn_run.disabled = False
            btn_run.description = "開始分析"
            return

        print(f"🚀 [系統啟動] 目標: {organism} | 資料庫: {database}")
        print("------------------------------------------------")
        
        # Step 1: Tm
        try:
            tm_f = mt.Tm_NN(Seq(f_seq.upper()))
            tm_r = mt.Tm_NN(Seq(r_seq.upper()))
            print(f"1. [序列檢查] ✅ 通過")
            print(f"   - Fwd Tm: {tm_f:.1f}°C")
            print(f"   - Rev Tm: {tm_r:.1f}°C")
        except Exception as e:
            print(f"❌ 序列錯誤: {e}")
            btn_run.disabled = False
            return

        # Step 2: Blast Fwd
        print(f"2. [NCBI 連線] 正在比對 Forward 引子... (需時約 30-60秒)")
        try:
            rh_f = NCBIWWW.qblast("blastn", database, f_seq, entrez_query=f'"{organism}"[Organism]', expect=20000, word_size=7)
            rec_f = NCBIXML.read(rh_f)
            print("   - Forward 比對完成！")
        except Exception as e:
            print(f"❌ 連線失敗 (Fwd): {e}")
            btn_run.disabled = False
            return

        # Step 3: Blast Rev
        print(f"3. [NCBI 連線] 正在比對 Reverse 引子... (需時約 30-60秒)")
        try:
            rh_r = NCBIWWW.qblast("blastn", database, r_seq, entrez_query=f'"{organism}"[Organism]', expect=20000, word_size=7)
            rec_r = NCBIXML.read(rh_r)
            print("   - Reverse 比對完成！")
        except Exception as e:
            print(f"❌ 連線失敗 (Rev): {e}")
            btn_run.disabled = False
            return

        # Step 4: Analyze
        print(f"4. [數據整合] 正在計算 PCR 產物...")
        products = []
        f_locs = [{'id': a.accession, 's': h.sbjct_start, 'str': h.frame[1], 't': a.title} for a in rec_f.alignments for h in a.hsps]
        r_locs = [{'id': a.accession, 's': h.sbjct_start, 'str': h.frame[1]} for a in rec_r.alignments for h in a.hsps]

        for f in f_locs:
            for r in r_locs:
                if f['id'] == r['id'] and f['str'] != r['str']:
                    dist = abs(r['s'] - f['s']) + 1
                    if 50 < dist < 5000 and not any(p['size'] == dist for p in products):
                        products.append({'size': dist, 'gene': f['t']})
        
        products.sort(key=lambda x: x['size'])

        # 輸出結果
        print("\n" + "█"*20 + " 分析結果報告 " + "█"*20)
        
        if not products:
            print("\n❌ 結果: 未發現任何 PCR 產物！")
            print("   > 可能原因: 引子距離太遠、跨越 Intron (若選mRNA資料庫)、或特異性太低。")
        else:
            target = products[0]
            cycle = calculate_pcr_conditions(target['size'], tm_f, tm_r, "Taq")
            
            print(f"\n✅ [主要產物] 長度: {target['size']} bp")
            print(f"   > 基因: {target['gene'][:60]}...")
            print(f"   > 建議 PCR 設定: {cycle}")
            
            if len(products) > 1:
                print(f"\n⚠️ [潛在風險] 發現 {len(products)-1} 個非專一性產物 (Off-target):")
                for p in products[1:4]:
                    print(f"   - {p['size']} bp ({p['gene'][:40]}...)")
            
            # 呼叫文字畫圖
            draw_ascii_gel(products)

    # 恢復按鈕
    btn_run.disabled = False
    btn_run.description = "開始分析"

btn_run.on_click(run_analysis_sync)

ui = widgets.VBox([
    widgets.HTML("<h3>🧬 PCR 戰情室 (文字極速版)</h3>"),
    widgets.HBox([w_org, w_db]),
    widgets.HBox([w_fwd, w_rev]),
    widgets.HTML("<hr>"),
    btn_run,
    output_area
])

