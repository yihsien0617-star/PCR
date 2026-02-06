import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import time

# 設定頁面資訊
st.set_page_config(page_title="PCR 引子戰情室", page_icon="🧬", layout="wide")

# ================= 核心計算邏輯 =================

def calculate_adjusted_tm(perfect_tm, identity_percent):
    """估算不完美結合時的 Tm 值"""
    mismatch_percent = 100 - identity_percent
    drop = mismatch_percent * 1.2 
    return perfect_tm - drop

def calculate_pcr_conditions(product_size, tm_f, tm_r, enzyme):
    tm_min = min(tm_f, tm_r)
    anneal_temp = tm_min - 3 
    if enzyme == "Taq":
        rate, name = 60, "General Taq"
    else:
        rate, name = 30, "High-Fidelity"
    ext_time = max(30, int((product_size / 1000) * rate))
    return name, anneal_temp, ext_time

# 使用 Streamlit 的快取功能，避免重複連線 NCBI
@st.cache_data(show_spinner=False)
def fetch_blast_data(f_seq, r_seq, organism, database):
    """連線 NCBI 抓取資料 (會被快取)"""
    
    # 1. 計算 Tm
    tm_f = mt.Tm_NN(Seq(f_seq.upper()))
    tm_r = mt.Tm_NN(Seq(r_seq.upper()))
    
    # 2. BLAST Forward
    rh_f = NCBIWWW.qblast("blastn", database, f_seq, entrez_query=f'"{organism}"[Organism]', expect=50000, word_size=7)
    rec_f = NCBIXML.read(rh_f)
    
    # 3. BLAST Reverse
    rh_r = NCBIWWW.qblast("blastn", database, r_seq, entrez_query=f'"{organism}"[Organism]', expect=50000, word_size=7)
    rec_r = NCBIXML.read(rh_r)
    
    # 4. 解析資料
    parsed_f = []
    for a in rec_f.alignments:
        for h in a.hsps:
            ident = (h.identities / h.align_length) * 100
            parsed_f.append({'id': a.accession, 's': h.sbjct_start, 'str': h.frame[1], 't': a.title, 'ident': ident})
            
    parsed_r = []
    for a in rec_r.alignments:
        for h in a.hsps:
            ident = (h.identities / h.align_length) * 100
            parsed_r.append({'id': a.accession, 's': h.sbjct_start, 'str': h.frame[1], 't': a.title, 'ident': ident})
            
    return parsed_f, parsed_r, tm_f, tm_r

def draw_gel_plot(products, temp_val, organism):
    """繪製 Matplotlib 圖表"""
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.set_facecolor('black')
    
    # Marker
    ladder = [100, 200, 300, 400, 500, 600, 800, 1000, 1500, 2000, 3000]
    for b in ladder:
        y = -np.log10(b)
        ax.hlines(y, 0.5, 1.5, colors='white', alpha=0.5)
        ax.text(0.1, y, f"{b}", color='white', fontsize=8, va='center')

    # Products
    visible_count = 0
    for i, p in enumerate(products):
        size = p['size']
        y = -np.log10(size)
        
        is_target = (p['f_ident'] == 100 and p['r_ident'] == 100)
        
        if is_target:
            color, alpha, lw = '#00FF00', 1.0, 3.5
            label = f"{size}bp (Target)"
        else:
            color, alpha, lw = 'yellow', 0.6, 2
            label = f"{size}bp (Non-specific)"
            
        ax.hlines(y, 2.5, 3.5, colors=color, alpha=alpha, linewidth=lw)
        if visible_count < 6: # 避免太多文字重疊
            ax.text(3.7, y, label, color=color, fontsize=9, va='center')
        visible_count += 1

    ax.set_xlim(0, 5)
    ax.set_ylim(-np.log10(3500), -np.log10(80))
    ax.set_xticks([1, 3])
    ax.set_xticklabels(['Marker', f'PCR @ {temp_val}°C'])
    ax.set_yticks([])
    ax.set_title(f"In-Silico PCR Simulation: {organism}", color='white')
    plt.tight_layout()
    return fig

# ================= Streamlit 介面佈局 =================

st.title("🧬 PCR 引子戰情室 (Streamlit Pro)")
st.markdown("實驗室專用引子設計驗證與熱模擬系統")

# Sidebar 設定區
with st.sidebar:
    st.header("⚙️ 參數設定")
    organism = st.selectbox("目標物種", ["Homo sapiens", "Mus musculus", "Rattus norvegicus"])
    database = st.selectbox("資料庫類型", ["refseq_rna", "nt"], index=0, help="mRNA用 refseq_rna, gDNA用 nt")
    enzyme = st.selectbox("酵素類型", ["Taq", "High-Fidelity"])
    st.info("💡 提示：輸入引子後按「開始分析」，完成後可調整溫度滑桿。")

# 主要輸入區
col1, col2 = st.columns(2)
with col1:
    f_seq = st.text_input("Forward Primer (5'->3')", placeholder="輸入序列...").strip()
with col2:
    r_seq = st.text_input("Reverse Primer (5'->3')", placeholder="輸入序列...").strip()

# 分析按鈕與邏輯
if st.button("🚀 開始分析 / 更新資料", type="primary"):
    if not f_seq or not r_seq:
        st.error("❌ 請輸入完整的 Fwd 和 Rev 引子序列！")
    else:
        try:
            with st.spinner(f"正在連線 NCBI 資料庫比對 {organism}... (這可能需要 1-2 分鐘)"):
                # 呼叫快取函數
                data_f, data_r, tm_f, tm_r = fetch_blast_data(f_seq, r_seq, organism, database)
                
                # 將結果存入 Session State (讓滑桿互動時不用重跑)
                st.session_state['blast_result'] = {
                    'f': data_f, 'r': data_r, 
                    'tm_f': tm_f, 'tm_r': tm_r,
                    'done': True
                }
            st.success("✅ 分析完成！請查看下方模擬結果。")
        except Exception as e:
            st.error(f"發生錯誤: {e}")

# ================= 結果顯示與熱模擬區 =================

if 'blast_result' in st.session_state and st.session_state['blast_result']['done']:
    res = st.session_state['blast_result']
    
    st.divider()
    st.subheader("🎛️ 動態熱模擬實驗室")
    
    # 顯示基礎資訊
    c1, c2, c3 = st.columns(3)
    c1.metric("Fwd Tm", f"{res['tm_f']:.1f} °C")
    c2.metric("Rev Tm", f"{res['tm_r']:.1f} °C")
    
    # 溫度滑桿
    default_temp = min(res['tm_f'], res['tm_r']) - 5
    temp_val = st.slider("調整 Annealing 溫度 (°C)", min_value=40.0, max_value=75.0, value=default_temp, step=0.5)
    
    # --- 即時計算產物 ---
    products = []
    f_hits = res['f']
    r_hits = res['r']
    
    for f in f_hits:
        hit_tm_f = calculate_adjusted_tm(res['tm_f'], f['ident'])
        if hit_tm_f < temp_val: continue # 溫度過高，引子脫落

        for r in r_hits:
            if f['id'] == r['id'] and f['str'] != r['str']:
                hit_tm_r = calculate_adjusted_tm(res['tm_r'], r['ident'])
                if hit_tm_r < temp_val: continue # 溫度過高，引子脫落
                
                dist = abs(r['s'] - f['s']) + 1
                if 50 < dist < 5000:
                    if not any(p['size'] == dist for p in products):
                        products.append({
                            'size': dist, 
                            'gene': f['t'],
                            'f_ident': f['ident'], 
                            'r_ident': r['ident']
                        })
    
    products.sort(key=lambda x: x['size'])
    
    # --- 顯示結果 ---
    col_res, col_img = st.columns([1, 2])
    
    with col_res:
        st.markdown(f"**模擬結果: 發現 {len(products)} 條產物**")
        if not products:
            st.warning("❄️ 此溫度下無產物 (引子無法結合)")
        else:
            # 找出主要產物給建議
            target = products[0]
            if target['f_ident'] == 100 and target['r_ident'] == 100:
                st.success(f"✅ 主要產物: {target['size']} bp")
                enz_name, ann, ext = calculate_pcr_conditions(target['size'], res['tm_f'], res['tm_r'], enzyme)
                st.info(f"🧪 **建議條件 ({enz_name})**\n\n"
                        f"- Anneal: {ann:.1f}°C\n"
                        f"- Extension: {ext}s")
            else:
                st.warning("⚠️ 主要產物似乎非完美結合")

            # 列出所有產物
            with st.expander("查看詳細產物列表", expanded=True):
                for p in products:
                    icon = "✅" if (p['f_ident']==100 and p['r_ident']==100) else "⚠️"
                    st.write(f"{icon} **{p['size']} bp** - {p['gene'][:30]}...")

    with col_img:
        if products:
            fig = draw_gel_plot(products, temp_val, organism)
            st.pyplot(fig)
