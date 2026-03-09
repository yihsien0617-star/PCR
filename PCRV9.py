import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

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
    tm_f = mt.Tm_NN(Seq(f_seq.upper()))
    tm_r = mt.Tm_NN(Seq(r_seq.upper()))
    
    rh_f = NCBIWWW.qblast("blastn", database, f_seq, entrez_query=f'"{organism}"[Organism]', expect=50000, word_size=7)
    rec_f = NCBIXML.read(rh_f)
    
    rh_r = NCBIWWW.qblast("blastn", database, r_seq, entrez_query=f'"{organism}"[Organism]', expect=50000, word_size=7)
    rec_r = NCBIXML.read(rh_r)
    
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
    # 將圖表稍微加寬，以配合滿版的排版
    fig, ax = plt.subplots(figsize=(10, 5))
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
        if visible_count < 6: 
            ax.text(3.7, y, label, color=color, fontsize=10, va='center')
        visible_count += 1

    ax.set_xlim(0, 5)
    ax.set_ylim(-np.log10(3500), -np.log10(80))
    ax.set_xticks([1, 3])
    ax.set_xticklabels(['Marker', f'PCR @ {temp_val}°C'], color='black')
    ax.set_yticks([])
    ax.set_title(f"Virtual Gel Simulation: {organism}", color='white', pad=15)
    plt.tight_layout()
    return fig

# ================= Streamlit 介面佈局 =================

st.title("🧬 PCR 引子戰情室")
st.markdown("實驗室專用：引子設計驗證與動態熱模擬系統")

# Sidebar 設定區
with st.sidebar:
    st.header("⚙️ 參數設定")
    organism = st.selectbox("目標物種", ["Homo sapiens", "Mus musculus", "Rattus norvegicus"])
    database = st.selectbox("資料庫類型", ["refseq_rna", "nt"], index=0, help="mRNA用 refseq_rna, gDNA用 nt")
    enzyme = st.selectbox("酵素類型", ["Taq", "High-Fidelity"])
    st.info("💡 提示：輸入引子後按「開始分析」，完成後可調整下方溫度滑桿。")

# 主要輸入區
col1, col2 = st.columns(2)
with col1:
    f_seq = st.text_input("Forward Primer (5'->3')", placeholder="輸入序列...").strip()
with col2:
    r_seq = st.text_input("Reverse Primer (5'->3')", placeholder="輸入序列...").strip()

# 分析按鈕與邏輯
if st.button("🚀 開始分析 / 更新資料", type="primary", use_container_width=True):
    if not f_seq or not r_seq:
        st.error("❌ 請輸入完整的 Fwd 和 Rev 引子序列！")
    else:
        try:
            with st.spinner(f"正在連線 NCBI 資料庫比對 {organism}... (約需 1-2 分鐘，請保持網頁開啟)"):
                data_f, data_r, tm_f, tm_r = fetch_blast_data(f_seq, r_seq, organism, database)
                
                st.session_state['blast_result'] = {
                    'f': data_f, 'r': data_r, 
                    'tm_f': tm_f, 'tm_r': tm_r,
                    'done': True
                }
            st.success("✅ 資料下載完成！請向下捲動查看分析與圖表。")
        except Exception as e:
            st.error(f"連線或分析發生錯誤: {e}")

# ================= 結果顯示與熱模擬區 (滿版設計) =================

if 'blast_result' in st.session_state and st.session_state['blast_result']['done']:
    res = st.session_state['blast_result']
    
    st.divider()
    st.header("🎛️ 動態熱模擬實驗室")
    
    # 顯示基礎資訊
    c1, c2, c3 = st.columns(3)
    c1.metric("Forward Tm", f"{res['tm_f']:.1f} °C")
    c2.metric("Reverse Tm", f"{res['tm_r']:.1f} °C")
    
    # 溫度滑桿 (放寬到滿版)
    default_temp = min(res['tm_f'], res['tm_r']) - 5
    temp_val = st.slider("🌡️ 調整 Annealing 溫度 (°C) 以模擬 PCR 嚴格度：", min_value=40.0, max_value=75.0, value=default_temp, step=0.5)
    
    # --- 即時計算產物 ---
    products = []
    for f in res['f']:
        hit_tm_f = calculate_adjusted_tm(res['tm_f'], f['ident'])
        if hit_tm_f < temp_val: continue

        for r in res['r']:
            if f['id'] == r['id'] and f['str'] != r['str']:
                hit_tm_r = calculate_adjusted_tm(res['tm_r'], r['ident'])
                if hit_tm_r < temp_val: continue
                
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
    
    # ================= 改變排版：先文字，後圖片 =================
    
    st.subheader(f"📝 診斷報告 (目前溫度: {temp_val}°C)")
    
    # --- 1. 文字敘述區 (滿版) ---
    if not products:
        st.warning("❄️ 此溫度下無任何 PCR 產物 (設定溫度過高，引子無法結合)。")
    else:
        # 主要產物與建議
        target = products[0]
        if target['f_ident'] == 100 and target['r_ident'] == 100:
            st.success(f"**🎯 預期主要產物發現： {target['size']} bp**")
            enz_name, ann, ext = calculate_pcr_conditions(target['size'], res['tm_f'], res['tm_r'], enzyme)
            st.info(f"🧪 **建議機器設定 ({enz_name})** ➜ Annealing: **{ann:.1f}°C** | Extension: **{ext} 秒**")
        else:
            st.warning(f"⚠️ 發現主要產物 ({target['size']} bp)，但並非 100% 完美結合，請注意非專一性風險。")

        # 列出所有產物細節 (字數顯示放寬到 80 字元)
        with st.expander(f"🔎 查看詳細產物列表 (共 {len(products)} 條)", expanded=True):
            for p in products:
                icon = "✅" if (p['f_ident']==100 and p['r_ident']==100) else "⚠️"
                st.write(f"{icon} **{p['size']} bp** ➜ {p['gene'][:80]}...")
    
    st.markdown("<br>", unsafe_allow_html=True) # 增加一點留白
    
    # --- 2. 視覺化圖表區 (滿版置中) ---
    if products:
        st.subheader("📊 虛擬電泳圖")
        fig = draw_gel_plot(products, temp_val, organism)
        # 顯示圖表，不強制填滿容器，保持適當比例
        st.pyplot(fig, use_container_width=False)
