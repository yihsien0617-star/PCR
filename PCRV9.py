"""
🧬 PCR 戰情室 — Streamlit 版 (v10.0)

在 Streamlit Community Cloud 上執行：
    1. 把本檔與 requirements.txt 一起放進 GitHub repo
    2. Streamlit Cloud 的 Main file path 填 streamlit_app.py

與 Colab/ipywidgets 版的差異
    - 移除 %matplotlib inline（Jupyter magic，不是合法 Python 語法）
    - 移除 ipywidgets，改用原生 Streamlit 元件
    - 移除背景執行緒；BLAST 改用 st.cache_data 快取，只有第一次要等
    - 反應條件改以 frozen dataclass 傳遞，不再依賴全域 widget

科學計算的修正全部保留：
    - Tm：SantaLucia & Hicks 2004 (DNA_NN4) + Owczarzy 2008 鹽校正 (saltcorr=7)
    - 錯配 Tm：由 BLAST 比對字串以 IMM/TMM/DE 熱力學表直接計算
    - identity 分母 = 引子全長（不是 HSP 長度）
    - amplicon 長度 = 兩條引子 5' 端距離（HSP 截斷時會外推）
    - 3' 端朝向檢查 + 3' 末端錯配封鎖
    - F-F / R-R 自我配對
"""

from __future__ import annotations

import io
from dataclasses import dataclass, replace

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st

from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt

st.set_page_config(page_title="PCR 戰情室", page_icon="🧬", layout="wide")

MIN_PRODUCT, MAX_PRODUCT = 50, 5000
THREE_PRIME_CRITICAL = 2      # 3' 末 2 nt 錯配 -> 幾乎不延伸
THREE_PRIME_WINDOW = 5        # 3' 末 5 nt 內的錯配額外註記


# =============================================================================
#  反應條件
# =============================================================================

@dataclass(frozen=True)
class Cond:
    """PCR 反應條件。frozen 才能當 st.cache_data 的 key。"""
    mg: float = 1.5           # mM
    dntp: float = 0.2         # mM (each)
    mono: float = 50.0        # mM K+
    tris: float = 10.0        # mM
    primer_nM: float = 200.0  # nM
    offset: float = 0.0       # °C 校正

    def key(self):
        return (self.mg, self.dntp, self.mono, self.tris,
                self.primer_nM, self.offset)

    def label(self):
        return (f"{self.mono:g} mM K⁺ / {self.tris:g} mM Tris / "
                f"{self.mg:.2f} mM Mg²⁺ / {self.dntp:g} mM dNTP / "
                f"{self.primer_nM:g} nM primer / offset {self.offset:+.1f}°C")


def nn_kwargs(cond: Cond) -> dict:
    """Bio.SeqUtils.MeltingTemp 的參數。

    saltcorr=7 (Owczarzy 2008) 是 Biopython 中唯一同時正確處理 Mg²⁺ 與
    dNTP 螯合的鹽校正法。用預設的 saltcorr=5 / Mg=0 會系統性低估 Tm 約 6-8°C。
    注意：Mg <= dNTP 時公式退回只算一價鹽，所以該區間 Tm 不隨 Mg 變動。
    """
    return dict(
        dnac1=cond.primer_nM, dnac2=0.0, selfcomp=False,
        Na=0.0, K=cond.mono, Tris=cond.tris,
        Mg=cond.mg, dNTPs=cond.dntp,
        nn_table=mt.DNA_NN4,
        tmm_table=mt.DNA_TMM1,
        imm_table=mt.DNA_IMM1,
        de_table=mt.DNA_DE1,
        saltcorr=7, strict=False,
    )


# =============================================================================
#  熱力學
# =============================================================================

def clean_seq(s: str) -> str:
    return "".join(ch for ch in s.upper() if not ch.isspace())


def valid_seq(s: str) -> bool:
    return bool(s) and set(s) <= set("ATGC")


def tm_perfect(seq: str, cond: Cond):
    seq = clean_seq(seq)
    if not valid_seq(seq):
        return None
    return mt.Tm_NN(Seq(seq), **nn_kwargs(cond)) + cond.offset


def tm_duplex(query_aln: str, sbjct_aln: str, cond: Cond):
    """由 BLAST 比對字串算含錯配的真實雙股 Tm。

    BLAST 的 hsp.query / hsp.sbjct 是「同義股」對齊（字母相同 = match），
    引子實際結合的是 sbjct 的『互補』股，左到右恰為 3'->5'，
    正是 Biopython c_seq 參數要求的方向。
    """
    q, s = query_aln.upper(), sbjct_aln.upper()
    if len(q) != len(s):
        return None, "length-mismatch"

    if "-" in q or "-" in s:
        # 含 indel：Tm_NN 不支援，退回保守經驗式
        bare = q.replace("-", "")
        base = tm_perfect(bare, cond)
        if base is None:
            return None, "gap-fallback"
        ident = sum(1 for a, b in zip(q, s) if a == b and a != "-")
        return base - (1 - ident / max(len(bare), 1)) * 100 * 1.2, "gap-fallback"

    try:
        c = str(Seq(s).complement())          # 3'->5' 的模板股
        return mt.Tm_NN(Seq(q), c_seq=Seq(c), **nn_kwargs(cond)) + cond.offset, "NN"
    except Exception:
        return None, "nn-failed"


def recommend_ta(tm_f, tm_r):
    """Taq / hot-start Taq 慣例：Ta = 較低的 Tm − 5，上限 72°C。"""
    if tm_f is None or tm_r is None:
        return None
    return min(min(tm_f, tm_r) - 5.0, 72.0)


def solve_mg(target_ta: float, f_seq: str, r_seq: str, cond: Cond):
    """由「實驗上跑得乾淨的 Ta」反解有效 Mg²⁺（二分法）。

    Tm 對 Mg²⁺ 在 Mg > dNTP 的區間內單調遞增，故二分法收斂。
    下界設在 dNTP+0.3 以避開 Owczarzy 公式的分支切換點。
    回傳 (mg 或 None, 說明字串)。
    """
    target_tm = target_ta + 5.0

    def limiting_tm(mg):
        c = replace(cond, mg=mg)
        return min(tm_perfect(f_seq, c), tm_perfect(r_seq, c))

    lo, hi = cond.dntp + 0.3, 10.0
    f_lo, f_hi = limiting_tm(lo), limiting_tm(hi)

    if target_tm < f_lo:
        return None, (f"實驗 Ta 偏低：即使 Mg²⁺ 只有 {lo:.1f} mM，限速引子 Tm 仍有 "
                      f"{f_lo:.1f}°C（對應 Ta {f_lo-5:.1f}°C）。"
                      f"差額請改用 offset = {target_tm - f_lo:+.1f}°C 補。")
    if target_tm > f_hi:
        return None, (f"實驗 Ta 偏高：即使 Mg²⁺ 拉到 {hi:.0f} mM，限速引子 Tm 只有 "
                      f"{f_hi:.1f}°C（對應 Ta {f_hi-5:.1f}°C）。"
                      f"差額請改用 offset = {target_tm - f_hi:+.1f}°C 補。")

    for _ in range(60):
        mid = (lo + hi) / 2.0
        if limiting_tm(mid) < target_tm:
            lo = mid
        else:
            hi = mid
    mg = (lo + hi) / 2.0
    ok = 0.8 <= mg <= 4.0
    note = ("（落在一般 master mix 範圍內 ✔）" if ok else
            "（⚠️ 超出一般 master mix 的合理範圍，差異來源可能不是 Mg²⁺，"
            "而是酶種類、master mix 內的增強劑或引子純度；建議改用 offset 校正）")
    return mg, f"有效 Mg²⁺ ≈ **{mg:.2f} mM** {note}"


# =============================================================================
#  3' 端規則
# =============================================================================

def three_prime_status(q_aln: str, s_aln: str, q_end: int, plen: int):
    """BLAST 會在錯配處截斷比對，所以 q_end < plen 通常就代表 3' 端錯配。"""
    if q_end < plen:
        return False, f"3'端未比對到（缺 {plen - q_end} nt）"

    tq, ts = q_aln.upper()[-THREE_PRIME_WINDOW:], s_aln.upper()[-THREE_PRIME_WINDOW:]
    mism = [i for i, (a, b) in enumerate(zip(tq, ts)) if a != b]
    if not mism:
        return True, "3'端完全配對"
    dist = min(len(tq) - 1 - i for i in mism)     # 0 = 最末鹼基
    if dist < THREE_PRIME_CRITICAL:
        return False, f"3'末端第 {dist+1} 位錯配 → 無法延伸"
    return True, f"3'端 -{dist+1} 位錯配（效率下降）"


# =============================================================================
#  BLAST（快取）
# =============================================================================

@st.cache_data(show_spinner=False, ttl=24 * 3600, max_entries=64)
def blast_xml(seq: str, organism: str, database: str) -> str:
    """向 NCBI 查詢並回傳 XML 字串。

    快取 key = (seq, organism, database)，所以調整溫度滑桿或反應條件
    都不會重新查詢——只有換引子/物種/資料庫才會重跑。
    """
    handle = NCBIWWW.qblast(
        "blastn", database, seq,
        entrez_query=f'"{organism}"[Organism]',
        expect=1000, word_size=7,
        hitlist_size=500, megablast=False,
    )
    try:
        return handle.read()
    finally:
        handle.close()


@st.cache_data(show_spinner=False, max_entries=128)
def parse_hits(xml: str, primer: str, tag: str, cond_key: tuple) -> list:
    """解析 XML 並算好每個結合位點的 Tm。"""
    cond = Cond(*cond_key)
    plen = len(primer)
    record = NCBIXML.read(io.StringIO(xml))

    hits = []
    for a in record.alignments:
        for h in a.hsps:
            s_start, s_end = h.sbjct_start, h.sbjct_end
            strand = 1 if s_end >= s_start else -1    # 負股 hit-from > hit-to

            ident = h.identities / plen * 100.0       # 分母必須是引子全長
            coverage = h.align_length / plen * 100.0

            # 把 HSP 外推到引子兩端（比對欄位左到右 = query 5'->3'）
            five_prime = s_start - strand * (h.query_start - 1)
            three_prime = s_end + strand * (plen - h.query_end)

            tm, method = tm_duplex(h.query, h.sbjct, cond)
            ext_ok, ext_note = three_prime_status(h.query, h.sbjct,
                                                  h.query_end, plen)
            hits.append({
                "tag": tag, "id": a.accession, "title": a.title,
                "strand": strand,
                "five_prime": five_prime, "three_prime": three_prime,
                "ident": ident, "coverage": coverage,
                "tm": tm, "tm_method": method,
                "ext_ok": ext_ok, "ext_note": ext_note,
                "perfect": ident >= 99.99 and coverage >= 99.99,
            })
    return hits


# =============================================================================
#  產物模擬
# =============================================================================

def pair_products(hits_f, hits_r, temp: float, allow_self: bool = True):
    """擴增條件：同一條序列、一正一負股、兩條 3' 端『面對面』、
    且在該溫度下都還結合著且 3' 端可延伸。
    amplicon 長度 = 兩條引子 5' 端之間的距離。"""
    def alive(h):
        return h["ext_ok"] and h["tm"] is not None and h["tm"] >= temp

    by_acc = {}
    for h in list(hits_f) + list(hits_r):
        if alive(h):
            by_acc.setdefault(h["id"], []).append(h)

    products, seen = [], set()
    for acc, group in by_acc.items():
        plus_hits = [h for h in group if h["strand"] == 1]
        minus_hits = [h for h in group if h["strand"] == -1]
        for plus in plus_hits:
            for minus in minus_hits:
                if not allow_self and plus["tag"] == minus["tag"]:
                    continue
                if plus["three_prime"] >= minus["three_prime"]:
                    continue                       # 3' 端背對背 → 不擴增
                size = minus["five_prime"] - plus["five_prime"] + 1
                if not (MIN_PRODUCT <= size <= MAX_PRODUCT):
                    continue
                key = (acc, plus["five_prime"], minus["five_prime"])
                if key in seen:
                    continue
                seen.add(key)
                products.append({
                    "size": size, "id": acc, "title": plus["title"],
                    "pair": f"{plus['tag']}/{minus['tag']}",
                    "f_tm": plus["tm"], "r_tm": minus["tm"],
                    "f_ident": plus["ident"], "r_ident": minus["ident"],
                    "perfect": plus["perfect"] and minus["perfect"]
                               and plus["tag"] != minus["tag"],
                    "margin": min(plus["tm"], minus["tm"]) - temp,
                })
    products.sort(key=lambda p: (not p["perfect"], p["size"]))
    return products


# =============================================================================
#  電泳圖
# =============================================================================

def draw_gel(products, temp: float, organism: str):
    fig, ax = plt.subplots(figsize=(7.5, 5.2), dpi=130)
    fig.patch.set_facecolor("#0E1117")
    ax.set_facecolor("black")

    for b in (100, 200, 300, 400, 500, 600, 800, 1000, 1500, 2000, 3000, 5000):
        y = -np.log10(b)
        ax.hlines(y, 0.5, 1.5, colors="white", alpha=0.45, linewidth=1.4)
        ax.text(0.38, y, f"{b}", color="white", fontsize=6.5, va="center", ha="right")

    order = sorted(range(len(products)),
                   key=lambda i: (not products[i]["perfect"], -products[i]["margin"]))
    labelled = 0
    for i in order:
        p = products[i]
        y = -np.log10(p["size"])
        if p["perfect"]:
            color, alpha, lw, tag = "#00FF41", 1.0, 4.0, "Target"
        else:
            alpha = float(np.clip(0.25 + p["margin"] / 20.0, 0.2, 0.85))
            color, lw, tag = "#FFD400", 2.0, f"Non-spec {p['pair']}"
        ax.hlines(y, 2.4, 3.5, colors=color, alpha=alpha, linewidth=lw)
        if labelled < 6:
            ax.text(3.65, y, f"{p['size']} bp ({tag})", color=color,
                    fontsize=7.5, va="center", alpha=max(alpha, 0.7))
            labelled += 1

    ax.set_xlim(0, 6.4)
    ax.set_ylim(-np.log10(MAX_PRODUCT * 1.15), -np.log10(MIN_PRODUCT * 0.85))
    ax.set_xticks([1, 2.95])
    ax.set_xticklabels(["Marker", f"PCR @ {temp:.1f}°C"], color="white", fontsize=8)
    ax.tick_params(colors="white")
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_color("#333333")
    ax.set_title(f"In-Silico PCR  ({organism})", color="white", fontsize=10)
    fig.tight_layout()
    return fig


# =============================================================================
#  側邊欄：反應條件
# =============================================================================

st.sidebar.header("🧪 反應條件")

mg_mode = st.sidebar.radio(
    "Mg²⁺ 濃度",
    ["① 已知濃度", "② 未知 → 區間估計", "③ 未知 → 由實驗 Ta 反推"],
    index=1,
    help="市售 master mix 幾乎都不公開 Mg²⁺ 濃度，預設用區間估計。",
)

mg_lo = mg_hi = None
if mg_mode.startswith("①"):
    mg_point = st.sidebar.number_input("Mg²⁺ (mM)", 0.0, 10.0, 1.5, 0.1)
elif mg_mode.startswith("②"):
    c1, c2 = st.sidebar.columns(2)
    mg_lo = c1.number_input("下限 (mM)", 0.0, 10.0, 1.0, 0.1)
    mg_hi = c2.number_input("上限 (mM)", 0.0, 10.0, 3.0, 0.1)
    if mg_lo > mg_hi:
        mg_lo, mg_hi = mg_hi, mg_lo
    mg_point = (mg_lo + mg_hi) / 2.0
    st.sidebar.caption("2X Taq master mix 稀釋成 1X 後，Mg²⁺ 多落在 1.5–2.5 mM；"
                       "1.0–3.0 是保守涵蓋範圍。Tm 與 Ta 會以區間顯示。")
else:
    mg_point = st.session_state.get("mg_solved", 1.5)
    st.sidebar.metric("目前採用的有效 Mg²⁺", f"{mg_point:.2f} mM")
    st.sidebar.caption("在下方「Mg²⁺ 反推」區塊輸入一組實測 Ta 即可求解。")

st.sidebar.divider()
c1, c2 = st.sidebar.columns(2)
dntp = c1.number_input("dNTP each (mM)", 0.0, 2.0, 0.2, 0.05)
mono = c2.number_input("K⁺ (mM)", 0.0, 200.0, 50.0, 5.0)
c3, c4 = st.sidebar.columns(2)
tris = c3.number_input("Tris (mM)", 0.0, 100.0, 10.0, 1.0)
primer_nM = c4.number_input("Primer (nM)", 10.0, 2000.0, 200.0, 10.0)
offset = st.sidebar.number_input(
    "校正 offset (°C)", -15.0, 15.0, 0.0, 0.5,
    help="模型與實驗的系統性差異。用一組已知可行的 Ta 校正後，其他引子就會準。")

cond = Cond(mg=mg_point, dntp=dntp, mono=mono, tris=tris,
            primer_nM=primer_nM, offset=offset)

st.sidebar.divider()
organism = st.sidebar.selectbox("物種", ["Homo sapiens", "Mus musculus"],
                                format_func=lambda x: {"Homo sapiens": "人類 (Homo sapiens)",
                                                       "Mus musculus": "小鼠 (Mus musculus)"}[x])
database = st.sidebar.selectbox("資料庫", ["refseq_rna", "nt"],
                                format_func=lambda x: {"refseq_rna": "mRNA (refseq_rna)",
                                                       "nt": "Genomic (nt)"}[x])


# =============================================================================
#  主畫面
# =============================================================================

st.title("🧬 PCR 戰情室")
st.caption("In-silico PCR ・ 熱力學 Tm ・ 非專一性產物模擬　|　"
           "Tm 模型：SantaLucia & Hicks 2004 + Owczarzy 2008 鹽校正")

c1, c2 = st.columns(2)
f_seq = clean_seq(c1.text_input("Forward primer (5'→3')",
                                "AGGTCAAAGAGGCTGCTTGG"))
r_seq = clean_seq(c2.text_input("Reverse primer (5'→3')",
                                "AACTGCATGGAATTGGTTGAC"))

bad = [n for n, s in (("Forward", f_seq), ("Reverse", r_seq)) if s and not valid_seq(s)]
if bad:
    st.error(f"{'、'.join(bad)} 含有 A/T/G/C 以外的字元，無法計算 Tm。")
    st.stop()
if not (f_seq and r_seq):
    st.info("請輸入兩條引子序列。")
    st.stop()


# ---------- Tm 面板 ----------
st.subheader("🌡️ Tm 與建議 Annealing 溫度")

tm_f = tm_perfect(f_seq, cond)
tm_r = tm_perfect(r_seq, cond)
ta = recommend_ta(tm_f, tm_r)


def gc_pct(s):
    return 100.0 * (s.count("G") + s.count("C")) / max(len(s), 1)


m1, m2, m3, m4 = st.columns(4)
m1.metric("Forward Tm", f"{tm_f:.1f} °C", f"{len(f_seq)} nt ・ GC {gc_pct(f_seq):.0f}%",
          delta_color="off")
m2.metric("Reverse Tm", f"{tm_r:.1f} °C", f"{len(r_seq)} nt ・ GC {gc_pct(r_seq):.0f}%",
          delta_color="off")
m3.metric("ΔTm", f"{abs(tm_f - tm_r):.1f} °C",
          "✔ 平衡" if abs(tm_f - tm_r) <= 5 else "⚠️ >5°C 建議重新設計",
          delta_color="off")
m4.metric("建議 Ta", f"{ta:.1f} °C", "Taq 慣例 Tm−5", delta_color="off")

if mg_lo is not None:
    lo_c, hi_c = replace(cond, mg=mg_lo), replace(cond, mg=mg_hi)
    ta_lo = recommend_ta(tm_perfect(f_seq, lo_c), tm_perfect(r_seq, lo_c))
    ta_hi = recommend_ta(tm_perfect(f_seq, hi_c), tm_perfect(r_seq, hi_c))
    st.info(f"**Mg²⁺ 不確定性** — Mg²⁺ 取 {mg_lo:.1f}–{mg_hi:.1f} mM 時，"
            f"建議 Ta 落在 **{ta_lo:.1f} – {ta_hi:.1f} °C**"
            f"（寬度僅 {ta_hi - ta_lo:.1f}°C）。"
            f"相較之下，若誤設 Mg²⁺ = 0 會低估約 10°C——"
            f"不知道 Mg²⁺ 是小問題，假設它是 0 才是大問題。")

st.caption(f"條件：{cond.label()}")

with st.expander("📉 Mg²⁺ 敏感度表"):
    rows = []
    for m in (0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0):
        cm = replace(cond, mg=m)
        a, b = tm_perfect(f_seq, cm), tm_perfect(r_seq, cm)
        rows.append({"Mg²⁺ (mM)": m, "Tm-F (°C)": round(a, 1),
                     "Tm-R (°C)": round(b, 1),
                     "建議 Ta (°C)": round(recommend_ta(a, b), 1)})
    st.dataframe(pd.DataFrame(rows), hide_index=True)
    st.caption("Mg²⁺ 的效應會飽和：0.5→1.5 掉 3.5°C，但 2.5→4.0 只差 0.8°C。")

with st.expander("🔧 Mg²⁺ 反推（用一組實測 Ta 校準整個模型）", expanded=mg_mode.startswith("③")):
    st.write("填入一組你**實際跑過、能得到乾淨單一 band** 的 Ta，"
             "系統會反解出符合該結果的有效 Mg²⁺。")
    sc1, sc2 = st.columns([1, 3])
    obs_ta = sc1.number_input("實測可行 Ta (°C)", 35.0, 75.0, 55.0, 0.5)
    if sc2.button("反推有效 Mg²⁺"):
        mg_sol, msg = solve_mg(obs_ta, f_seq, r_seq, cond)
        if mg_sol is None:
            st.warning(msg)
        else:
            st.session_state["mg_solved"] = round(mg_sol, 2)
            st.success(msg + "　已存入，請把左側切到模式 ③ 套用。")

st.divider()


# ---------- BLAST ----------
st.subheader("🔍 引子結合位點（NCBI BLAST）")

sig = (f_seq, r_seq, organism, database)
run = st.button("🚀 查詢 NCBI", type="primary")

if run:
    st.session_state["blast_sig"] = sig

if st.session_state.get("blast_sig") != sig:
    st.info("按上方按鈕向 NCBI 查詢。第一次約需 1–3 分鐘（NCBI 佇列），"
            "結果會快取 24 小時——之後調整溫度或反應條件都是瞬間反應，不會重查。")
    st.stop()

try:
    with st.status("向 NCBI 查詢中…", expanded=True) as status:
        st.write("⏳ Forward 引子…")
        xml_f = blast_xml(f_seq, organism, database)
        st.write("⏳ Reverse 引子…")
        xml_r = blast_xml(r_seq, organism, database)
        st.write("🧮 解析比對結果並計算各位點 Tm…")
        hits_f = parse_hits(xml_f, f_seq, "F", cond.key())
        hits_r = parse_hits(xml_r, r_seq, "R", cond.key())
        status.update(label=f"完成：F {len(hits_f)} 個 / R {len(hits_r)} 個結合位點",
                      state="complete", expanded=False)
except Exception as e:
    st.error(f"BLAST 失敗：{type(e).__name__}: {e}\n\n"
             "NCBI 偶爾會因流量限制拒絕連線，稍候再試。")
    st.stop()


# ---------- 溫度模擬 ----------
st.subheader("🎛️ 溫度模擬")

sl1, sl2 = st.columns([3, 1])
temp = sl1.slider("Annealing 溫度 Ta (°C)", 35.0, 78.0,
                  float(np.clip(round(ta * 2) / 2, 35.0, 78.0)), 0.5)
allow_self = sl2.checkbox("計入 F-F / R-R 產物", value=True,
                          help="同一條引子分別落在正負股時也會產生真實的非專一性條帶。")

products = pair_products(hits_f, hits_r, temp, allow_self)

left, right = st.columns([3, 2])

with left:
    if not products:
        st.warning("❄️ 此溫度下無任何產物 —— 引子已脫落，或 3' 端無法延伸。試著調低溫度。")
    else:
        n_t = sum(1 for p in products if p["perfect"])
        a, b = st.columns(2)
        a.metric("完美配對產物", n_t)
        b.metric("非專一性產物", len(products) - n_t)

        df = pd.DataFrame([{
            "大小 (bp)": p["size"],
            "類型": "✅ Target" if p["perfect"] else f"⚠️ 非專一 {p['pair']}",
            "Tm-F": round(p["f_tm"], 1),
            "Tm-R": round(p["r_tm"], 1),
            "餘裕 (°C)": round(p["margin"], 1),
            "Accession": p["id"],
            "描述": p["title"][:70],
        } for p in products])
        st.dataframe(df, hide_index=True, height=340)

with right:
    if products:
        st.pyplot(draw_gel(products, temp, organism))

with st.expander("🔬 所有結合位點明細"):
    hd = pd.DataFrame([{
        "引子": h["tag"], "Accession": h["id"],
        "股": "+" if h["strand"] == 1 else "−",
        "5'位置": h["five_prime"], "3'位置": h["three_prime"],
        "Identity (%)": round(h["ident"], 1),
        "Coverage (%)": round(h["coverage"], 1),
        "Tm (°C)": round(h["tm"], 1) if h["tm"] is not None else None,
        "可延伸": "✔" if h["ext_ok"] else "✘",
        "3'端狀態": h["ext_note"],
    } for h in hits_f + hits_r])
    st.dataframe(hd.sort_values(["引子", "Tm (°C)"], ascending=[True, False]),
                 hide_index=True, height=380)
    st.caption("Identity 的分母是引子全長（不是 HSP 長度），所以只對到一半的短 HSP "
               "不會再被誤判成 100% 的完美結合。")
