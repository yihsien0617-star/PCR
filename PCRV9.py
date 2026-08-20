# @title 🧬 PCR 戰情室 (v9.2 — Mg²⁺ 未知模式)
# =============================================================================
#  v9.1 -> v9.2
#  商業 master mix 幾乎都不公開 Mg²⁺ 濃度，所以 Mg²⁺ 改為三種模式：
#    A. 已知濃度      —— 直接輸入（自行配製 buffer 時用）
#    B. 未知：區間估計 —— 給一個合理範圍，Tm / Ta 以「區間」呈現
#    C. 未知：反推     —— 輸入一組實驗上跑得乾淨的 Ta，反解有效 Mg²⁺
#  另附 Mg²⁺ 敏感度表，讓你一眼看出這個未知數到底影響多大。
#
#  v9.0 -> v9.1 的修正仍全部保留（Tm 熱力學、identity 分母、amplicon 座標、
#  3' 端朝向與延伸規則、F-F/R-R 配對、widget 重複綁定…）
# =============================================================================

%matplotlib inline
import ipywidgets as widgets
from IPython.display import display, clear_output
import matplotlib.pyplot as plt
import numpy as np
import threading

try:
    from Bio.Blast import NCBIWWW, NCBIXML
    from Bio.SeqUtils import MeltingTemp as mt
    from Bio.Seq import Seq
except ImportError:
    raise SystemExit("⚠️ 請先安裝 Biopython：pip install biopython")

# ================= 全域狀態 =================
raw_blast_data = {"f": [], "r": []}
primer_info = {"tm_f": None, "tm_r": None, "seq_f": "", "seq_r": "", "ta": None}
_busy = threading.Lock()


# =============================================================================
#  1. Mg²⁺ 模式
# =============================================================================

def mg_bounds():
    """回傳 (下限, 代表值, 上限)。三種模式共用的單一入口。"""
    mode = w_mg_mode.value
    if mode == "known":
        v = float(w_mg.value)
        return v, v, v
    if mode == "solve":
        v = float(w_mg_solved.value)
        return v, v, v
    lo, hi = float(w_mg_lo.value), float(w_mg_hi.value)
    if lo > hi:
        lo, hi = hi, lo
    return lo, (lo + hi) / 2.0, hi


def effective_mg():
    """點估計時使用的 Mg²⁺（區間模式取中點）。"""
    return mg_bounds()[1]


def nn_kwargs(mg=None):
    """Bio.SeqUtils.MeltingTemp 參數。

    saltcorr=7 (Owczarzy 2008) 是 Biopython 中唯一同時正確處理
    Mg²⁺ 與 dNTP 螯合的鹽校正法。注意：當 Mg <= dNTP 時公式會退回
    只算一價鹽，所以 Mg 在 0 ~ dNTP 之間 Tm 不會變動。
    """
    return dict(
        dnac1=float(w_primer_nM.value),
        dnac2=0.0,
        selfcomp=False,
        Na=0.0,
        K=float(w_mono.value),
        Tris=float(w_tris.value),
        Mg=float(effective_mg() if mg is None else mg),
        dNTPs=float(w_dntp.value),
        nn_table=mt.DNA_NN4,
        tmm_table=mt.DNA_TMM1,
        imm_table=mt.DNA_IMM1,
        de_table=mt.DNA_DE1,
        saltcorr=7,
        strict=False,
    )


# =============================================================================
#  2. 熱力學核心
# =============================================================================

def tm_perfect(seq, mg=None):
    seq = seq.strip().upper().replace(" ", "")
    if not seq:
        return None
    return mt.Tm_NN(Seq(seq), **nn_kwargs(mg)) + float(w_offset.value)


def tm_span(seq):
    """回傳 (下限, 代表值, 上限) 的 Tm。非區間模式時三者相同。"""
    lo, mid, hi = mg_bounds()
    return (tm_perfect(seq, lo), tm_perfect(seq, mid), tm_perfect(seq, hi))


def tm_duplex(query_aln, sbjct_aln, mg=None):
    """由 BLAST 比對字串直接算含錯配的真實雙股 Tm。

    BLAST 的 hsp.query / hsp.sbjct 是「同義股」對齊（字母相同 = match），
    引子實際結合的是 sbjct 的『互補』股，左到右恰為 3'->5'，
    正是 Biopython c_seq 要求的方向。
    """
    q, s = query_aln.upper(), sbjct_aln.upper()
    if len(q) != len(s):
        return None, "length-mismatch"

    if "-" in q or "-" in s:
        # 含 indel：Tm_NN 不支援，退回保守經驗式
        ident = sum(1 for a, b in zip(q, s) if a == b and a != "-")
        base = tm_perfect(q.replace("-", ""), mg)
        if base is None:
            return None, "gap-fallback"
        n = max(len(q.replace("-", "")), 1)
        return base - (1 - ident / n) * 100 * 1.2, "gap-fallback"

    try:
        c = str(Seq(s).complement())
        return mt.Tm_NN(Seq(q), c_seq=Seq(c), **nn_kwargs(mg)) \
               + float(w_offset.value), "NN-mismatch"
    except Exception:
        return None, "nn-failed"


def recommend_ta(tm_f, tm_r):
    """Taq / hot-start Taq 慣例：Ta = 較低的 Tm − 5，上限 72°C。"""
    if tm_f is None or tm_r is None:
        return None
    return min(min(tm_f, tm_r) - 5.0, 72.0)


def solve_mg(target_ta, f_seq, r_seq):
    """由「實驗上跑得乾淨的 Ta」反解有效 Mg²⁺（二分法）。

    Tm 對 Mg²⁺ 在 Mg > dNTP 的區間內單調遞增，所以二分法收斂。
    搜尋下界刻意設在 dNTP+0.3 以避開 Owczarzy 公式的分支切換點。
    回傳 (mg, 說明字串)；無解時 mg 為 None。
    """
    target_tm = target_ta + 5.0          # 反推 recommend_ta

    def limiting_tm(mg):
        a, b = tm_perfect(f_seq, mg), tm_perfect(r_seq, mg)
        return min(a, b)

    lo = float(w_dntp.value) + 0.3
    hi = 10.0
    f_lo, f_hi = limiting_tm(lo), limiting_tm(hi)

    if target_tm < f_lo:
        return None, (f"實驗 Ta 偏低：即使 Mg²⁺ 只有 {lo:.1f} mM，"
                      f"限速引子 Tm 仍有 {f_lo:.1f}°C（對應 Ta {f_lo-5:.1f}°C）。"
                      f"差額請用 offset 補（建議 offset = {target_tm - f_lo:+.1f}°C）。")
    if target_tm > f_hi:
        return None, (f"實驗 Ta 偏高：即使 Mg²⁺ 拉到 {hi:.0f} mM，"
                      f"限速引子 Tm 只有 {f_hi:.1f}°C（對應 Ta {f_hi-5:.1f}°C）。"
                      f"差額請用 offset 補（建議 offset = {target_tm - f_hi:+.1f}°C）。")

    for _ in range(60):
        mid = (lo + hi) / 2.0
        if limiting_tm(mid) < target_tm:
            lo = mid
        else:
            hi = mid
    mg = (lo + hi) / 2.0
    return mg, (f"有效 Mg²⁺ ≈ {mg:.2f} mM"
                + ("　（合理範圍內 ✔）" if 0.8 <= mg <= 4.0 else
                   "　⚠️ 超出一般 master mix 範圍，可能是其他因素（酶種類、"
                   "增強劑、引子純度）造成，建議改用 offset 校正"))


# =============================================================================
#  3. 3' 端規則
# =============================================================================

THREE_PRIME_CRITICAL = 2   # 末 2 nt 錯配 -> 幾乎完全不延伸
THREE_PRIME_WINDOW = 5     # 末 5 nt 內的錯配額外註記


def three_prime_status(query_aln, sbjct_aln, q_end, plen):
    """q_end < plen 代表 HSP 沒蓋到引子 3' 端 —— BLAST 會在錯配處截斷比對，
    所以這通常就是「3' 端錯配」，必須視為不可延伸。"""
    if q_end < plen:
        return False, f"3'端未比對到 (缺 {plen - q_end} nt)"

    q, s = query_aln.upper(), sbjct_aln.upper()
    tq, ts = q[-THREE_PRIME_WINDOW:], s[-THREE_PRIME_WINDOW:]
    mism = [i for i, (a, b) in enumerate(zip(tq, ts)) if a != b]
    if not mism:
        return True, "3'端完全配對"
    dist = [len(tq) - 1 - i for i in mism]     # 0 = 最末鹼基
    if min(dist) < THREE_PRIME_CRITICAL:
        return False, f"3'末端第 {min(dist)+1} 位錯配 → 無法延伸"
    return True, f"3'端 -{min(dist)+1} 位錯配 (效率下降)"


# =============================================================================
#  4. BLAST 解析
# =============================================================================

def parse_hits(record, primer_seq, tag):
    plen = len(primer_seq)
    hits = []
    for a in record.alignments:
        for h in a.hsps:
            s_start, s_end = h.sbjct_start, h.sbjct_end
            strand = 1 if s_end >= s_start else -1      # 負股 hit-from > hit-to

            ident_full = h.identities / plen * 100.0    # 分母必須是引子全長
            coverage = h.align_length / plen * 100.0

            # 把 HSP 外推到引子兩端（比對欄位左到右 = query 5'->3'）
            five_prime = s_start - strand * (h.query_start - 1)
            three_prime = s_end + strand * (plen - h.query_end)

            tm, method = tm_duplex(h.query, h.sbjct)
            ext_ok, ext_note = three_prime_status(h.query, h.sbjct,
                                                  h.query_end, plen)
            hits.append({
                "tag": tag, "id": a.accession, "title": a.title,
                "aln": (h.query, h.sbjct),
                "strand": strand,
                "five_prime": five_prime, "three_prime": three_prime,
                "ident": ident_full, "coverage": coverage,
                "tm": tm, "tm_method": method,
                "ext_ok": ext_ok, "ext_note": ext_note,
                "perfect": (ident_full >= 99.99 and coverage >= 99.99),
            })
    return hits


# =============================================================================
#  5. 產物模擬
# =============================================================================

MIN_PRODUCT, MAX_PRODUCT = 50, 5000


def pair_products(hits_a, hits_b, temp, allow_self=True):
    """擴增條件：同一條序列、一正一負股、兩條 3' 端『面對面』、
    且在該溫度下都還結合著且 3' 端可延伸。
    amplicon 長度 = 兩條引子 5' 端之間的距離。"""
    def alive(h):
        return h["ext_ok"] and h["tm"] is not None and h["tm"] >= temp

    live = [h for h in hits_a if alive(h)] + [h for h in hits_b if alive(h)]

    by_acc = {}
    for h in live:
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
                    continue                      # 3' 端背對背 -> 不擴增
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
                    "f_ident": plus["ident"], "r_ident": minus["ident"],
                    "f_tm": plus["tm"], "r_tm": minus["tm"],
                    "perfect": plus["perfect"] and minus["perfect"]
                               and plus["tag"] != minus["tag"],
                    "margin": min(plus["tm"], minus["tm"]) - temp,
                })
    products.sort(key=lambda p: (-p["perfect"], p["size"]))
    return products


# =============================================================================
#  6. 繪圖
# =============================================================================

def draw_gel(products, temp, organism):
    plt.close("all")
    fig, ax = plt.subplots(figsize=(9, 5.5), dpi=100)
    fig.patch.set_facecolor("#111111")
    ax.set_facecolor("black")

    for b in [100, 200, 300, 400, 500, 600, 800, 1000, 1500, 2000, 3000, 5000]:
        y = -np.log10(b)
        ax.hlines(y, 0.5, 1.5, colors="white", alpha=0.45, linewidth=1.5)
        ax.text(0.35, y, f"{b}", color="white", fontsize=7, va="center", ha="right")

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
            color, lw, tag = "#FFD400", 2.0, f"Non-specific {p['pair']}"
        ax.hlines(y, 2.5, 3.6, colors=color, alpha=alpha, linewidth=lw)
        if labelled < 6:
            ax.text(3.75, y, f"{p['size']} bp  ({tag})", color=color,
                    fontsize=8.5, va="center", alpha=max(alpha, 0.7))
            labelled += 1

    ax.set_xlim(0, 6.6)
    ax.set_ylim(-np.log10(MAX_PRODUCT * 1.15), -np.log10(MIN_PRODUCT * 0.85))
    ax.set_xticks([1, 3.05])
    ax.set_xticklabels(["Marker", f"PCR @ {temp:.1f}°C"], color="white")
    ax.tick_params(colors="white")
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_color("#444444")
    ax.set_title(f"In-Silico PCR  ({organism})", color="white", fontsize=11)
    fig.tight_layout()
    display(fig)
    plt.close(fig)


# =============================================================================
#  7. 模擬 / Tm 面板
# =============================================================================

def gc_pct(s):
    s = s.upper()
    return 100.0 * (s.count("G") + s.count("C")) / max(len(s), 1)


def _cond_line():
    lo, mid, hi = mg_bounds()
    mg_txt = f"{mid:.2f} mM" if lo == hi else f"{lo:.1f}–{hi:.1f} mM (取中點 {mid:.2f})"
    return (f"{w_mono.value:g} mM K⁺ / {w_tris.value:g} mM Tris / Mg²⁺ {mg_txt} / "
            f"{w_dntp.value:g} mM dNTP / {w_primer_nM.value:g} nM primer"
            f"  (offset {w_offset.value:+.1f}°C)")


def update_simulation(temp=None):
    if temp is None:
        temp = w_slider.value
    temp = float(temp)

    with sim_output:
        clear_output(wait=True)
        f_hits, r_hits = raw_blast_data["f"], raw_blast_data["r"]
        if not f_hits and not r_hits:
            print("尚無 BLAST 資料，請先按上方按鈕取得引子結合位點。")
            return

        products = pair_products(f_hits, r_hits, temp, allow_self=w_selfpair.value)
        tm_f, tm_r, ta = primer_info["tm_f"], primer_info["tm_r"], primer_info["ta"]

        print(f"🌡️  設定 Ta = {temp:.1f}°C" + (f"   (建議值 {ta:.1f}°C)" if ta else ""))
        if tm_f and tm_r:
            print(f"    完美結合 Tm：F = {tm_f:.1f}°C，R = {tm_r:.1f}°C")
        print(f"    條件：{_cond_line()}")
        print("-" * 72)

        if not products:
            print("❄️  此溫度下無任何產物 —— 引子已脫落或 3' 端無法延伸。")
            return

        n_t = sum(1 for p in products if p["perfect"])
        print(f"📊 共 {len(products)} 條產物"
              f"（{n_t} 條完美配對，{len(products)-n_t} 條非專一性）")
        for p in products[:15]:
            mark = "✅ Target" if p["perfect"] else f"⚠️  非專一 {p['pair']}"
            print(f"   {p['size']:>5} bp  {mark}"
                  f"   Tm(F/R)={p['f_tm']:.1f}/{p['r_tm']:.1f}"
                  f"   餘裕 +{p['margin']:.1f}°C")
            print(f"          {p['id']}  {p['title'][:58]}")
        if len(products) > 15:
            print(f"   ... 另有 {len(products)-15} 條未列出")

        draw_gel(products, temp, w_org.value)


def refresh_tm_panel(*_):
    """引子或反應條件改變時重算 Tm。不需要重跑 BLAST。"""
    f, r = w_fwd.value.strip().upper(), w_rev.value.strip().upper()
    with tm_output:
        clear_output(wait=True)
        if not f or not r:
            print("填入兩條引子後，這裡會即時顯示 Tm 與建議 Ta。")
            return
        try:
            f_lo, f_mid, f_hi = tm_span(f)
            r_lo, r_mid, r_hi = tm_span(r)
        except Exception as e:
            print(f"⚠️ 序列無法計算 Tm（只接受 A/T/G/C）：{e}")
            return

        ta = recommend_ta(f_mid, r_mid)
        primer_info.update(tm_f=f_mid, tm_r=r_mid, seq_f=f, seq_r=r, ta=ta)
        ranged = (f_lo != f_hi)

        def fmt(lo, mid, hi):
            return f"{mid:5.1f}°C" + (f"  [{lo:.1f} – {hi:.1f}]" if ranged else "")

        print(f"Forward  {len(f):2d} nt  GC {gc_pct(f):4.1f}%   Tm = {fmt(f_lo,f_mid,f_hi)}")
        print(f"Reverse  {len(r):2d} nt  GC {gc_pct(r):4.1f}%   Tm = {fmt(r_lo,r_mid,r_hi)}")
        dtm = abs(f_mid - r_mid)
        print(f"ΔTm = {dtm:.1f}°C" +
              ("   ⚠️ 相差 >5°C，建議重新設計" if dtm > 5 else "   ✔"))

        if ranged:
            ta_lo = recommend_ta(f_lo, r_lo)
            ta_hi = recommend_ta(f_hi, r_hi)
            print(f"👉 建議 Ta：{ta:.1f}°C　【Mg²⁺ 不確定性造成的範圍 "
                  f"{ta_lo:.1f} – {ta_hi:.1f}°C，寬度 {ta_hi-ta_lo:.1f}°C】")
        else:
            print(f"👉 建議 Ta（Taq 慣例 Tm−5）：{ta:.1f}°C")

        # ---- Mg²⁺ 敏感度表：讓「未知」這件事變得可量化 ----
        print("\n   Mg²⁺ 敏感度（其餘條件不變）")
        print("   Mg(mM)   Tm-F    Tm-R   建議Ta")
        for m in (0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0):
            a, b = tm_perfect(f, m), tm_perfect(r, m)
            print(f"   {m:5.1f}  {a:6.1f}  {b:6.1f}   {recommend_ta(a,b):5.1f}"
                  + ("   ← 目前" if abs(m - effective_mg()) < 0.26 else ""))

    if raw_blast_data["f"] or raw_blast_data["r"]:
        for key in ("f", "r"):
            for h in raw_blast_data[key]:
                h["tm"], h["tm_method"] = tm_duplex(*h["aln"])
        update_simulation()


def on_solve(_):
    f, r = w_fwd.value.strip().upper(), w_rev.value.strip().upper()
    with solve_output:
        clear_output(wait=True)
        if not f or not r:
            print("請先填入兩條引子序列。")
            return
        mg, msg = solve_mg(float(w_obs_ta.value), f, r)
        print(f"實驗 Ta = {w_obs_ta.value:.1f}°C  →  {msg}")
        if mg is not None:
            w_mg_solved.value = round(mg, 2)
            print("已套用。之後所有引子都會沿用這個有效 Mg²⁺。")
    refresh_tm_panel()


def on_mg_mode(_):
    m = w_mg_mode.value
    box_known.layout.display = "" if m == "known" else "none"
    box_range.layout.display = "" if m == "range" else "none"
    box_solve.layout.display = "" if m == "solve" else "none"
    refresh_tm_panel()


# =============================================================================
#  8. BLAST 執行緒
# =============================================================================

def run_blast_thread(f_seq, r_seq, organism, database):
    if not _busy.acquire(blocking=False):
        return
    try:
        with log_output:
            clear_output(wait=True)
            try:
                w_progress.value, w_progress.bar_style = 5, "info"
                f_seq, r_seq = f_seq.strip().upper(), r_seq.strip().upper()

                tm_f, tm_r = tm_perfect(f_seq), tm_perfect(r_seq)
                ta = recommend_ta(tm_f, tm_r)
                primer_info.update(tm_f=tm_f, tm_r=tm_r,
                                   seq_f=f_seq, seq_r=r_seq, ta=ta)
                print(f"🧮 Tm  F={tm_f:.1f}°C  R={tm_r:.1f}°C   建議 Ta={ta:.1f}°C")
                print(f"   條件：{_cond_line()}")

                bkw = dict(entrez_query=f'"{organism}"[Organism]',
                           expect=1000, word_size=7,
                           hitlist_size=500, megablast=False)

                w_progress.value = 20
                print("⏳ Forward 引子 BLAST 中…（NCBI 佇列可能要 1–3 分鐘）")
                rec_f = NCBIXML.read(NCBIWWW.qblast("blastn", database, f_seq, **bkw))

                w_progress.value = 55
                print("⏳ Reverse 引子 BLAST 中…")
                rec_r = NCBIXML.read(NCBIWWW.qblast("blastn", database, r_seq, **bkw))

                w_progress.value = 85
                raw_blast_data["f"] = parse_hits(rec_f, f_seq, "F")
                raw_blast_data["r"] = parse_hits(rec_r, r_seq, "R")
                print(f"✅ 取得 F {len(raw_blast_data['f'])} 個 / "
                      f"R {len(raw_blast_data['r'])} 個結合位點")

                w_progress.value, w_progress.bar_style = 100, "success"

                lo, hi = w_slider.min, w_slider.max
                w_slider.value = float(np.clip(round(ta * 2) / 2, lo, hi))
                if not (lo <= ta <= hi):
                    print(f"⚠️ 建議 Ta {ta:.1f}°C 超出滑桿範圍 [{lo}, {hi}]，已夾擠")
                update_simulation()

            except Exception as e:
                w_progress.bar_style = "danger"
                print(f"❌ 錯誤：{type(e).__name__}: {e}")
    finally:
        btn_run.disabled = False
        _busy.release()


def on_run(_):
    if not w_fwd.value.strip() or not w_rev.value.strip():
        with log_output:
            clear_output(wait=True)
            print("請先填入兩條引子序列。")
        return
    btn_run.disabled = True
    w_progress.value, w_progress.bar_style = 0, "info"
    sim_output.clear_output()
    threading.Thread(target=run_blast_thread,
                     args=(w_fwd.value, w_rev.value, w_org.value, w_db.value),
                     daemon=True).start()


# =============================================================================
#  9. 介面
# =============================================================================

style = {"description_width": "initial"}
wide = widgets.Layout(width="98%")
half = widgets.Layout(width="48%")
num = widgets.Layout(width="215px")

w_org = widgets.Dropdown(
    options=[("人類", "Homo sapiens"), ("小鼠", "Mus musculus")],
    value="Homo sapiens", description="物種:", style=style, layout=half)
w_db = widgets.Dropdown(
    options=[("mRNA (refseq_rna)", "refseq_rna"), ("Genomic (nt)", "nt")],
    value="refseq_rna", description="資料庫:", style=style, layout=half)
w_fwd = widgets.Text(placeholder="Forward primer 5'->3'", continuous_update=False,
                     description="Fwd:", style=style, layout=wide)
w_rev = widgets.Text(placeholder="Reverse primer 5'->3'", continuous_update=False,
                     description="Rev:", style=style, layout=wide)

# ---- Mg²⁺ 三模式 ----
w_mg_mode = widgets.RadioButtons(
    options=[("① 已知濃度（自配 buffer）", "known"),
             ("② 未知 — 用區間估計（推薦：商業 master mix）", "range"),
             ("③ 未知 — 由實驗 Ta 反推有效濃度（最準）", "solve")],
    value="range", description="Mg²⁺:", style=style,
    layout=widgets.Layout(width="98%"))

w_mg = widgets.BoundedFloatText(value=1.5, min=0, max=10, step=0.1,
                                description="Mg²⁺ (mM):", style=style, layout=num)
box_known = widgets.HBox([w_mg])

w_mg_lo = widgets.BoundedFloatText(value=1.0, min=0, max=10, step=0.1,
                                   description="下限 (mM):", style=style, layout=num)
w_mg_hi = widgets.BoundedFloatText(value=3.0, min=0, max=10, step=0.1,
                                   description="上限 (mM):", style=style, layout=num)
box_range = widgets.VBox([
    widgets.HBox([w_mg_lo, w_mg_hi]),
    widgets.HTML("<small>市售 2X Taq master mix 稀釋成 1X 後，Mg²⁺ 幾乎都落在 "
                 "1.5–2.5 mM；1.0–3.0 是保守的涵蓋範圍。Tm 與建議 Ta 會以區間顯示。"
                 "</small>")])

w_obs_ta = widgets.BoundedFloatText(value=55.0, min=35, max=75, step=0.5,
                                    description="實驗可行 Ta (°C):", style=style,
                                    layout=widgets.Layout(width="260px"))
btn_solve = widgets.Button(description="反推有效 Mg²⁺", button_style="info",
                           icon="calculator", layout=widgets.Layout(width="180px"))
w_mg_solved = widgets.BoundedFloatText(value=1.5, min=0, max=10, step=0.01,
                                       description="反推結果 (mM):", style=style,
                                       disabled=True, layout=num)
solve_output = widgets.Output()
box_solve = widgets.VBox([
    widgets.HTML("<small>填入一組你<b>實際跑過、能得到乾淨單一 band</b> 的引子與 Ta，"
                 "系統會反解出符合該結果的有效 Mg²⁺。一組已知答案就能校準整個模型。"
                 "</small>"),
    widgets.HBox([w_obs_ta, btn_solve, w_mg_solved]),
    solve_output])

w_dntp = widgets.BoundedFloatText(value=0.2, min=0, max=2, step=0.05,
                                  description="dNTP each (mM):", style=style, layout=num)
w_mono = widgets.BoundedFloatText(value=50, min=0, max=200, step=5,
                                  description="K⁺ (mM):", style=style, layout=num)
w_tris = widgets.BoundedFloatText(value=10, min=0, max=100, step=1,
                                  description="Tris (mM):", style=style, layout=num)
w_primer_nM = widgets.BoundedFloatText(value=200, min=10, max=2000, step=10,
                                       description="Primer (nM):", style=style, layout=num)
w_offset = widgets.BoundedFloatText(value=0.0, min=-15, max=15, step=0.5,
                                    description="校正 offset (°C):", style=style, layout=num)
w_selfpair = widgets.Checkbox(value=True, description="計入 F-F / R-R 自我配對產物",
                              style=style, indent=False)

tm_output = widgets.Output()
btn_run = widgets.Button(description="連線 NCBI 取得結合位點", button_style="primary",
                         icon="cloud-download", layout=widgets.Layout(width="100%"))
w_progress = widgets.IntProgress(value=0, max=100, description="進度:",
                                 bar_style="info", layout=wide)
log_output = widgets.Output()
w_slider = widgets.FloatSlider(value=58.0, min=35.0, max=78.0, step=0.5,
                               description="Annealing Ta (°C):",
                               continuous_update=False, readout_format=".1f",
                               style=style, layout=wide)
sim_output = widgets.Output()

# ---- 事件綁定（每個只綁一次）----
btn_run.on_click(on_run)
btn_solve.on_click(on_solve)
w_slider.observe(lambda ch: update_simulation(ch["new"]), names="value")
w_mg_mode.observe(on_mg_mode, names="value")
w_selfpair.observe(lambda ch: update_simulation(), names="value")
for _w in (w_mg, w_mg_lo, w_mg_hi, w_mg_solved, w_dntp, w_mono,
           w_tris, w_primer_nM, w_offset, w_fwd, w_rev):
    _w.observe(refresh_tm_panel, names="value")

ui = widgets.VBox([
    widgets.HTML("<h3>🧬 PCR 戰情室 <small>v9.2 — Mg²⁺ 未知模式</small></h3>"),
    widgets.HBox([w_org, w_db]),
    w_fwd, w_rev,
    widgets.HTML("<hr><b>🧪 反應條件</b>"),
    w_mg_mode, box_known, box_range, box_solve,
    widgets.HTML("<br>"),
    widgets.HBox([w_dntp, w_mono, w_tris]),
    widgets.HBox([w_primer_nM, w_offset]),
    tm_output,
    widgets.HTML("<hr>"),
    btn_run, w_progress, log_output,
    widgets.HTML("<hr><b>🎛️ 溫度模擬</b>"),
    w_slider, w_selfpair, sim_output,
])

on_mg_mode(None)     # 初始化顯示/隱藏
display(ui)
