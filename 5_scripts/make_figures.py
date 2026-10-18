"""Figures 1-4 for the J Phycol manuscript (v4).
Run with a python that has matplotlib + openpyxl:  python make_figures.py
Inputs: algal_telomere_study/data/processed/{figure_data,tr_search/rnaseq}, Table_S1_assemblies.xlsx.
Outputs: figures/Figure_<n>.png (600 dpi) and .pdf"""
import os, re, csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch

D = os.path.dirname(os.path.abspath(__file__))
P = os.path.join(D, "../../algal_telomere_study/data/processed")
OUT = os.path.join(D, "figures"); os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.family": "Arial", "font.size": 8, "pdf.fonttype": 42, "axes.linewidth": 0.6})

# one colour per repeat, used in every figure (colour-blind-safe Okabe-Ito)
COL = {2: "#E69F00", 3: "#0072B2", 4: "#999999", 5: "#009E73", 6: "#CC79A7"}
REP = {2: "TTAGGG", 3: "TTTAGGG", 4: "TTTTAGGG", 5: "TTTTTAGGG", 6: "TTTTTTAGGG"}
KEYLAB = {2: "T2 = TTAGGG", 3: "T3 = TTTAGGG", 4: "T4 = TTTTAGGG", 5: "T5 = TTTTTAGGG", 6: "T6 = TTTTTTAGGG"}


def colour_key(target, items, ncol, title="Color key", **kw):
    """Draw a colour key (list of (colour, label)) as a legend on a figure or axes."""
    h = [Patch(facecolor=c, edgecolor="none", label=l) for c, l in items]
    return target.legend(handles=h, ncol=ncol, frameon=False, fontsize=6.5, title=title, title_fontsize=6.8,
                         handlelength=1.4, handleheight=0.9, columnspacing=1.2, **kw)


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"), dpi=600, bbox_inches="tight")
    plt.close(fig)
    print("wrote", name)


# ---------------- Figure 1: five-order tree with repeats ----------------
def figure1():
    fig, ax = plt.subplots(figsize=(8.8, 3.6))
    ax.set_xlim(0, 13.4); ax.set_ylim(0.3, 7.0); ax.axis("off")
    # leaves (y) and topology after Hess et al. (2022), Fig. 2
    leaf = {"Spirogloeales": 6.3, "Serritaeniales": 5.3, "Zygnematales_a": 4.4, "Zygnematales_b": 3.7,
            "Spirogyrales": 2.5, "Desmidiales": 1.3}
    xl = 3.2
    lw = 1.1
    def h(x0, x1, y, c="k"): ax.plot([x0, x1], [y, y], c=c, lw=lw, solid_capstyle="butt")
    def v(x, y0, y1): ax.plot([x, x], [y0, y1], c="k", lw=lw)
    # (Spirogyrales, Desmidiales)
    xa = 2.4; v(xa, leaf["Desmidiales"], leaf["Spirogyrales"])
    h(xa, xl, leaf["Spirogyrales"], COL[2]); h(xa, xl, leaf["Desmidiales"])
    # Zygnema pair
    xz = 2.7; v(xz, leaf["Zygnematales_b"], leaf["Zygnematales_a"])
    h(xz, xl, leaf["Zygnematales_a"]); h(xz, xl, leaf["Zygnematales_b"], COL[6])
    yz = (leaf["Zygnematales_a"] + leaf["Zygnematales_b"]) / 2
    # Zygnematales + (Spirogyrales, Desmidiales)
    ya = (leaf["Spirogyrales"] + leaf["Desmidiales"]) / 2
    xb = 1.7; v(xb, ya, yz); h(xb, xa, ya); h(xb, xz, yz, COL[5])
    yb = (ya + yz) / 2
    # + Serritaeniales
    xc = 1.1; v(xc, yb, leaf["Serritaeniales"]); h(xc, xb, yb); h(xc, xl, leaf["Serritaeniales"])
    yc = (yb + leaf["Serritaeniales"]) / 2
    # + Spirogloeales
    xd = 0.5; v(xd, yc, leaf["Spirogloeales"]); h(xd, xc, yc); h(xd, xl, leaf["Spirogloeales"])
    yd = (yc + leaf["Spirogloeales"]) / 2
    h(0.15, xd, yd)
    # change markers
    ax.text((xb + xz) / 2, yz + 0.12, "T3→T5", color=COL[5], ha="center", fontsize=7, weight="bold")
    ax.text((xz + xl) / 2, leaf["Zygnematales_b"] - 0.32, "T5→T6?", color=COL[6], ha="center", fontsize=7, weight="bold")
    ax.text((xa + xl) / 2, leaf["Spirogyrales"] + 0.12, "T3→T2", color=COL[2], ha="center", fontsize=7, weight="bold")
    ax.text(0.12, yd + 0.08, "T3", color=COL[3], fontsize=7, weight="bold", va="bottom")
    # rows
    rows = [("Spirogloeales", "Spirogloeales", "Spirogloea muscicola", 3, "T3 (published TR)", "complete"),
            ("Serritaeniales", "Serritaeniales", "“Mesotaenium endlicherianum”", 3, "T3 (published TR)", "complete"),
            ("Zygnematales_a", "Zygnematales", "Z. circumcarinatum (3 strains)", 5, "T4; expressed", "complete"),
            ("Zygnematales_b", "", "Z. cf. cylindricum", 6, "T4", "split (2 scaffolds)"),
            ("Spirogyrales", "Spirogyrales", "Spirogyra pratensis", 2, "not found", "partial"),
            ("Desmidiales", "Desmidiales", "Closterium sp. (8 strains)", 3, "not found", "complete")]
    xs = {"order": 3.3, "taxon": 4.75, "rep": 8.0, "tr": 9.75, "tert": 11.6}
    hy = 6.85
    for k, t in [("order", "Order"), ("taxon", "Taxon"), ("rep", "Telomeric repeat"), ("tr", "TR template"),
                 ("tert", "TERT")]:
        ax.text(xs[k], hy, t, weight="bold", fontsize=7.5, va="center")
    for key, order, taxon, t, tr, tert in rows:
        y = leaf[key]
        ax.text(xs["order"], y, order, va="center", fontsize=7.5)
        ax.text(xs["taxon"], y, taxon, va="center", fontsize=7.5, style="italic")
        ax.add_patch(Rectangle((xs["rep"] - 0.05, y - 0.2), 1.35, 0.4, color=COL[t], lw=0))
        ax.text(xs["rep"] + 0.62, y, REP[t], va="center", ha="center", fontsize=7, color="white", weight="bold",
                family="DejaVu Sans Mono")
        ax.text(xs["tr"], y, tr, va="center", fontsize=7, color="#B00000" if tr.startswith("T4") else "k")
        ax.text(xs["tert"], y, tert, va="center", fontsize=7)
    ax.text(0.15, 0.45, "Topology after Hess et al. (2022). Colored branches and labels: inferred change of the repeat. "
            "Red TR template: template does not match the telomere.", fontsize=6.5, color="#444444")
    colour_key(fig, [(COL[t], KEYLAB[t]) for t in (3, 5, 6, 2)], ncol=4, title="Telomeric repeat (T$n$ = number of T residues)",
               loc="upper center", bbox_to_anchor=(0.5, 0.10))
    save(fig, "Figure_1")


# ---------------- Figure 2: T-tract length distributions ----------------
def read_long(label):
    out = {}
    for l in open(os.path.join(P, "figure_data", f"long_{label}.txt")):
        f = l.rstrip("\n").split("\t")
        if re.match(r"T\d$", f[0]): out[int(f[0][1:])] = float(f[4])
    return out


def read_short(label):
    out = {}
    for l in open(os.path.join(P, "figure_data", f"short_{label}.txt")):
        f = l.rstrip("\n").split("\t")
        if re.match(r"T\d$", f[0]): out[int(f[0][1:])] = float(f[4])
    return out


def figure2():
    panels = [("Spirogloea\nmuscicola", "Illumina", read_short("Spmu_Illumina")),
              ("“Mesotaenium\nendlicherianum”", "Illumina", read_short("Mend_Illumina")),
              ("Closterium\nNIES-64", "PacBio CLR", read_long("Clos_NIES-64")),
              ("Closterium\nNIES-65", "PacBio CLR", read_long("Clos_NIES-65")),
              ("Closterium\nNIES-67", "PacBio CLR", read_long("Clos_NIES-67")),
              ("Closterium\nNIES-68", "PacBio CLR", read_long("Clos_NIES-68")),
              ("Closterium\nNaga37s-1", "PacBio CLR", read_long("Clos_Naga37s-1")),
              ("Closterium\nYama58-4", "PacBio CLR", read_long("Clos_Yama58-4")),
              ("Z. circumcarinatum\nSAG 698-1b", "Illumina", read_short("Zcir_SAG698-1b_Illumina")),
              ("Z. circumcarinatum\nUTEX 1559", "Illumina", read_short("Zcir_UTEX1559_Illumina")),
              ("Z. circumcarinatum\nUTEX 1560", "PacBio HiFi", read_long("Zcir_UTEX1560_HiFi")),
              ("Z. cf. cylindricum\nSAG 698-1a_XF", "PacBio HiFi", read_long("Zcyl_HiFi")),
              ("Z. cf. cylindricum\nSAG 698-1a_XF", "Illumina", read_short("Zcyl_Illumina")),
              ("Spirogyra pratensis\nMZCH 10213", "PacBio HiFi", read_long("Spra_HiFi"))]
    fig, axs = plt.subplots(3, 5, figsize=(7.0, 5.0), sharex=True, sharey=True)
    axs = axs.ravel()
    for i, (name, plat, d) in enumerate(panels):
        ax = axs[i]
        xs = list(range(1, 10)); ys = [d.get(x, 0) for x in xs]
        mode = max(d, key=d.get)
        ax.bar(xs, ys, color=[COL.get(x, "#BBBBBB") if x == mode else "#CCCCCC" for x in xs], width=0.8)
        ax.set_title(name, fontsize=6.6, style="italic", pad=2)
        ax.text(0.97, 0.93, f"{plat}\nT{mode} = {d[mode]:.1f}%", transform=ax.transAxes, ha="right", va="top",
                fontsize=6)
        ax.set_ylim(0, 100); ax.set_xticks(xs)
        ax.tick_params(labelsize=6, length=2)
        for s in ("top", "right"): ax.spines[s].set_visible(False)
    for ax in axs[len(panels):]: ax.axis("off")
    colour_key(axs[-1], [(COL[t], KEYLAB[t]) for t in (2, 3, 5, 6)] + [("#CCCCCC", "other unit lengths")], ncol=1,
               title="Modal unit (colored bar)", loc="center")
    fig.supxlabel("Number of T residues per telomeric unit (T$n$AGGG)", fontsize=7.5, y=0.02)
    fig.supylabel("Units (%)", fontsize=7.5, x=0.02)
    fig.tight_layout(rect=(0.02, 0.03, 1, 1))
    save(fig, "Figure_2")


# ---------------- Figure 3: TR templates + RNA depth ----------------
def depth(tsv, L):
    cov = [0] * L
    for r in csv.DictReader(open(tsv), delimiter="\t"):
        s, e = int(r["start"]), int(r["end"])
        for i in range(max(0, s), min(L, e)): cov[i] += 1
    return cov


def figure3():
    fig = plt.figure(figsize=(7.0, 5.0))
    ax = fig.add_axes([0.02, 0.535, 0.96, 0.445]); ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 5.4)
    ax.text(0, 5.2, "a", weight="bold", fontsize=10)
    ax.text(3.3, 4.75, "Template region of the telomerase RNA (5′→3′)", weight="bold", fontsize=7.5)
    ax.text(7.45, 4.75, "Predicted", weight="bold", fontsize=7.5); ax.text(8.6, 4.75, "Telomere", weight="bold", fontsize=7.5)
    rows = [("Spirogloea muscicola", "CCGAAA", "CCCUAAACCC", "ACCA", 3, 3),
            ("“Mesotaenium endlicherianum”", "CACCCU", "CAACCCUAAACA", "CUU", 3, 3),
            ("Zygnema circumcarinatum", "GCAAGCUAAA", "CCCUAAAACCC", "UACU", 4, 5),
            ("Zygnema cf. cylindricum", "GCAUGCCAAAA", "CCCUAAAACCC", "AUCA", 4, 6)]
    mono = {"family": "DejaVu Sans Mono", "fontsize": 7.5}
    cw = 0.112
    for i, (name, left, core, right, pred, tel) in enumerate(rows):
        y = 4.0 - i * 1.0
        ax.text(0, y, name, style="italic", fontsize=7.5, va="center")
        x = 4.75
        ax.text(x - 0.03, y, "…" + left, color="#777777", va="center", ha="right", **mono)
        # colour the A-run inside CCCUA(n)CCC
        m = re.search(r"CCCU(A+)", core)
        for j, ch in enumerate(core):
            inrun = m and m.start(1) <= j < m.end(1)
            ax.text(x, y, ch, va="center", color=("#B00000" if inrun and pred == 4 else COL[3] if inrun else "k"),
                    weight="bold" if inrun else "normal", **mono)
            x += cw
        ax.text(x, y, right + "…", color="#777777", va="center", **mono)
        for xx, t in ((7.45, pred), (8.6, tel)):
            ax.add_patch(Rectangle((xx - 0.05, y - 0.2), 1.05, 0.4, color=COL[t], lw=0))
            ax.text(xx + 0.47, y, REP[t], ha="center", va="center", color="white", weight="bold",
                    family="DejaVu Sans Mono", fontsize=6.5)
    ax.text(0, 0.2, "Colored letters: the A-run of the template (blue, three A, matching the telomere; red, four A, not matching). "
            "Spirogloea and M. endlicherianum TRs from Fajkus et al. (2021).", fontsize=6.5, color="#444444")
    colour_key(fig, [(COL[t], KEYLAB[t]) for t in (3, 4, 5, 6)], ncol=4, title=None,
               loc="center", bbox_to_anchor=(0.5, 0.505))
    # panel b
    reg = {"Zcir": (500, 782, 1282)}
    s, e, L = reg["Zcir"]
    axb = fig.add_axes([0.08, 0.075, 0.88, 0.30])
    for run, c in (("SRR24205694", "#0072B2"), ("SRR24286561", "#D55E00")):
        cov = depth(os.path.join(P, "tr_search/rnaseq", f"Zcir_{run}.tsv"), L)
        axb.plot([i - s for i in range(L)], cov, color=c, lw=0.9, label=f"{run} (total RNA, rRNA-depleted)")
    axb.axvspan(0, e - s, color="#DDDDDD", zorder=0)
    axb.text((e - s) / 2, axb.get_ylim()[1] * 0.92, "telomerase RNA", ha="center", fontsize=7)
    axb.set_xlabel("Position relative to the 5′ end of the telomerase RNA gene (nt)", fontsize=7.5)
    axb.set_ylabel("RNA-seq read depth", fontsize=7.5)
    axb.legend(fontsize=6.5, frameon=False, loc="upper right")
    axb.tick_params(labelsize=6.5)
    for sp in ("top", "right"): axb.spines[sp].set_visible(False)
    fig.text(0.02, 0.405, "b", weight="bold", fontsize=10)
    fig.text(0.08, 0.42, "Zygnema circumcarinatum SAG 698-1b", style="italic", fontsize=7.5)
    save(fig, "Figure_3")


# ---------------- Figure 4: survey overview ----------------
def figure4():
    import openpyxl
    rows = list(openpyxl.load_workbook(os.path.join(D, "Table_S1_assemblies.xlsx")).active.iter_rows(values_only=True))[2:]
    rows = [r for r in rows if r and r[1]]
    def group(lin):
        l = lin.lower()
        for k, g in [("cyanidio", "Cyanidiophyceae"), ("rhodophyta", "Other Rhodophyta"),
                     ("glauco", "Glaucophyta"), ("zygnemato", "Streptophyta"), ("charophyceae", "Streptophyta"),
                     ("chloro", "Chlorophyta"), ("treboux", "Chlorophyta"), ("ulvo", "Chlorophyta"),
                     ("phaeo", "Phaeophyceae"), ("bacillar", "Bacillariophyta"), ("hapto", "Haptophyta"),
                     ("dino", "Dinophyceae"), ("symbiod", "Dinophyceae"), ("eugleno", "Euglenophyceae"),
                     ("chromerid", "Chromerida")]:
            if k in l: return g
        return lin
    order = ["Chlorophyta", "Streptophyta", "Glaucophyta", "Cyanidiophyceae", "Other Rhodophyta",
             "Phaeophyceae", "Bacillariophyta", "Haptophyta", "Dinophyceae", "Euglenophyceae", "Chromerida"]
    recs = sorted(rows, key=lambda r: (order.index(group(r[0])), r[1]))
    def repcol(rep):
        r = str(rep)
        if "TTTTTTAGGG" in r: return COL[6], "TTTTTTAGGG"
        if "TTTTTAGGG" in r: return COL[5], "TTTTTAGGG"
        if "TTTATTAGGG" in r: return "#56B4E9", "TTTATTAGGG"
        if "AATGGGGGG" in r: return "#56B4E9", "AATGGGGGG"
        if "TATAGGG" in r: return "#F0E442", "TATAGGG (1–3 copies)"
        if "TTGGGG" in r: return "#E69F00", "TTAGGG/TTGGGG"
        if "TTTAGGG" in r: return COL[3], "TTTAGGG"
        if "TTAGGG" in r: return COL[2], "TTAGGG"
        return None, "none at ends" if "none" in r else "not tested"
    n = len(recs)
    fig, ax = plt.subplots(figsize=(7.0, 0.16 * n + 0.9))
    ax.set_xlim(0, 10); ax.set_ylim(-0.8, n + 0.8); ax.axis("off"); ax.invert_yaxis()
    hdr = {0.0: "Lineage", 2.1: "Species", 5.35: "Telomeric repeat", 7.4: "Vs published record", 9.55: "TERT"}
    for x, t in hdr.items(): ax.text(x, -0.5, t, weight="bold", fontsize=7)
    last = None
    for i, r in enumerate(recs):
        y = i + 0.5
        g = group(r[0])
        if g != last:
            ax.text(0, y, g, fontsize=6.5, va="center", weight="bold")
            if last is not None: ax.plot([0, 10], [i, i], color="#DDDDDD", lw=0.5)
            last = g
        w = re.sub(r"\s*[('].*", "", str(r[1])).split()  # drop strain notes in () or quotes
        name = " ".join(w[:3] if len(w) > 1 and w[1] in ("sp.", "cf.") else w[:2])
        ax.text(2.1, y, name, fontsize=6.3, va="center", style="italic")
        c, lab = repcol(r[4])
        if c:
            ax.add_patch(Rectangle((5.3, y - 0.38), 1.95, 0.76, color=c, lw=0))
        ax.text(5.35, y, lab, fontsize=6, va="center", family="DejaVu Sans Mono",
                color="white" if c in (COL[3], COL[5], COL[6]) else "k")
        ax.text(7.4, y, str(r[7]), fontsize=6.3, va="center")
        tert = str(r[8])
        ax.text(9.6, y, "●" if tert.startswith("Detected") else ("–" if "Not searched" in tert else "○"),
                fontsize=7, va="center")
    ax.text(0, n + 0.6, "● TERT detected; ○ not detected (not evidence of absence); – not searched. "
           , fontsize=6.2, color="#444444")
    colour_key(ax, [(COL[3], "TTTAGGG"), (COL[2], "TTAGGG (incl. TTAGGG/TTGGGG)"), (COL[5], "TTTTTAGGG"),
                     (COL[6], "TTTTTTAGGG"), ("#56B4E9", "Cyanidiophyceae repeats"), ("#F0E442", "TATAGGG, tentative")],
               ncol=3, title="Telomeric repeat (no color: none at assembly ends or not tested)",
               loc="upper center", bbox_to_anchor=(5, n + 1.0), bbox_transform=ax.transData)
    save(fig, "Figure_4")


if __name__ == "__main__":
    figure1(); figure2(); figure3(); figure4()
