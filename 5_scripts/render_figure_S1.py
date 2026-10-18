"""Figure S1: 18S cladogram of the survey taxa (rebuilt 2026-10-14). Based on phase3_render_figure.py.
UFBoot from the unconstrained ML run where the bipartition is shared, plus
telomere-motif and TERT columns. Writes Figure_S1.{svg,png}.
"""
import csv
from pathlib import Path
from Bio import Phylo
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parents[1] / "data/processed/phylogeny"
OUTGROUP = ("Euglena_gracilis", "Eutreptiella_sp")

def load(f):
    t = Phylo.read(OUT / f, "newick"); t.root_with_outgroup(*OUTGROUP); t.ladderize(reverse=True)
    return t

con, ml = load("constrained.treefile"), load("tree.treefile")
alltips = frozenset(x.name for x in ml.get_terminals())
def split(c):
    s = frozenset(x.name for x in c.get_terminals())
    return min(s, alltips - s, key=lambda z: sorted(z))
sup = {split(c): (c.name.split("/")[-1] if c.name else (f"{c.confidence:.0f}" if c.confidence is not None else None)) for c in con.get_nonterminals()}
ENFORCED = {split(con.common_ancestor(*g)) for g in (
    ("Galdieria_yellowstonensis", "Cyanidium_caldarium", "Cyanidiococcus_yangmingshanensis"),
    ("Gracilaria_vermiculophylla", "Porphyra_dioica", "Porphyra_linearis"),
    ("Chlorella_sorokiniana", "Chlorella_pyrenoidosa", "Auxenochlorella_pyrenoidosa", "Picochlorum_sp_SENEW3"))}

traits = {r["tip"]: r for r in csv.DictReader(open(OUT / "trait_table.tsv"), delimiter="\t")}
COL = {"TTTAGGG": "#2a78d6", "TTAGGG": "#eb6834", "T5AGGG": "#1baf7a", "WEAK_T5": "#1baf7a",
       "TTTATTAGGG": "#eda100", "AATG6": "#e87ba4", "ULTRASHORT": "#4a3aa7"}
TERT = {"Complete": ("●", "#30302e"), "Partial": ("◐", "#30302e"), "Proteome": ("◆", "#30302e"),
        "Not detected": ("○", "#8a8984"), "Not searched": ("–", "#8a8984")}
GROUP_COL = {"Chlorophyta": "#cfe8d9", "Streptophyta": "#e3f0cf", "Rhodophyta": "#f6d9dc",
             "Glaucophyta": "#d9e7f2", "Stramenopiles": "#efe2cc", "Alveolata": "#dcdaf2",
             "Haptophyta": "#f3e7f1", "Euglenozoa": "#e6e5e0"}

tips = con.get_terminals()
y = {t: i for i, t in enumerate(tips)}
depth = {}
def setx(c, d=0):
    depth[c] = d
    for k in c.clades: setx(k, d + 1)
setx(con.root)
maxd = max(depth[t] for t in tips)
X = lambda c: maxd if c.is_terminal() else depth[c]
def sety(c):
    if c.is_terminal(): return y[c]
    y[c] = sum(sety(k) for k in c.clades) / len(c.clades); return y[c]
sety(con.root)

n = len(tips)
plt.rcParams["font.family"] = "DejaVu Sans"
fig, ax = plt.subplots(figsize=(11, 0.24 * n + 1.3))
ax.set_xlim(-0.3, maxd + 19); ax.set_ylim(n, -2.2); ax.axis("off")
plt.rcParams["font.family"] = "DejaVu Sans"
for c in con.get_nonterminals():
    ys = [y[k] for k in c.clades]
    ax.plot([X(c)] * 2, [min(ys), max(ys)], color="#5f5e5a", lw=1)
    for k in c.clades:
        ax.plot([X(c), X(k)], [y[k]] * 2, color="#5f5e5a", lw=1)
    if c is con.root: continue
    s = split(c)
    lab = sup.get(s)
    txt = "C" if s in ENFORCED else (lab if lab and float(lab) >= 70 else "")
    if txt:
        ax.text(X(c) - 0.08, y[c] - 0.18, txt, fontsize=5.5, ha="right", va="bottom",
                color="#b23a3a" if txt == "C" else "#6b6a65")
L0 = maxd + 0.4
for t in tips:
    r = traits[t.name]
    ax.add_patch(plt.Rectangle((L0 - 0.2, y[t] - 0.5), 8.4, 1, color=GROUP_COL[r["group"]], lw=0, zorder=0))
    ax.text(L0, y[t], (t.name.replace("_", " ") + " ").replace(" sp ", " sp. ").strip(), fontsize=7, va="center", style="italic")
    cls = r["motif_class"]
    if cls in COL:
        weak = "weak" in r["motif_evidence"] or cls == "WEAK_T5"
        ax.add_patch(plt.Rectangle((L0 + 8.4, y[t] - 0.38), 4.6, 0.76, facecolor=COL[cls] if not weak else "white",
                                   edgecolor=COL[cls], lw=1.2, alpha=1 if not weak else 1))
        if weak:
            ax.add_patch(plt.Rectangle((L0 + 8.4, y[t] - 0.38), 4.6, 0.76, facecolor=COL[cls], alpha=0.25, lw=0))
        ax.text(L0 + 10.7, y[t], r["motif_label"], fontsize=6.3, ha="center", va="center",
                color="white" if not weak else "#30302e", family="monospace")
    else:
        ax.text(L0 + 10.7, y[t], r["motif_label"], fontsize=6.3, ha="center", va="center", color="#8a8984")
    ax.text(L0 + 13.4, y[t], r["motif_evidence"], fontsize=5.8, va="center", color="#5f5e5a")
    sym, sc = TERT[r["tert"]]
    ax.text(L0 + 17.4, y[t], sym, fontsize=8.5, ha="center", va="center", color=sc)
    if r["tert_note"].startswith("BUSCO"):
        ax.text(L0 + 18.0, y[t], r["tert_note"].replace("BUSCO ", ""), fontsize=5.3, va="center", color="#8a8984")
for x, h in ((L0, "Taxon (lineage band)"), (L0 + 10.7, "Telomere motif"), (L0 + 13.4, "Evidence"), (L0 + 17.4, "TERT")):
    ax.text(x, -1.2, h, fontsize=7.5, weight="bold", ha="center" if x in (L0 + 10.7, L0 + 17.4) else "left")
ax.text(0, n + 0.2, "Cladogram (not to scale). Numbers: UFBoot ≥70 (1000 replicates); C = node enforced by "
        "published-topology constraint (AU p = 0.23). Rooted on Euglenozoa for display.\nOutlined motif = weak support. "
        "TERT: ● complete  ◐ partial  ◆ proteome only  ○ not detected (BUSCO % shown)  – not searched.",
        fontsize=6, va="top", color="#5f5e5a")
fig.savefig(OUT / "Figure_S1.svg", bbox_inches="tight")
fig.savefig(OUT / "Figure_S1.png", dpi=220, bbox_inches="tight")
