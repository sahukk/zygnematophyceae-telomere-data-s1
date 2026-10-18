"""Build Table_S1_assemblies.xlsx (46 survey entries) from the master tracker plus the
per-entry status/evidence decided in the 2026-10 re-examination. Run with telomere_env python."""
import os, sys, openpyxl
from openpyxl.styles import Font, Alignment
D = os.path.dirname(os.path.abspath(__file__))
TR = os.path.join(D, "../../algal_telomere_study/data/processed/master_species_tracker.xlsx")
# key (species prefix) -> (earlier survey call, status, evidence)
M = {
 "Symbiodinium microadriaticum": ("none", "No array at ends", "assembly ends; Hi-C pseudo-chromosomes"),
 "Nitzschia inconspicua": ("TTAGGGTTGGGG (Oliver et al. 2021)", "Confirmed (weak)", "mixed TTAGGG/TTGGGG at ~33/114 assembly ends; PacBio CLR read-end 6-mers ~1.2-1.3x"),
 "Scenedesmus obliquus": ("TTTAGGG", "Confirmed (control)", "T2T assembly; 33/34 ends"),
 "Zygnema circumcarinatum": ("none", "New", "Illumina 66-69% T5; HiFi 75% T5 (3 strains)"),
 "Porphyra dioica": ("none", "No array at ends", "assembly ends"),
 "Porphyra linearis": ("none", "No array at ends", "assembly ends"),
 "Phaeocystis globosa": ("TTAGGG", "Confirmed", "11-16/46 assembly ends"),
 "Gracilaria vermiculophylla": ("none", "No array at ends", "assembly ends (weak)"),
 "Chloromonas typhlos": ("none", "No array at ends", "assembly ends"),
 "Galdieria yellowstonensis": ("AACCCT (TTAGGG family)", "Changed", "111/144 assembly ends; 5,533 PacBio read ends, 0 interiors; as Cho et al. 2023"),
 "Cyanophora paradoxa": ("none", "No array at ends", "contig-level assembly"),
 "Tetraselmis helgolandica": ("not searched", "Not tested", "no public assembly"),
 "Auxenochlorella pyrenoidosa": ("none", "New", "16-17/26 assembly ends"),
 "Penium margaritaceum": ("none", "No array at ends", "fragmented scaffold assembly"),
 "Spirogloea muscicola": ("none", "New", "Illumina SRR9052848: 94.5% T3 (1,345 reads); published TR template T3"),
 "Mesotaenium endlicherianum": ("none", "New (agrees with TRAP)", "Illumina SRR9023141: 95.4% T3; TRAP TTTAGGG (Fulneckova et al. 2013)"),
 "Saccharina japonica": ("none", "New (conflicts with Yang et al. 2017)", "25/114 assembly ends; 430 PacBio read ends, 0 interiors"),
 "Saccharina sculpera": ("none", "New", "37/136 assembly ends"),
 "Undaria pinnatifida": ("TTAGGG", "Confirmed", "16-23/84 assembly ends"),
 "Coccomyxa elongata": ("AACCCT", "Changed", "53 ONT read ends, 0 interiors"),
 "Chlorella sorokiniana": ("TTTAGGG", "Confirmed (control)", "near-T2T; 23/24 ends"),
 "Chlorella pyrenoidosa": ("TTTAGGG", "Confirmed (control)", "near-T2T; 18-26/26 ends"),
 "Picochlorum sp. SENEW3": ("AAAC", "Changed (tentative)", "no end enrichment of AAAC; 291 of 12,468 adapter-marked ONT read ends start with 1-3 TATAGGG (TATAGGG predicted for P. costavermella, Fajkus et al. 2021)"),
 "Ulva compressa": ("none", "No array at ends", "fragmented assembly"),
 "Breviolum minutum": ("none", "New", "141/220 assembly ends; HiFi read-end enrichment ~2.3x"),
 "Effrenium voratum": ("none", "New (weak)", "12/276 assembly ends; also interstitial"),
 "Freudenthalidium endolithicum": ("none", "New (weak)", "38/342 assembly ends"),
 "Symbiodinium pilosum": ("none", "New (weak)", "24/882 assembly ends"),
 "Symbiodinium tridacnidarum": ("none", "New (weak)", "29/314 assembly ends"),
 "Isochrysis galbana": ("none", "No array at ends", "assembly ends"),
 "Prymnesium parvum": ("none", "New", "37/68 assembly ends; 159 ONT read ends, 0 interiors"),
 "Rebecca billardiae": ("none", "New", "43/68 assembly ends"),
 "Euglena gracilis": ("none", "No array at ends", "assembly ends"),
 "Eutreptiella sp.": ("none", "New (weak)", "33/386 assembly ends"),
 "Polarella glacialis": ("none", "No array at ends", "fragmented assembly"),
 "Prorocentrum cordatum": ("none", "No array at ends", "fragmented assembly"),
 "Chromera velia": ("not tested", "Not tested", "fragmented; TERT search only"),
 "Vitrella brassicaformis": ("not tested", "Not tested", "fragmented; TERT search only"),
 "Bryopsis sp.": ("TTTAGGG", "Confirmed", "15/36 assembly ends"),
 "Spirogyra pratensis": ("TTAGGG", "Confirmed", "32/62 assembly ends; HiFi 93% T2, 47 arrays at read ends"),
 "Cyanidium caldarium": ("ACCC", "Changed", "40/40 ends (tolerant); PacBio CLR read-end 6-mers ~2.4x; as Cho et al. 2023"),
 "Cyanidiococcus yangmingshanensis": ("CCCCCCTT", "Changed", "23/40 ends (tolerant); PacBio CLR read-end 6-mers ~1.7x; as Cho et al. 2023"),
 "Zygnema cf. cylindricum": ("(added in this study)", "Added", "HiFi 69% T6; Illumina 71% T6"),
 "Closterium sp.": ("(added in this study)", "Added", "PacBio 63-71% T3 in 6 strains (NIES-64, -65, -67, -68, Naga37s-1, Yama58-4); ONT 86-93% T3 at read ends in NIES-53 and NIES-54; 44/1,674 NIES-67 assembly ends"),
}
PUB = {  # key -> (published record (source), comparison)
 "Symbiodinium microadriaticum": ("TTTAGGG (read-based; TeloBase)", "No array at ends"),
 "Nitzschia inconspicua": ("TTAGGGTTGGGG (Oliver et al. 2021)", "Agrees"),
 "Scenedesmus obliquus": ("TTTAGGG in other Scenedesmus spp. (Fulneckova et al. 2012, 2013)", "New for species"),
 "Zygnema circumcarinatum": ("TTTAGGG by TRAP under this name, strain SAG 698-1a (Fulneckova et al. 2013)", "Conflicts"),
 "Porphyra dioica": ("-", "No array at ends"), "Porphyra linearis": ("-", "No array at ends"),
 "Phaeocystis globosa": ("TTAGGG (read-based; TeloBase)", "Agrees"),
 "Gracilaria vermiculophylla": ("-", "No array at ends"), "Chloromonas typhlos": ("-", "No array at ends"),
 "Galdieria yellowstonensis": ("(TTTATTAGGG)n (Cho et al. 2023; Galdieria: Fulneckova et al. 2013)", "Agrees"),
 "Cyanophora paradoxa": ("-", "No array at ends"), "Tetraselmis helgolandica": ("-", "Not tested"),
 "Auxenochlorella pyrenoidosa": ("TTTAGGG in A. protothecoides (TRAP; Fulneckova et al. 2013)", "New for species"),
 "Penium margaritaceum": ("-", "No array at ends"),
 "Spirogloea muscicola": ("none (TR template predicts TTTAGGG; Fajkus et al. 2021)", "New for genus"),
 "Mesotaenium endlicherianum": ("TTTAGGG (TRAP; Fulneckova et al. 2013)", "Agrees"),
 "Saccharina japonica": ("TTTAGGG (Yang et al. 2017; Liu et al. 2022); TTAGGG read-based in TeloBase", "Conflicts"),
 "Saccharina sculpera": ("none; genus record S. japonica TTTAGGG (Yang et al. 2017), which conflicts with our read-level TTAGGG for S. japonica and with the read-based TeloBase prediction", "New for species"),
 "Undaria pinnatifida": ("TTAGGG (read-based; TeloBase)", "Agrees"),
 "Coccomyxa elongata": ("TTTAGGG in C. subellipsoidea (Fulneckova et al. 2013)", "New for species"),
 "Chlorella sorokiniana": ("TTTAGGG (read-based; TeloBase)", "Agrees"),
 "Chlorella pyrenoidosa": ("TTTAGGG in other Chlorella spp. (TeloBase)", "New for species"),
 "Picochlorum sp. SENEW3": ("TATAGGG in P. costavermella (TR template; Fajkus et al. 2021)", "New for species (tentative)"),
 "Ulva compressa": ("-", "No array at ends"),
 "Breviolum minutum": ("none", "New for genus"), "Effrenium voratum": ("none", "New for genus (weak)"),
 "Freudenthalidium endolithicum": ("none", "New for genus (weak)"),
 "Symbiodinium pilosum": ("TTTAGGG in Symbiodinium (TeloBase)", "New for species (weak)"),
 "Symbiodinium tridacnidarum": ("TTTAGGG in Symbiodinium (TeloBase)", "New for species (weak)"),
 "Isochrysis galbana": ("-", "No array at ends"),
 "Prymnesium parvum": ("none", "New for genus"), "Rebecca billardiae": ("none", "New for genus"),
 "Euglena gracilis": ("TTAGGG (Dooijes et al. 2000; TeloBase)", "No array at ends"),
 "Eutreptiella sp.": ("none (TTAGGG in Euglena)", "New for genus (weak)"),
 "Polarella glacialis": ("-", "No array at ends"), "Prorocentrum cordatum": ("-", "No array at ends"),
 "Chromera velia": ("-", "Not tested"), "Vitrella brassicaformis": ("-", "Not tested"),
 "Bryopsis sp.": ("none", "New for genus"),
 "Spirogyra pratensis": ("none", "New for genus"),
 "Cyanidium caldarium": ("(AATGGGGGG)n (Cho et al. 2023)", "Agrees"),
 "Cyanidiococcus yangmingshanensis": ("(AATGGGGGG)n (Cho et al. 2023)", "Agrees"),
 "Zygnema cf. cylindricum": ("TTTAGGG by TRAP for parent culture SAG 698-1a (Fulneckova et al. 2013)", "Conflicts"),
 "Closterium sp.": ("none (desmid Micrasterias: TTTAGGG by TRAP)", "New for genus"),
}

# taxa excluded from this table
EXCLUDE = ("Mesotaenium kramstae",)
REPEAT_OVERRIDE = {"Spirogloea muscicola": "TTTAGGG"}
LINEAGE_OVERRIDE = {"Penium margaritaceum": "Zygnematophyceae (Streptophyta)"}  # desmid (Desmidiales; Hess et al. 2022)
wb = openpyxl.load_workbook(TR, read_only=True)
rows = list(wb.worksheets[0].iter_rows(values_only=True))[1:]
out = openpyxl.Workbook(); ws = out.active; ws.title = "Table S1"
hdr = ["Lineage", "Species (strain)", "Assembly", "Assembly level", "Telomeric repeat (this study)",
       "Evidence", "Published record (source)", "Comparison with published record", "TERT (domain-based search)"]
ws.append(["Table S1. Telomeric repeats and TERT in the 44 algal genome assemblies after read-level testing, compared with published records."])
ws.append(hdr)
used = set()
for r in rows:
    if r[0].startswith(EXCLUDE): continue
    key = next(k for k in M if r[0].startswith(k)); used.add(key)
    earlier, status, ev = M[key]
    rep = REPEAT_OVERRIDE.get(key, r[7] or "-")
    if str(rep).startswith("No telomeric array"): rep = "none at assembly ends"
    tert = (r[10] or "").split(" (")[0].split(" --")[0].split(";")[0]
    if "split" in (r[10] or ""): tert = "Detected (split across two scaffolds)"
    pubrec, comp = PUB[key]
    name = "Closterium sp. (8 strains; NIES-67 assembly analyzed)" if key == "Closterium sp." else r[0]
    ws.append([LINEAGE_OVERRIDE.get(key, r[1]), name, str(r[2]).split(" (")[0], r[3], rep, ev, pubrec, comp, tert])
assert len(used) == len(M) - 2 == 44, (len(used), set(M) - used)
for c in ws[2]: c.font = Font(bold=True)
ws["A1"].font = Font(bold=True)
for col, w in zip("ABCDEFGHI", [30, 40, 26, 26, 30, 60, 50, 28, 34]):
    ws.column_dimensions[col].width = w
for row in ws.iter_rows(min_row=2):
    for c in row: c.alignment = Alignment(wrap_text=True, vertical="top")
ws.append([]); ws.append(["Status counts are given in the Results; 'none at assembly ends' is not evidence that canonical telomeres are absent."])
p = os.path.join(D, "Table_S1_assemblies.xlsx"); out.save(p)
import collections; print(p); print(collections.Counter(v[1].split(" (")[0] for k,v in PUB.items()))
