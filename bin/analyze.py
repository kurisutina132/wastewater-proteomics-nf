#!/usr/bin/env python3
"""
Downstream analysis of Sage output for wastewater proteomics.

Steps
  1. QC per sample (PSMs / peptides at FDR, mass error, missed cleavages)
  2. Peptide taxonomy: organism-specific / group-specific / ambiguous
  3. Protein quantification (Top3 of unique peptides from Sage LFQ)
  4. Taxonomic composition per sample (share of quantified signal)
  5. Human biomarker + immunoglobulin panel
  6. PCA and, if there are >=2 groups with >=2 samples, a Welch t-test + BH
  7. Self-contained HTML report
"""
import argparse
import base64
import io
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

PALETTE = {"Human": "#4C72B0", "Gut bacteria": "#55A868", "Livestock": "#C44E52",
           "Rodent": "#8172B2", "Contaminant": "#937860", "Ambiguous": "#BBBBBB",
           "Other": "#DD8452"}


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def read_fasta_headers(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                h = line[1:].strip()
                pid = h.split()[0]
                os_ = re.search(r"OS=(.+?)(?: OX=| GN=| PE=| SV=|$)", h)
                gn = re.search(r"GN=(\S+)", h)
                rows.append((pid, os_.group(1) if os_ else "unknown",
                             gn.group(1) if gn else ""))
    return pd.DataFrame(rows, columns=["protein_id", "os", "gene"])


def default_group(os_):
    if os_.startswith("Homo sapiens"):
        return "Human"
    if os_ == "contaminant":
        return "Contaminant"
    return "Other"


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(q, 0, 1)
    return out


def fig_to_b64(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def save(fig, outdir, name):
    fig.savefig(outdir / f"{name}.png", dpi=150, bbox_inches="tight")
    return fig_to_b64(fig)


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--psms", required=True)
    ap.add_argument("--lfq", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--samplesheet", required=True, help="sample,group,mzml")
    ap.add_argument("--taxmap", default=None)
    ap.add_argument("--panel", required=True)
    ap.add_argument("--fdr", type=float, default=0.01)
    ap.add_argument("--min_peptides", type=int, default=2)
    ap.add_argument("--title", default="Wastewater proteomics report")
    ap.add_argument("--outdir", default="results")
    a = ap.parse_args()

    out = Path(a.outdir)
    (out / "figures").mkdir(parents=True, exist_ok=True)
    figs = out / "figures"
    imgs = {}

    ss = pd.read_csv(a.samplesheet)
    file2sample = dict(zip(ss.mzml, ss["sample"]))
    sample2group = dict(zip(ss["sample"], ss["group"]))
    samples = list(ss["sample"])

    # ---------------- protein annotation --------------------------------
    ann = read_fasta_headers(a.fasta)
    if a.taxmap and Path(a.taxmap).exists():
        tm = pd.read_csv(a.taxmap, sep="\t")
        ann = ann.merge(tm, on="protein_id", how="left")
        ann["organism"] = ann["organism"].fillna(ann["os"])
        ann["group"] = ann["group"].fillna(ann["os"].map(default_group))
    else:
        ann["organism"] = ann["os"]
        ann["group"] = ann["os"].map(default_group)
    ann = ann.set_index("protein_id")

    # ---------------- PSMs / QC -----------------------------------------
    psm = pd.read_csv(a.psms, sep="\t")
    psm["sample"] = psm["filename"].map(file2sample).fillna(psm["filename"])
    tgt = psm[(psm.label == 1)]
    ok = tgt[(tgt.spectrum_q <= a.fdr) & (tgt.peptide_q <= a.fdr)]
    qc = (ok.groupby("sample")
            .agg(psms=("psm_id", "size"), peptides=("peptide", "nunique"),
                 median_precursor_ppm=("precursor_ppm", "median"),
                 pct_missed_cleavage=("missed_cleavages", lambda x: 100 * (x > 0).mean()))
            .reindex(samples).fillna(0))
    qc["decoy_psms_reported"] = psm[psm.label == -1].groupby("sample").size().reindex(samples).fillna(0)
    qc.to_csv(out / "qc_per_sample.tsv", sep="\t")

    fig, axes = plt.subplots(1, 2, figsize=(10, 3.5))
    qc[["psms", "peptides"]].plot.bar(ax=axes[0], color=["#4C72B0", "#55A868"])
    axes[0].set_title(f"Identifications at {a.fdr:.0%} FDR"); axes[0].set_xlabel("")
    axes[1].hist(ok.precursor_ppm, bins=40, color="#4C72B0")
    axes[1].set_title("Precursor mass error (ppm)")
    imgs["qc"] = save(fig, figs, "qc")

    # ---------------- peptide taxonomy ----------------------------------
    peps = ok[["peptide", "proteins"]].drop_duplicates("peptide").copy()

    def classify(prot_str):
        ids = [p for p in str(prot_str).split(";") if p and not p.startswith("rev_")]
        orgs = {ann.organism.get(p, "unknown") for p in ids}
        grps = {ann.group.get(p, "Other") for p in ids}
        if len(orgs) == 1:
            return "organism-specific", next(iter(orgs)), next(iter(grps))
        if len(grps) == 1:
            return "group-specific", "multiple", next(iter(grps))
        return "ambiguous", "multiple", "Ambiguous"

    cls = peps["proteins"].map(classify)
    peps["specificity"] = [c[0] for c in cls]
    peps["organism"] = [c[1] for c in cls]
    peps["group"] = [c[2] for c in cls]
    peps["unique_protein"] = peps["proteins"].map(lambda s: len(str(s).split(";")) == 1)
    peps.to_csv(out / "peptide_taxonomy.tsv", sep="\t", index=False)

    # ---------------- LFQ -------------------------------------------------
    lfq = pd.read_csv(a.lfq, sep="\t")
    qcols = [c for c in lfq.columns if c in file2sample]
    if not qcols:
        raise SystemExit("LFQ columns do not match samplesheet mzml names: "
                         f"{list(lfq.columns)[6:]} vs {list(file2sample)}")
    lfq = lfq[lfq.q_value <= a.fdr]
    lfq = lfq[lfq.peptide.isin(set(ok.peptide))]
    lfq = lfq.rename(columns=file2sample)
    # one row per peptide sequence (sum charge states if present)
    pep_int = lfq.groupby(["peptide", "proteins"], as_index=False)[samples].sum()
    pep_int[samples] = pep_int[samples].replace(0, np.nan)
    pep_int = pep_int.merge(peps[["peptide", "group", "organism", "specificity"]],
                            on="peptide", how="left")
    pep_int.to_csv(out / "peptide_intensities.tsv", sep="\t", index=False)

    # taxonomic composition (all quantified peptides, by group)
    comp = pep_int.groupby("group")[samples].sum(min_count=1).fillna(0)
    comp_pct = 100 * comp / comp.sum()
    comp_pct.to_csv(out / "taxonomic_composition_pct.tsv", sep="\t")
    org = (pep_int[pep_int.specificity == "organism-specific"]
           .groupby("organism")[samples].sum(min_count=1).fillna(0))
    org.to_csv(out / "organism_intensity.tsv", sep="\t")

    fig, ax = plt.subplots(figsize=(max(5, 0.8 * len(samples) + 2), 4))
    bottom = np.zeros(len(samples))
    for g in comp_pct.index:
        ax.bar(samples, comp_pct.loc[g].values, bottom=bottom, label=g,
               color=PALETTE.get(g, "#999999"))
        bottom += comp_pct.loc[g].values
    ax.set_ylabel("% of quantified peptide signal")
    ax.set_title("Source of the protein signal")
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", frameon=False)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    imgs["taxonomy"] = save(fig, figs, "taxonomic_composition")

    # ---------------- protein quantification (Top3 unique) --------------
    uniq = pep_int[pep_int.proteins.map(lambda s: len(str(s).split(";")) == 1)].copy()
    uniq["protein_id"] = uniq["proteins"]
    rows = []
    for pid, g in uniq.groupby("protein_id"):
        if g.peptide.nunique() < a.min_peptides:
            continue
        vals = {}
        for s in samples:
            v = g[s].dropna().sort_values(ascending=False).head(3)
            vals[s] = v.mean() if len(v) else np.nan
        vals["protein_id"] = pid
        vals["n_peptides"] = g.peptide.nunique()
        rows.append(vals)
    prot = pd.DataFrame(rows)
    if prot.empty:
        raise SystemExit("No protein passed the filters; lower --min_peptides or check FDR.")
    prot = prot.set_index("protein_id")
    prot = prot.join(ann[["gene", "organism", "group"]], how="left")
    prot[["gene", "organism", "group", "n_peptides"] + samples].to_csv(
        out / "protein_matrix.tsv", sep="\t")
    log2 = np.log2(prot[samples])

    # ---------------- biomarker / immunoglobulin panel ------------------
    panel = pd.read_csv(a.panel, sep="\t")
    human = prot[prot.group == "Human"]
    bm = human[human.gene.isin(panel.gene)].copy()
    bm = bm.merge(panel, left_on="gene", right_on="gene", how="left")
    bm = bm[["gene", "category", "rationale", "n_peptides"] + samples]
    bm.to_csv(out / "biomarker_panel.tsv", sep="\t", index=False)
    missing = sorted(set(panel.gene) - set(bm.gene))
    if len(bm):
        mat = np.log2(bm.set_index("gene")[samples])
        mat = mat.sub(mat.mean(axis=1), axis=0)  # row-centred
        fig, ax = plt.subplots(figsize=(max(5, 0.7 * len(samples) + 2), 0.35 * len(mat) + 1.5))
        im = ax.imshow(mat.values, aspect="auto", cmap="RdBu_r",
                       vmin=-np.nanmax(abs(mat.values)), vmax=np.nanmax(abs(mat.values)))
        ax.set_yticks(range(len(mat))); ax.set_yticklabels(mat.index)
        ax.set_xticks(range(len(samples))); ax.set_xticklabels(samples, rotation=45, ha="right")
        ax.set_title("Human biomarker & Ig panel (log2, row-centred)")
        fig.colorbar(im, ax=ax, shrink=0.7)
        imgs["panel"] = save(fig, figs, "biomarker_panel")

    # ---------------- PCA -------------------------------------------------
    keep = log2.notna().mean(axis=1) >= 0.5
    X = log2[keep].copy()
    X = X.apply(lambda r: r.fillna(np.nanmin(r) - 1), axis=1)  # left-censored imputation
    pca_txt = "not enough proteins/samples"
    if X.shape[0] >= 3 and X.shape[1] >= 3:
        Z = (X.T - X.T.mean()).values
        U, S, Vt = np.linalg.svd(Z, full_matrices=False)
        pcs = U * S
        var = S ** 2 / np.sum(S ** 2)
        fig, ax = plt.subplots(figsize=(5, 4))
        groups = [sample2group[s] for s in samples]
        for g in sorted(set(groups)):
            idx = [i for i, x in enumerate(groups) if x == g]
            ax.scatter(pcs[idx, 0], pcs[idx, 1], label=g, s=60)
        for i, s in enumerate(samples):
            ax.annotate(s, (pcs[i, 0], pcs[i, 1]), fontsize=7, xytext=(3, 3),
                        textcoords="offset points")
        ax.set_xlabel(f"PC1 ({var[0]:.0%})"); ax.set_ylabel(f"PC2 ({var[1]:.0%})")
        ax.legend(frameon=False); ax.set_title("PCA of protein abundances")
        imgs["pca"] = save(fig, figs, "pca")
        pca_txt = f"{X.shape[0]} proteins, PC1 {var[0]:.0%}, PC2 {var[1]:.0%}"

    # ---------------- group comparison -----------------------------------
    grp_counts = pd.Series(sample2group).value_counts()
    valid = grp_counts[grp_counts >= 2].index.tolist()
    diff_txt = "skipped (needs >=2 groups with >=2 samples each)"
    if len(valid) >= 2:
        g1, g2 = valid[:2]
        s1 = [s for s in samples if sample2group[s] == g1]
        s2 = [s for s in samples if sample2group[s] == g2]
        D = log2[(log2[s1].notna().sum(axis=1) >= 2) & (log2[s2].notna().sum(axis=1) >= 2)]
        if len(D) >= 3:
            t, p = stats.ttest_ind(D[s2], D[s1], axis=1, equal_var=False)
            res = pd.DataFrame({"log2FC": D[s2].mean(axis=1) - D[s1].mean(axis=1),
                                "p": p}, index=D.index)
            res["q_BH"] = bh(res.p.fillna(1))
            res = res.join(prot[["gene", "organism", "group"]]).sort_values("p")
            res.to_csv(out / f"diff_{g2}_vs_{g1}.tsv", sep="\t")
            fig, ax = plt.subplots(figsize=(5, 4))
            c = res.group.map(PALETTE).fillna("#999999")
            ax.scatter(res.log2FC, -np.log10(res.p), c=c, s=18, alpha=0.8)
            ax.axhline(-np.log10(0.05), ls="--", c="grey", lw=0.8)
            ax.set_xlabel(f"log2 FC ({g2} / {g1})"); ax.set_ylabel("-log10 p")
            ax.set_title("Welch t-test (exploratory)")
            imgs["volcano"] = save(fig, figs, "volcano")
            diff_txt = (f"{g2} vs {g1}: {len(res)} proteins tested, "
                        f"{(res.q_BH <= 0.05).sum()} at BH q<=0.05")

    # ---------------- report ---------------------------------------------
    def table(df, n=30):
        return df.head(n).to_html(float_format=lambda x: f"{x:,.3g}", border=0,
                                  classes="tbl")

    summary = {
        "Samples": len(samples),
        "PSMs at FDR": int(len(ok)),
        "Peptides at FDR": int(ok.peptide.nunique()),
        "Proteins quantified (>= %d unique peptides)" % a.min_peptides: int(len(prot)),
        "Organism-specific peptides": int((peps.specificity == "organism-specific").sum()),
        "Ambiguous peptides": int((peps.specificity == "ambiguous").sum()),
        "PCA": pca_txt,
        "Group comparison": diff_txt,
        "Panel proteins not detected": ", ".join(missing) if missing else "none",
    }
    img_html = "".join(
        f'<figure><img src="data:image/png;base64,{b}"/><figcaption>{k}</figcaption></figure>'
        for k, b in imgs.items())
    html = f"""<!DOCTYPE html><html><head><meta charset="utf-8"><title>{a.title}</title>
<style>
body{{font-family:system-ui,sans-serif;max-width:1100px;margin:2rem auto;padding:0 1rem;color:#222}}
h1{{font-size:1.6rem}} h2{{font-size:1.2rem;border-bottom:1px solid #ddd;padding-bottom:.2rem;margin-top:2rem}}
.tbl{{border-collapse:collapse;font-size:.8rem}} .tbl td,.tbl th{{padding:3px 8px;border-bottom:1px solid #eee;text-align:right}}
figure{{display:inline-block;margin:.5rem}} img{{max-width:520px}} figcaption{{font-size:.75rem;color:#666;text-align:center}}
.note{{background:#f6f8fa;padding:.6rem 1rem;border-left:3px solid #4C72B0;font-size:.85rem}}
</style></head><body>
<h1>{a.title}</h1>
<p class="note">Search: Sage, target-decoy at {a.fdr:.0%} PSM and peptide FDR.
Quantification: Sage LFQ (MS1), protein = mean of Top3 unique peptides.
Taxonomy: peptide-level assignment against the search database (no LCA beyond database content).</p>
<h2>Summary</h2>{pd.Series(summary).to_frame("value").to_html(border=0, classes="tbl")}
<h2>Figures</h2>{img_html}
<h2>QC per sample</h2>{table(qc)}
<h2>Taxonomic composition (% signal)</h2>{table(comp_pct)}
<h2>Human biomarker & immunoglobulin panel</h2>{table(bm)}
<h2>Top proteins by mean intensity</h2>{table(prot.assign(mean=prot[samples].mean(axis=1)).sort_values('mean', ascending=False)[['gene','organism','group','n_peptides','mean']])}
</body></html>"""
    (out / "report.html").write_text(html)
    print(pd.Series(summary).to_string())


if __name__ == "__main__":
    main()
