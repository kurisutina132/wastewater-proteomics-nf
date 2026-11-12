#!/usr/bin/env python3
"""
Build docs/index.html (project page for GitHub Pages) from the files in
assets/ and nextflow.config, so the page never drifts from the pipeline.

    python docs/build_page.py
"""
import csv
import html
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
REPO = "https://github.com/kurisutina132/wastewater-proteomics-nf"


def read_tsv(p):
    return list(csv.DictReader(open(p), delimiter="\t"))


def config_params():
    txt = (ROOT / "nextflow.config").read_text()
    block = re.search(r"^params\s*\{(.*?)^\}", txt, re.S | re.M).group(1)
    out = {}
    for line in block.splitlines():
        m = re.match(r"\s*(\w+)\s*=\s*(.+?)\s*$", line)
        if m:
            out[m.group(1)] = m.group(2).strip("'\"")
    return out


def esc(x):
    return html.escape(str(x))


def table(rows, cols, labels=None, cls=""):
    labels = labels or cols
    head = "".join(f"<th>{esc(l)}</th>" for l in labels)
    body = "".join("<tr>" + "".join(f"<td>{c(r) if callable(c) else esc(r[c])}</td>"
                                     for c in cols) + "</tr>" for r in rows)
    return f'<table class="{cls}"><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>'


def main():
    prot = read_tsv(ROOT / "assets/proteomes.tsv")
    panel = read_tsv(ROOT / "assets/biomarker_panel.tsv")
    p = config_params()

    n_min = sum(r["set"] == "minimal" for r in prot)
    groups = {}
    for r in prot:
        groups.setdefault(r["group"], []).append(r)

    group_color = {"Human": "var(--human)", "Gut bacteria": "var(--bact)",
                   "Livestock": "var(--live)", "Rodent": "var(--rod)"}
    db_rows = [dict(r, link=f'<a href="https://www.uniprot.org/proteomes/{esc(r["proteome"])}">'
                            f'{esc(r["proteome"])}</a>') for r in prot]
    db_table = table(
        db_rows,
        [lambda r: r["link"], lambda r: f"<i>{esc(r['label'])}</i>",
         lambda r: f'<span class="dot" style="background:{group_color.get(r["group"], "#999")}"></span>{esc(r["group"])}',
         lambda r: '<span class="tag on">minimal</span>' if r["set"] == "minimal"
         else '<span class="tag">full</span>',
         lambda r: f"<code>{esc(r['query'])}</code>"],
        ["Proteome", "Organism", "Group", "Set", "UniProt query"])

    cats = {}
    for r in panel:
        cats.setdefault(r["category"], []).append(r)
    panel_html = ""
    for cat, rows in cats.items():
        chips = "".join(f'<div class="chip"><b>{esc(r["gene"])}</b><span>{esc(r["rationale"])}</span></div>'
                        for r in rows)
        panel_html += f'<div class="cat"><h4>{esc(cat)}</h4><div class="chips">{chips}</div></div>'

    params = [
        ("Search engine", "Sage (open source, fragment-index search)"),
        ("Enzyme", f"Trypsin/P, up to {p['missed_cleavages']} missed cleavages, 7–30 aa"),
        ("Precursor / fragment tolerance", f"±{p['precursor_ppm']} ppm / ±{p['fragment_ppm']} ppm"),
        ("Modifications", "Carbamidomethyl C (fixed, optional) · Oxidation M (variable)"),
        ("FDR", f"Target-decoy (reversed), {float(p['fdr']):.0%} at PSM and peptide level"),
        ("Quantification", "MS1 label-free (Sage LFQ); protein = mean of Top3 unique peptides"),
        ("Protein reporting", f"≥ {p['min_peptides']} unique peptides, no razor peptides"),
        ("Taxonomy", "Peptide → organism / group when all matches agree; otherwise ambiguous"),
        ("Statistics", "PCA (left-censored imputation) · Welch t-test + Benjamini–Hochberg"),
        ("Default subset", f"{p['max_files']} RAW files from PRIDE, ≤ {p['max_file_mb']} MB each"),
    ]
    params_html = "".join(f"<tr><th>{esc(k)}</th><td>{esc(v)}</td></tr>" for k, v in params)

    outputs = [
        ("report.html", "Self-contained report with every figure and table"),
        ("qc_per_sample.tsv", "PSMs and peptides at FDR, mass error, missed cleavages"),
        ("peptide_taxonomy.tsv", "Specificity of each peptide and the organism or group it points to"),
        ("taxonomic_composition_pct.tsv", "% of the quantified signal by source"),
        ("protein_matrix.tsv", "Top3 intensities per protein and sample"),
        ("biomarker_panel.tsv", "Human biomarker and immunoglobulin panel"),
        ("diff_&lt;g2&gt;_vs_&lt;g1&gt;.tsv", "log2FC, p-value, BH q-value"),
        ("database/db_manifest.tsv", "Entries per proteome actually downloaded"),
    ]
    out_html = "".join(f"<tr><td><code>{k}</code></td><td>{esc(v)}</td></tr>" for k, v in outputs)

    steps = [
        ("PRIDE", "RAW subset", "fetch_pride.py"),
        ("Convert", "RAW → mzML", "ThermoRawFileParser"),
        ("Database", f"{n_min} proteomes + cRAP", "build_db.py"),
        ("Search", "PSMs · FDR · LFQ", "Sage"),
        ("Taxonomy", "peptide → source", "analyze.py"),
        ("Quant", "Top3 proteins", "analyze.py"),
        ("Panel", "biomarkers · Ig", "analyze.py"),
        ("Report", "PCA · stats · HTML", "analyze.py"),
    ]
    flow = "".join(f'<div class="step"><div class="n">{i+1}</div><b>{a}</b><span>{b}</span>'
                   f'<code>{c}</code></div>' for i, (a, b, c) in enumerate(steps))

    figs = [
        ("taxonomic_composition", "Source of the protein signal",
         "The synthetic design puts 3× more human protein in site B; the pipeline recovers it."),
        ("biomarker_panel", "Human biomarker & immunoglobulin panel", "Row-centred log2 intensities."),
        ("pca", "PCA of protein abundances", "Samples separate by site."),
        ("qc", "Identification QC", "PSMs and peptides at 1% FDR, precursor mass error."),
        ("volcano", "Site B vs site A", "Exploratory Welch t-test, coloured by source."),
    ]
    fig_html = "".join(
        f'<figure><img src="img/{f}.png" alt="{esc(t)}" loading="lazy"/>'
        f'<figcaption><b>{esc(t)}</b> {esc(c)}</figcaption></figure>' for f, t, c in figs)

    page = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"/>
<meta name="viewport" content="width=device-width,initial-scale=1"/>
<title>wastewater-proteomics-nf</title>
<style>
:root{{--ink:#1d2433;--mute:#5d6678;--line:#e3e7ee;--bg:#fafbfc;--acc:#2f5d8a;
--human:#4C72B0;--bact:#55A868;--live:#C44E52;--rod:#8172B2}}
*{{box-sizing:border-box}}
body{{margin:0;font:15px/1.6 system-ui,-apple-system,Segoe UI,Roboto,sans-serif;color:var(--ink);background:#fff}}
header{{background:linear-gradient(135deg,#1f3b5a,#2f5d8a);color:#fff;padding:3.2rem 1.5rem 2.6rem}}
header .w,main{{max-width:1080px;margin:0 auto}}
header h1{{margin:0;font-size:2rem;letter-spacing:-.01em}}
header p{{margin:.6rem 0 0;max-width:720px;opacity:.92}}
.btns{{margin-top:1.2rem;display:flex;gap:.6rem;flex-wrap:wrap}}
.btn{{color:#fff;border:1px solid rgba(255,255,255,.5);padding:.4rem .9rem;border-radius:6px;text-decoration:none;font-size:.9rem}}
.btn:hover{{background:rgba(255,255,255,.12)}}
main{{padding:0 1.5rem 4rem}}
section{{margin-top:2.8rem}}
h2{{font-size:1.25rem;margin:0 0 .3rem}}
.lead{{color:var(--mute);margin:0 0 1rem;max-width:760px}}
.flow{{display:grid;grid-template-columns:repeat(auto-fit,minmax(120px,1fr));gap:.6rem}}
.step{{border:1px solid var(--line);border-radius:8px;padding:.7rem .7rem .6rem;background:var(--bg);display:flex;flex-direction:column;gap:.15rem;position:relative}}
.step .n{{position:absolute;top:.5rem;right:.6rem;font-size:.75rem;color:var(--mute)}}
.step span{{font-size:.85rem;color:var(--mute)}}
.step code{{font-size:.72rem;color:var(--acc);overflow-wrap:anywhere}}
table{{border-collapse:collapse;width:100%;font-size:.86rem}}
th,td{{text-align:left;padding:.45rem .6rem;border-bottom:1px solid var(--line);vertical-align:top}}
thead th{{font-weight:600;color:var(--mute);font-size:.78rem;text-transform:uppercase;letter-spacing:.03em}}
.kv th{{width:240px;color:var(--mute);font-weight:500}}
code{{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:.8rem}}
.dot{{display:inline-block;width:9px;height:9px;border-radius:50%;margin-right:6px}}
.tag{{font-size:.72rem;border:1px solid var(--line);border-radius:10px;padding:.05rem .5rem;color:var(--mute)}}
.tag.on{{border-color:var(--acc);color:var(--acc)}}
.cards{{display:grid;grid-template-columns:repeat(auto-fit,minmax(300px,1fr));gap:1rem}}
.card{{border:1px solid var(--line);border-radius:10px;padding:1rem 1.1rem}}
.card h3{{margin:0 0 .3rem;font-size:1rem}}
.card .meta{{font-size:.8rem;color:var(--mute)}}
.card p{{margin:.5rem 0 0;font-size:.9rem}}
.cat h4{{margin:1rem 0 .4rem;font-size:.85rem;color:var(--mute);text-transform:uppercase;letter-spacing:.04em}}
.chips{{display:grid;grid-template-columns:repeat(auto-fit,minmax(240px,1fr));gap:.5rem}}
.chip{{border:1px solid var(--line);border-radius:8px;padding:.5rem .7rem;background:var(--bg)}}
.chip b{{display:block;font-family:ui-monospace,Menlo,monospace}}
.chip span{{font-size:.82rem;color:var(--mute)}}
.warn{{border-left:4px solid #d9a400;background:#fff8e1;padding:.7rem 1rem;border-radius:4px;font-size:.9rem}}
.figs{{display:grid;grid-template-columns:repeat(auto-fit,minmax(320px,1fr));gap:1.2rem}}
figure{{margin:0;border:1px solid var(--line);border-radius:10px;padding:.8rem;background:#fff}}
figure img{{width:100%;height:auto}}
figcaption{{font-size:.82rem;color:var(--mute);margin-top:.4rem}}
pre{{background:#0f1720;color:#e6edf3;padding:1rem;border-radius:8px;overflow:auto;font-size:.82rem}}
ul.lim li{{margin:.3rem 0}}
footer{{border-top:1px solid var(--line);color:var(--mute);font-size:.8rem;padding:1.5rem;text-align:center}}
a{{color:var(--acc)}}
</style></head><body>
<header><div class="w">
<h1>wastewater-proteomics-nf</h1>
<p>A lightweight, reproducible Nextflow pipeline for shotgun proteomics of municipal wastewater:
from public PRIDE data to the source of the protein signal, human population biomarkers and
immunoglobulins, in a single command on a laptop.</p>
<div class="btns"><a class="btn" href="{REPO}">Code on GitHub</a>
<a class="btn" href="example_report_synthetic.html">Example report</a>
<a class="btn" href="#run">How to run</a></div>
</div></header>
<main>

<section><h2>Workflow</h2>
<p class="lead">Every step runs as a Nextflow process. The database and the biomarker panel are plain TSV files in <code>assets/</code>, so they can be changed without touching code.</p>
<div class="flow">{flow}</div></section>

<section><h2>Questions it answers</h2>
<div class="cards">
<div class="card"><h3>Where does the protein come from?</h3><p>Share of the quantified signal from human, gut bacteria, livestock and rodent proteins, counting only peptides whose matches all agree on the source.</p></div>
<div class="card"><h3>Which human markers are measurable?</h3><p>Uromodulin and albumin (urine), calprotectin and lactoferrin (gut inflammation), and a secretory IgA / immunoglobulin block for immunoproteomics.</p></div>
<div class="card"><h3>Do sites or dates differ?</h3><p>PCA and an exploratory Welch t-test with Benjamini–Hochberg correction whenever there are at least two groups of two samples.</p></div>
</div></section>

<section><h2>Public data</h2>
<p class="lead">Default input: a small, evenly spaced subset of RAW files, so the analysis stays light.</p>
<div class="cards">
<div class="card"><h3>PXD059455 <span class="tag on">default</span></h3>
<div class="meta">PRIDE · Tugui et al., ISME Communications 2025 · TU Delft</div>
<p>Untreated municipal wastewater from two urban locations, sampled over several months (Nov 2023 – Feb 2024), measured on an Orbitrap. Microbial and human proteins.</p>
<p><a href="https://www.ebi.ac.uk/pride/archive/projects/PXD059455">PRIDE</a> ·
<a href="https://doi.org/10.1093/ismeco/ycaf243">paper</a></p></div>
<div class="card"><h3>Your own data</h3>
<div class="meta">samplesheet.csv · sample,file,group</div>
<p>Any Thermo <code>.raw</code> or <code>.mzML</code>, local or by URL. You can also bring your own FASTA and taxonomy map with <code>--fasta</code> and <code>--taxmap</code>.</p></div>
</div></section>

<section><h2>Search database</h2>
<p class="lead">UniProt reference proteomes plus cRAP contaminants. <code>minimal</code> (default) = {n_min} proteomes; <code>full</code> = {len(prot)}. The build stops if any query returns zero entries.</p>
{db_table}</section>

<section><h2>Biomarker &amp; immunoglobulin panel</h2>
<p class="lead">Curated human proteins reported per sample (<code>assets/biomarker_panel.tsv</code>).</p>
{panel_html}</section>

<section><h2>Method parameters</h2>
<table class="kv"><tbody>{params_html}</tbody></table></section>

<section id="example"><h2>Example output</h2>
<p class="warn"><b>Synthetic data.</b> The figures below come from <code>-profile test</code>, which builds random protein sequences and simulated spectra. They show the pipeline working end to end, not real wastewater biology.</p>
<div class="figs">{fig_html}</div></section>

<section><h2>Outputs</h2>
<table><thead><tr><th>File</th><th>Content</th></tr></thead><tbody>{out_html}</tbody></table></section>

<section id="run"><h2>How to run</h2>
<pre># synthetic test (no download)
nextflow run kurisutina132/wastewater-proteomics-nf -profile test,conda

# public wastewater data: 6 RAW files from PRIDE
nextflow run kurisutina132/wastewater-proteomics-nf -profile conda \\
    --pride PXD059455 --max_files {p['max_files']}

# your own samples
nextflow run kurisutina132/wastewater-proteomics-nf -profile conda \\
    --samplesheet samplesheet.csv</pre></section>

<section><h2>Limitations</h2>
<ul class="lim">
<li>A small curated database separates human, animal and bacterial signal but under-represents microbial diversity. Bacterial identifications are a lower bound.</li>
<li>"Organism-specific" means specific within this database, not across all of UniProt.</li>
<li>With a few files per group the statistics are exploratory and generate hypotheses, not findings.</li>
</ul></section>

</main>
<footer>MIT licence · Sage (Lazear, J Proteome Res 2023) · ThermoRawFileParser · UniProt · cRAP ·
generated by <code>docs/build_page.py</code></footer>
</body></html>"""
    (DOCS / "index.html").write_text(page)
    print("wrote", DOCS / "index.html")


if __name__ == "__main__":
    main()
