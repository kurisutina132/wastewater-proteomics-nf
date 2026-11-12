#!/usr/bin/env python3
"""
Generate a small *synthetic* dataset to test the pipeline end to end
without downloading anything.

Output:
  test_db.fasta          random proteins with UniProt-style headers
                         (human biomarker gene names + gut bacteria)
  <sample>.mzML          MS1 isotope envelopes + HCD-like MS2 (b/y ions)
  test_samplesheet.csv   sample,file,group

The sequences are random: they do NOT correspond to real proteins.
They only exist so the pipeline wiring (search -> FDR -> LFQ -> taxonomy ->
stats -> report) can be checked in CI in a few seconds.
"""
import argparse
import csv
import random
from pathlib import Path

import numpy as np
from psims.mzml.writer import MzMLWriter
from pyteomics import mass, parser

AA = "ADEFGHIKLNPQRSTVWY"  # no C / M -> no mods needed in the test
PROTON = 1.007276

ORGANISMS = [
    # (OS, OX, gene names)
    ("Homo sapiens", 9606, ["UMOD", "PIGR", "IGHA1", "JCHAIN", "IGKC", "S100A8",
                            "S100A9", "LTF", "ALB", "LYZ", "CEACAM5", "KRT1"]),
    ("Escherichia coli (strain K12)", 83333, ["tuf", "groL", "atpA", "ompA"]),
    ("Bacteroides thetaiotaomicron", 818, ["tuf", "groL", "susC"]),
    ("Sus scrofa", 9823, ["ALB", "HBB"]),
]


def random_protein(rng, n_peptides=6):
    """Concatenate random tryptic peptides (each ends in K/R)."""
    peps = []
    for _ in range(n_peptides):
        length = rng.randint(8, 16)
        core = "".join(rng.choice(AA.replace("K", "").replace("R", "").replace("P", ""))
                       for _ in range(length - 1))
        peps.append(core + rng.choice("KR"))
    return "M" + "".join(peps)


def peaks_by(peptide, rng):
    mzs, ints = [], []
    for i in range(1, len(peptide)):
        b = mass.fast_mass(peptide[:i], ion_type="b", charge=1)
        y = mass.fast_mass(peptide[i:], ion_type="y", charge=1)
        for m in (b, y):
            if 150 < m < 2000:
                mzs.append(m)
                ints.append(rng.uniform(2e3, 1e5))
    # a few noise peaks
    for _ in range(10):
        mzs.append(rng.uniform(150, 1500))
        ints.append(rng.uniform(100, 3e3))
    order = np.argsort(mzs)
    return np.array(mzs)[order], np.array(ints)[order]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="test_data")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--background", type=int, default=400,
                    help="proteins in the FASTA that are absent from the samples")
    ap.add_argument("--noise_rate", type=float, default=0.6,
                    help="probability of an unexplained MS2 scan per cycle")
    args = ap.parse_args()
    rng = random.Random(args.seed)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    # ---- FASTA -----------------------------------------------------------
    proteins = []
    acc_n = 0
    for os_, ox, genes in ORGANISMS:
        for g in genes:
            acc_n += 1
            acc = f"T{acc_n:05d}"
            seq = random_protein(rng)
            header = (f">tr|{acc}|{acc}_SYN Synthetic {g} OS={os_} OX={ox} GN={g} PE=1 SV=1")
            proteins.append((header, seq, os_, g))
    # background proteins: in the database but never in the samples
    # (gives the target-decoy competition something realistic to work with)
    n_present = len(proteins)
    for i in range(args.background):
        acc_n += 1
        os_, ox, _ = ORGANISMS[i % len(ORGANISMS)]
        acc = f"T{acc_n:05d}"
        proteins.append((f">tr|{acc}|{acc}_SYN Background protein OS={os_} OX={ox} "
                         f"GN=BG{i} PE=4 SV=1", random_protein(rng, 10), os_, f"BG{i}"))
    groups = {"Homo sapiens": "Human", "Sus scrofa": "Livestock"}
    with open(out / "test_taxmap.tsv", "w") as fh:
        fh.write("protein_id\torganism\tgroup\n")
        for h, _, os_, _ in proteins:
            fh.write(f"{h[1:].split()[0]}\t{os_}\t{groups.get(os_, 'Gut bacteria')}\n")
    with open(out / "test_db.fasta", "w") as fh:
        for h, s, _, _ in proteins:
            fh.write(h + "\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")

    # ---- peptides --------------------------------------------------------
    peptides = []
    for h, s, os_, g in proteins[:n_present]:
        for p in parser.cleave(s, parser.expasy_rules["trypsin"], 0):
            if 7 <= len(p) <= 30 and "M" not in p:
                peptides.append((p, os_, g))
    peptides = sorted(set(peptides))

    # ---- samples: 2 "sites" x 2 timepoints --------------------------------
    samples = [("siteA_t1", "siteA"), ("siteA_t2", "siteA"),
               ("siteB_t1", "siteB"), ("siteB_t2", "siteB")]
    base_abund = {p: rng.uniform(1e5, 1e7) for p, _, _ in peptides}
    rt_of = {p: rng.uniform(5, 55) for p, _, _ in peptides}  # minutes

    with open(out / "test_samplesheet.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sample", "file", "group"])
        for name, grp in samples:
            w.writerow([name, f"{name}.mzML", grp])

    for name, grp in samples:
        srng = random.Random(args.seed + sum(map(ord, name)))
        # site B has more human-derived signal (e.g., more population / less dilution)
        abund = {}
        for p, os_, g in peptides:
            f = 3.0 if (grp == "siteB" and os_ == "Homo sapiens") else 1.0
            abund[p] = base_abund[p] * f * srng.uniform(0.7, 1.3)

        scans = []
        scan_no = 0
        rt_grid = np.arange(1.0, 60.0, 0.05)  # 3 s cycle
        for rt in rt_grid:
            # MS1
            mzs, ints = [], []
            for p, _, _ in peptides:
                d = rt - rt_of[p]
                if abs(d) < 0.4:
                    w_ = np.exp(-0.5 * (d / 0.12) ** 2)
                    mono = (mass.fast_mass(p) + 2 * PROTON) / 2
                    for iso, rel in enumerate((1.0, 0.8, 0.4, 0.15)):
                        mzs.append(mono + iso * 1.003355 / 2)
                        ints.append(abund[p] * w_ * rel)
            for _ in range(20):
                mzs.append(srng.uniform(350, 1500)); ints.append(srng.uniform(1e3, 5e4))
            o = np.argsort(mzs)
            scan_no += 1
            ms1_id = f"controllerType=0 controllerNumber=1 scan={scan_no}"
            scans.append(dict(id=ms1_id, mz=np.array(mzs)[o], it=np.array(ints)[o],
                              level=1, rt=rt))
            # MS2 for peptides at apex in this cycle
            for p, _, _ in peptides:
                if abs(rt - rt_of[p]) < 0.025:
                    mono = (mass.fast_mass(p) + 2 * PROTON) / 2
                    fm, fi = peaks_by(p, srng)
                    scan_no += 1
                    scans.append(dict(id=f"controllerType=0 controllerNumber=1 scan={scan_no}",
                                      mz=fm, it=fi, level=2, rt=rt + 0.001,
                                      prec=dict(mz=mono, charge=2, intensity=abund[p],
                                                scan_id=ms1_id)))
            # unexplained MS2 spectra (random precursor, noise fragments)
            if srng.random() < args.noise_rate:
                pmz = srng.uniform(400, 1200)
                fm = np.sort(np.array([srng.uniform(150, 2 * pmz) for _ in range(150)]))
                fi = np.array([srng.uniform(1e2, 5e4) for _ in range(150)])
                scan_no += 1
                scans.append(dict(id=f"controllerType=0 controllerNumber=1 scan={scan_no}",
                                  mz=fm, it=fi, level=2, rt=rt + 0.002,
                                  prec=dict(mz=pmz, charge=2, intensity=1e4, scan_id=ms1_id)))

        with MzMLWriter(open(out / f"{name}.mzML", "wb"), close=True) as w:
            w.controlled_vocabularies()
            w.file_description(["MS1 spectrum", "MSn spectrum"])
            w.software_list([{"id": "psims-writer", "version": "1", "params": ["python-psims"]}])
            w.instrument_configuration_list([w.InstrumentConfiguration(
                id="IC1", component_list=w.ComponentList([
                    w.Source(1, ["electrospray ionization"]),
                    w.Analyzer(2, ["orbitrap"]),
                    w.Detector(3, ["inductive detector"])]),
                params=["Q Exactive"])])
            w.data_processing_list([w.DataProcessing(
                [dict(software_reference="psims-writer", order=0,
                      params=["Conversion to mzML"])], id="DP1")])
            with w.run(id=name, instrument_configuration="IC1"):
                with w.spectrum_list(count=len(scans)):
                    for s in scans:
                        params = [{"ms level": s["level"]},
                                  "MS1 spectrum" if s["level"] == 1 else "MSn spectrum",
                                  "centroid spectrum",
                                  {"total ion current": float(s["it"].sum())}]
                        kw = {}
                        if s["level"] == 2:
                            pr = s["prec"]
                            kw["precursor_information"] = {
                                "mz": pr["mz"], "intensity": pr["intensity"],
                                "charge": pr["charge"], "scan_id": pr["scan_id"],
                                "activation": ["beam-type collision-induced dissociation",
                                               {"collision energy": 28}],
                                "isolation_window": [0.7, pr["mz"], 0.7]}
                        w.write_spectrum(s["mz"], s["it"], id=s["id"], params=params,
                                         centroided=True,
                                         scan_start_time=float(s["rt"]), **kw)
    print(f"proteins={len(proteins)} peptides={len(peptides)} samples={len(samples)} -> {out}")


if __name__ == "__main__":
    main()
