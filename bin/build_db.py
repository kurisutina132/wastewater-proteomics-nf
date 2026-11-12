#!/usr/bin/env python3
"""
Build the search database from UniProt reference proteomes + cRAP contaminants.

Writes:
  db.fasta        target database (Sage generates decoys itself)
  taxmap.tsv      protein_id -> organism label, group   (used for peptide taxonomy)
  db_manifest.tsv one row per proteome with the number of entries downloaded

Fails loudly if any proteome query returns zero entries, so a wrong
proteome ID never silently produces a smaller database.
"""
import argparse
import csv
import gzip
import io
import sys
import time
import urllib.parse
import urllib.request

UNIPROT = "https://rest.uniprot.org/uniprotkb/stream"
CRAP = "https://ftp.thegpm.org/fasta/cRAP/crap.fasta"


def fetch(url, retries=4):
    for i in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=300) as r:
                data = r.read()
            if data[:2] == b"\x1f\x8b":
                data = gzip.decompress(data)
            return data.decode()
        except Exception as e:  # noqa: BLE001
            wait = 5 * (i + 1)
            print(f"[build_db] {e} -> retry in {wait}s", file=sys.stderr)
            time.sleep(wait)
    raise SystemExit(f"[build_db] could not download {url}")


def parse_fasta(text):
    header, seq = None, []
    for line in io.StringIO(text):
        line = line.rstrip()
        if line.startswith(">"):
            if header:
                yield header, "".join(seq)
            header, seq = line[1:], []
        elif line:
            seq.append(line)
    if header:
        yield header, "".join(seq)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--proteomes", required=True)
    ap.add_argument("--set", choices=["minimal", "full"], default="minimal")
    ap.add_argument("--contaminants", choices=["true", "false"], default="true")
    ap.add_argument("--out_prefix", default="db")
    a = ap.parse_args()

    rows = list(csv.DictReader(open(a.proteomes), delimiter="\t"))
    if a.set == "minimal":
        rows = [r for r in rows if r["set"] == "minimal"]

    seen = set()
    manifest = []
    with open(f"{a.out_prefix}.fasta", "w") as fa, open("taxmap.tsv", "w") as tm:
        tm.write("protein_id\torganism\tgroup\n")
        for r in rows:
            url = UNIPROT + "?" + urllib.parse.urlencode(
                {"query": r["query"], "format": "fasta", "compressed": "true"})
            print(f"[build_db] {r['label']}: {r['query']}", file=sys.stderr)
            n = 0
            for h, s in parse_fasta(fetch(url)):
                pid = h.split()[0]
                if pid in seen:
                    continue
                seen.add(pid)
                fa.write(f">{h}\n{s}\n")
                tm.write(f"{pid}\t{r['label']}\t{r['group']}\n")
                n += 1
            if n == 0:
                raise SystemExit(f"[build_db] 0 entries for {r['proteome']} ({r['label']}). "
                                 "Check the proteome ID in assets/proteomes.tsv")
            manifest.append((r["proteome"], r["label"], r["group"], n))

        if a.contaminants == "true":
            n = 0
            for h, s in parse_fasta(fetch(CRAP)):
                pid = "CONT_" + h.split()[0].replace("|", "_")
                if pid in seen:
                    continue
                seen.add(pid)
                fa.write(f">{pid} {h} OS=contaminant\n{s}\n")
                tm.write(f"{pid}\tcRAP contaminant\tContaminant\n")
                n += 1
            manifest.append(("cRAP", "cRAP contaminants", "Contaminant", n))

    with open("db_manifest.tsv", "w") as fh:
        fh.write("proteome\tlabel\tgroup\tn_entries\n")
        for m in manifest:
            fh.write("\t".join(map(str, m)) + "\n")
    print(f"[build_db] total entries: {len(seen)}", file=sys.stderr)


if __name__ == "__main__":
    main()
