#!/usr/bin/env python3
"""
List the RAW files of a PRIDE project and write a samplesheet with a
small, evenly spaced subset of them (keeps the demo light).

Output samplesheet columns: sample,file,group,size_mb
`file` is an https URL; Nextflow downloads it when the process needs it.
"""
import argparse
import csv
import json
import re
import sys
import urllib.request

ENDPOINTS = [
    "https://www.ebi.ac.uk/pride/ws/archive/v3/projects/{acc}/files/all",
    "https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{acc}/files?pageSize=1000&page=0",
    "https://www.ebi.ac.uk/pride/ws/archive/v2/files/byProject?accession={acc}",
]


def get_json(url):
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.loads(r.read().decode())


def list_files(acc):
    errors = []
    for ep in ENDPOINTS:
        url = ep.format(acc=acc)
        try:
            data = get_json(url)
        except Exception as e:  # noqa: BLE001
            errors.append(f"{url}: {e}")
            continue
        if isinstance(data, dict):
            data = data.get("_embedded", {}).get("files", data.get("files", []))
        if data:
            return data
        errors.append(f"{url}: empty")
    raise SystemExit("[fetch_pride] no PRIDE endpoint answered:\n  " + "\n  ".join(errors) +
                     "\nWrite the samplesheet by hand (see README).")


def https_location(f):
    for loc in f.get("publicFileLocations", []):
        v = loc.get("value", "")
        if v.startswith("ftp://"):
            return "https://" + v[len("ftp://"):]
        if v.startswith("http"):
            return v
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--accession", required=True)
    ap.add_argument("--max_files", type=int, default=6)
    ap.add_argument("--include_regex", default=None,
                    help="only keep RAW files whose name matches")
    ap.add_argument("--group_regex", default=None,
                    help="regex with one capture group applied to the file name to define groups")
    ap.add_argument("--max_file_mb", type=float, default=3000)
    ap.add_argument("--out", default="samplesheet.csv")
    ap.add_argument("--all_out", default="pride_files.tsv")
    a = ap.parse_args()

    files = list_files(a.accession)
    raws = []
    with open(a.all_out, "w") as fh:
        fh.write("file\tcategory\tsize_mb\turl\n")
        for f in files:
            name = f.get("fileName", "")
            cat = (f.get("fileCategory") or {}).get("value", "")
            size = (f.get("fileSizeBytes") or 0) / 1e6
            url = https_location(f)
            fh.write(f"{name}\t{cat}\t{size:.1f}\t{url}\n")
            if name.lower().endswith(".raw") and url:
                raws.append((name, size, url))

    raws.sort()
    if a.include_regex:
        raws = [r for r in raws if re.search(a.include_regex, r[0])]
    raws = [r for r in raws if r[1] <= a.max_file_mb]
    if not raws:
        raise SystemExit("[fetch_pride] no RAW files left after filtering")

    n = min(a.max_files, len(raws))
    step = len(raws) / n
    chosen = [raws[int(i * step)] for i in range(n)]

    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sample", "file", "group", "size_mb"])
        for name, size, url in chosen:
            stem = re.sub(r"\.raw$", "", name, flags=re.I)
            grp = "all"
            if a.group_regex:
                m = re.search(a.group_regex, name)
                grp = m.group(1) if m else "unassigned"
            w.writerow([stem, url, grp, f"{size:.1f}"])
    total = sum(c[1] for c in chosen)
    print(f"[fetch_pride] {len(raws)} RAW files available, selected {n} "
          f"({total:.0f} MB to download)", file=sys.stderr)


if __name__ == "__main__":
    main()
