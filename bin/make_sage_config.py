#!/usr/bin/env python3
"""Write a Sage JSON config from pipeline parameters."""
import argparse
import json


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", default="sage_config.json")
    ap.add_argument("--missed_cleavages", type=int, default=2)
    ap.add_argument("--precursor_ppm", type=float, default=10.0)
    ap.add_argument("--fragment_ppm", type=float, default=20.0)
    ap.add_argument("--carbamidomethyl", choices=["true", "false"], default="true")
    ap.add_argument("--oxidation", choices=["true", "false"], default="true")
    ap.add_argument("--min_charge", type=int, default=2)
    ap.add_argument("--max_charge", type=int, default=4)
    ap.add_argument("--min_matched_peaks", type=int, default=4)
    ap.add_argument("--lfq", choices=["true", "false"], default="true")
    a = ap.parse_args()

    static_mods = {"C": 57.021464} if a.carbamidomethyl == "true" else {}
    variable_mods = {"M": [15.994915]} if a.oxidation == "true" else {}

    cfg = {
        "database": {
            "bucket_size": 8192,
            "enzyme": {"missed_cleavages": a.missed_cleavages, "min_len": 7,
                       "max_len": 30, "cleave_at": "KR", "restrict": "P"},
            "fragment_min_mz": 150.0,
            "fragment_max_mz": 2000.0,
            "peptide_min_mass": 500.0,
            "peptide_max_mass": 5000.0,
            "ion_kinds": ["b", "y"],
            "min_ion_index": 2,
            "static_mods": static_mods,
            "variable_mods": variable_mods,
            "max_variable_mods": 2,
            "decoy_tag": "rev_",
            "generate_decoys": True,
            "fasta": a.fasta,
        },
        "precursor_tol": {"ppm": [-a.precursor_ppm, a.precursor_ppm]},
        "fragment_tol": {"ppm": [-a.fragment_ppm, a.fragment_ppm]},
        "precursor_charge": [a.min_charge, a.max_charge],
        "isotope_errors": [-1, 3],
        "deisotope": False,
        "chimera": False,
        "wide_window": False,
        "predict_rt": True,
        "min_peaks": 15,
        "max_peaks": 150,
        "min_matched_peaks": a.min_matched_peaks,
        "report_psms": 1,
        "quant": {
            "lfq": a.lfq == "true",
            "lfq_settings": {"peak_scoring": "Hybrid", "integration": "Sum",
                             "spectral_angle": 0.7, "ppm_tolerance": 5.0,
                             "combine_charge_states": True},
        },
    }
    with open(a.out, "w") as fh:
        json.dump(cfg, fh, indent=2)


if __name__ == "__main__":
    main()
