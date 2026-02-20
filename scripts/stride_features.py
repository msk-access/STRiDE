#!/usr/bin/env python3
import argparse
import os

from stride.pipeline import run_feature_generation
from stride.utils import safe_name, strip_ext

def main():
    p = argparse.ArgumentParser(description="Generate STRiDE MSI feature TSV from tumor/normal BAMs.")
    p.add_argument("--site-list", required=True, help="Microsatellite site list TSV")
    p.add_argument("--tumor-bam", required=True, help="Tumor BAM")
    p.add_argument("--normal-bam", required=True, help="Normal BAM")
    p.add_argument("--output-tsv", default=None, help="Optional explicit output TSV path.")
    p.add_argument("--out-dir", default="out", help="Base output directory if --output-tsv is not set.")
    p.add_argument("--sample-id", default=None, help="Sample ID (used for output filename).")
    p.add_argument("--min-coverage", type=int, default=20)
    p.add_argument("--max-repeat-bins", type=int, default=100)
    args = p.parse_args()

    if args.output_tsv:
        out_tsv = args.output_tsv
    else:
        sid = args.sample_id.strip() if args.sample_id else strip_ext(args.tumor_bam)
        features_dir = os.path.join(args.out_dir, "features")
        os.makedirs(features_dir, exist_ok=True)
        out_tsv = os.path.join(features_dir, f"msi_features_{safe_name(sid)}.tsv")

    run_feature_generation(
        sites_file=args.site_list,
        tumor_bam=args.tumor_bam,
        normal_bam=args.normal_bam,
        out_features_tsv=out_tsv,
        min_coverage=args.min_coverage,
        max_repeat_bins=args.max_repeat_bins
    )
    print(out_tsv)

if __name__ == "__main__":
    main()
