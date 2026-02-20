#!/usr/bin/env python3
import argparse
from stride.pipeline import run_end_to_end_single, run_end_to_end_batch

def main():
    p = argparse.ArgumentParser(description="End-to-end STRiDE MSI pipeline: BAMs -> feature TSV -> per-sample MSI predictions.")
    p.add_argument("--site-list", required=True, help="Microsatellite site list TSV")
    p.add_argument("--model-joblib", required=True, help="Trained model .joblib")
    p.add_argument("--out-dir", required=True)

    # batch mode
    p.add_argument("--samples-list", default=None,
                   help="CSV/TSV with columns sample_id,tumor_bam,normal_bam. If provided, runs batch mode.")

    # single mode
    p.add_argument("--sample-id", default=None)
    p.add_argument("--tumor-bam", default=None)
    p.add_argument("--normal-bam", default=None)

    p.add_argument("--min-coverage", type=int, default=20)
    p.add_argument("--max-repeat-bins", type=int, default=100)
    p.add_argument(
        "--delete-features",
        action="store_true",
        help="Delete generated feature TSVs after prediction (default: keep)."
    )
    args = p.parse_args()

    if args.samples_list:
        results = run_end_to_end_batch(
            samples_list_path=args.samples_list,
            sites_file=args.site_list,
            model_joblib=args.model_joblib,
            out_dir=args.out_dir,
            min_coverage=args.min_coverage,
            max_repeat_bins=args.max_repeat_bins,
            keep_features=not args.delete_features
        )
        print(f"Completed batch for {len(results)} samples.")
        for r in results:
            print(f"{r['sample_id']}\t{r['prediction_txt']}")
        return

    # single mode requires tumor+normal BAM
    if not args.tumor_bam or not args.normal_bam:
        raise SystemExit("Provide either --samples-list (batch) OR (--tumor-bam and --normal-bam) for single mode.")

    res = run_end_to_end_single(
        sites_file=args.site_list,
        tumor_bam=args.tumor_bam,
        normal_bam=args.normal_bam,
        model_joblib=args.model_joblib,
        out_dir=args.out_dir,
        sample_id=args.sample_id,
        min_coverage=args.min_coverage,
        max_repeat_bins=args.max_repeat_bins,
        keep_features=not args.delete_features
    )
    print(f"Completed sample {res['sample_id']}")
    print(res["prediction_txt"])

if __name__ == "__main__":
    main()
