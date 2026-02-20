#!/usr/bin/env python3
import argparse
import glob
import os

from stride.predictor import predict_from_feature_tsvs, write_one_output_per_sample

def main():
    p = argparse.ArgumentParser(description="Evaluate MSI/MSS from STRiDE feature TSV file(s).")
    p.add_argument("--model-joblib", required=True, help="Trained model .joblib")
    p.add_argument("--features-dir", default=None, help="Directory containing feature TSVs (searched recursively).")
    p.add_argument("--feature-files", nargs="*", default=None, help="One or more feature TSV files.")
    p.add_argument("--list-file", default=None, help="Text file with one feature TSV path per line.")
    p.add_argument("--out-dir", required=True, help="Directory where per-sample *_msi.txt will be written.")
    # Backward-compatible aliases
    p.add_argument("--samples-dir", default=None, help=argparse.SUPPRESS)
    p.add_argument("--sample-files", nargs="*", default=None, help=argparse.SUPPRESS)
    args = p.parse_args()

    features_dir = args.features_dir or args.samples_dir
    feature_files = args.feature_files or args.sample_files

    # Build feature TSV list from inputs.
    feature_tsvs = []
    if features_dir:
        feature_tsvs.extend(glob.glob(os.path.join(features_dir, "**", "*.tsv"), recursive=True))
    if feature_files:
        feature_tsvs.extend(feature_files)
    if args.list_file:
        with open(args.list_file, "r") as f:
            for ln in f:
                ln = ln.strip()
                if ln:
                    feature_tsvs.append(ln)

    df_preds = predict_from_feature_tsvs(args.model_joblib, feature_tsvs)
    paths = write_one_output_per_sample(df_preds, args.out_dir)
    print(f"Wrote {len(paths)} prediction files to {args.out_dir}")

if __name__ == "__main__":
    main()
