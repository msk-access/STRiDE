import pysam
import numpy as np
import logging
from datetime import datetime
from scipy.stats import chi2_contingency, entropy, wasserstein_distance

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("MSIProfileGenerator")


class MSIProfileGenerator:
    def __init__(self, sites_file, tumor_bam_path, normal_bam_path, min_coverage=20, max_repeat_bins=100):
        self.sites_file = sites_file
        self.tumor_bam_path = tumor_bam_path
        self.normal_bam_path = normal_bam_path
        self.min_coverage = min_coverage
        self.max_repeat_bins = max_repeat_bins
        self.tumor_bam = pysam.AlignmentFile(self.tumor_bam_path, "rb")
        self.normal_bam = pysam.AlignmentFile(self.normal_bam_path, "rb")

    def parse_sites(self):
        with open(self.sites_file, 'r') as f:
            _ = next(f)  # header
            for line in f:
                yield line.strip().split('\t')

    def analyze_site(self, site):
        chrom, start, unit_len, _, times, *_, repeat_unit, left_flank, right_flank = site
        start = int(start)
        unit_len = int(unit_len)
        times = int(times)
        end = start + unit_len * times

        logger.info(f"Processing site {chrom}:{start}-{end}")

        if (unit_len == 1 and not (10 <= times <= 50)) or (unit_len > 1 and not (5 <= times <= 40)):
            logger.info(f"Skipped: repeat count filter (unit_len={unit_len}, times={times})")
            return None

        if chrom not in self.tumor_bam.references or chrom not in self.normal_bam.references:
            logger.info(f"Skipped: chromosome {chrom} not in BAM references")
            return None

        tumor_cov = self.tumor_bam.count(chrom, start, end)
        normal_cov = self.normal_bam.count(chrom, start, end)
        if tumor_cov < self.min_coverage or normal_cov < self.min_coverage:
            logger.info(f"Skipped: low coverage (tumor={tumor_cov}, normal={normal_cov})")
            return None

        try:
            tumor_freq, tumor_mapq, tumor_bq, tumor_insert_all, tumor_insert_ref, tumor_insert_alt = self.get_repeat_frequencies(
                self.tumor_bam.fetch(chrom, start, end), repeat_unit, left_flank, right_flank, start, end, times
            )
            normal_freq, normal_mapq, normal_bq, normal_insert_all, normal_insert_ref, normal_insert_alt = self.get_repeat_frequencies(
                self.normal_bam.fetch(chrom, start, end), repeat_unit, left_flank, right_flank, start, end, times
            )

            if sum(tumor_freq) < self.min_coverage or sum(normal_freq) < self.min_coverage:
                logger.info(f"Skipped: summed frequency below min_coverage")
                return None

            idx = np.union1d(np.nonzero(tumor_freq)[0], np.nonzero(normal_freq)[0])
            contingency = [tumor_freq[idx], normal_freq[idx]]
            _, p_value, _, _ = chi2_contingency(contingency, correction=False)

            def shannon_entropy(freqs):
                f = np.array(freqs, dtype=float)
                f = f[f > 0] / np.sum(f)
                return entropy(f, base=2) if len(f) > 0 else 0

            def count_alleles_above_thresh(freqs, threshold=5):
                return np.sum(np.array(freqs) >= threshold)

            tumor_entropy_val = shannon_entropy(tumor_freq)
            normal_entropy_val = shannon_entropy(normal_freq)

            tumor_total = np.sum(tumor_freq)
            normal_total = np.sum(normal_freq)

            norm_tumor = tumor_freq / tumor_total if tumor_total > 0 else np.zeros_like(tumor_freq)
            norm_normal = normal_freq / normal_total if normal_total > 0 else np.zeros_like(normal_freq)

            l1_dist = np.sum(np.abs(norm_tumor - norm_normal))
            l2_dist = np.sqrt(np.sum((norm_tumor - norm_normal) ** 2))
            wass_dist = wasserstein_distance(norm_tumor, norm_normal)

            allele_thresholds = [1, 5, 10, 15, 20, 25, 30]
            allele_diffs = {
                f"n_alleles_diff_{thresh}": count_alleles_above_thresh(tumor_freq, threshold=thresh) -
                                          count_alleles_above_thresh(normal_freq, threshold=thresh)
                for thresh in allele_thresholds
            }

            norm_allele_diffs = {
                f"n_alleles_diff_norm_{int(thresh * 100)}": np.sum(norm_tumor >= thresh) - np.sum(norm_normal >= thresh)
                for thresh in [0.01, 0.02, 0.03, 0.04, 0.05, 0.06]
            }

            return {
                "chrom": chrom,
                "start": start,
                "unit_len": unit_len,
                "repeat_unit": repeat_unit,
                "repeat_count": times,
                "left_flank": left_flank,
                "right_flank": right_flank,
                "tumor_freqs": tumor_freq.tolist(),
                "normal_freqs": normal_freq.tolist(),
                "tumor_norm_freqs": norm_tumor.tolist(),
                "normal_norm_freqs": norm_normal.tolist(),
                "tumor_total_coverage": tumor_total,
                "normal_total_coverage": normal_total,
                "tumor_mapq_mean": np.mean(tumor_mapq) if tumor_mapq else 0,
                "normal_mapq_mean": np.mean(normal_mapq) if normal_mapq else 0,
                "tumor_bq_mean": np.mean(tumor_bq) if tumor_bq else 0,
                "normal_bq_mean": np.mean(normal_bq) if normal_bq else 0,
                "tumor_insert_mean_all": np.mean(tumor_insert_all) if tumor_insert_all else 0,
                "normal_insert_mean_all": np.mean(normal_insert_all) if normal_insert_all else 0,
                "tumor_insert_mean_ref": np.mean(tumor_insert_ref) if tumor_insert_ref else 0,
                "normal_insert_mean_ref": np.mean(normal_insert_ref) if normal_insert_ref else 0,
                "tumor_insert_mean_alt": np.mean(tumor_insert_alt) if tumor_insert_alt else 0,
                "normal_insert_mean_alt": np.mean(normal_insert_alt) if normal_insert_alt else 0,
                "tumor_entropy": tumor_entropy_val,
                "normal_entropy": normal_entropy_val,
                "entropy_diff": tumor_entropy_val - normal_entropy_val,
                "l1_distance": l1_dist,
                "l2_distance": l2_dist,
                "wasserstein_distance": wass_dist,
                "p_value": p_value,
                **allele_diffs,
                **norm_allele_diffs
            }
        except Exception as e:
            logger.error(f"Error processing site {site}: {e}")
            return None

    def get_repeat_frequencies(self, reads, repeat_unit, left_flank, right_flank, start, end, ref_repeat_count):
        freq = np.zeros(self.max_repeat_bins, dtype=int)
        map_quals = []
        base_quals = []
        insert_sizes_all = []
        insert_sizes_ref = []
        insert_sizes_alt = []

        for read in reads:
            if read.is_unmapped:
                continue
            count = self.count_repeats(read.query_sequence, repeat_unit, left_flank, right_flank)
            if count > 0 and count - 1 < self.max_repeat_bins:
                freq[count - 1] += 1
                if 10 < abs(read.template_length) < 1000:
                    insert_size = abs(read.template_length)
                    insert_sizes_all.append(insert_size)
                    if count == ref_repeat_count:
                        insert_sizes_ref.append(insert_size)
                    else:
                        insert_sizes_alt.append(insert_size)
                map_quals.append(read.mapping_quality)
                for qpos, rpos in read.get_aligned_pairs(with_seq=False):
                    if rpos is not None and start <= rpos < end and qpos is not None:
                        base_quals.append(read.query_qualities[qpos])

        return freq, map_quals, base_quals, insert_sizes_all, insert_sizes_ref, insert_sizes_alt

    def count_repeats(self, seq, repeat_unit, left, right):
        start = 0
        while True:
            start = seq.find(left, start)
            if start == -1:
                break
            pos = start + len(left)
            count = 0
            while seq.startswith(repeat_unit, pos):
                count += 1
                pos += len(repeat_unit)
                if seq.startswith(right, pos):
                    if (len(repeat_unit) == 1 and count >= 5) or (len(repeat_unit) > 1 and count >= 3):
                        return count
                    else:
                        return 0
            start += 1
        return 0

    def run(self, output_tsv):
        logger.info("Starting MSI analysis...")
        start_time = datetime.now()

        results = [self.analyze_site(site) for site in self.parse_sites()]

        header = [
            "chrom", "start", "unit_len", "repeat_unit", "repeat_count", "left_flank", "right_flank",
            "normal_freqs", "tumor_freqs",
            "normal_norm_freqs", "tumor_norm_freqs",
            "normal_total_coverage", "tumor_total_coverage",
            "normal_mapq_mean", "tumor_mapq_mean",
            "normal_bq_mean", "tumor_bq_mean",
            "normal_insert_mean_all", "tumor_insert_mean_all",
            "normal_insert_mean_ref", "tumor_insert_mean_ref",
            "normal_insert_mean_alt", "tumor_insert_mean_alt",
            "normal_entropy", "tumor_entropy", "entropy_diff",
            "l1_distance", "l2_distance", "wasserstein_distance", "p_value"
        ]
        header += [f"n_alleles_diff_{thresh}" for thresh in [1, 5, 10, 15, 20, 25, 30]]
        header += [f"n_alleles_diff_norm_{i}" for i in [1, 2, 3, 4, 5, 6]]

        with open(output_tsv, 'w') as out_file:
            out_file.write("\t".join(header) + "\n")
            for result in filter(None, results):
                line = [
                    result["chrom"], str(result["start"]), str(result["unit_len"]), result["repeat_unit"],
                    str(result["repeat_count"]), result["left_flank"], result["right_flank"],
                    " ".join(map(str, result["normal_freqs"])),
                    " ".join(map(str, result["tumor_freqs"])),
                    " ".join(f"{v:.5f}" for v in result["normal_norm_freqs"]),
                    " ".join(f"{v:.5f}" for v in result["tumor_norm_freqs"]),
                    str(result["normal_total_coverage"]), str(result["tumor_total_coverage"]),
                    f"{result['normal_mapq_mean']:.2f}", f"{result['tumor_mapq_mean']:.2f}",
                    f"{result['normal_bq_mean']:.2f}", f"{result['tumor_bq_mean']:.2f}",
                    f"{result['normal_insert_mean_all']:.2f}", f"{result['tumor_insert_mean_all']:.2f}",
                    f"{result['normal_insert_mean_ref']:.2f}", f"{result['tumor_insert_mean_ref']:.2f}",
                    f"{result['normal_insert_mean_alt']:.2f}", f"{result['tumor_insert_mean_alt']:.2f}",
                    f"{result['normal_entropy']:.4f}", f"{result['tumor_entropy']:.4f}", f"{result['entropy_diff']:.4f}",
                    f"{result['l1_distance']:.4f}", f"{result['l2_distance']:.4f}", f"{result['wasserstein_distance']:.4f}",
                    f"{result['p_value']:.6g}"
                ]
                line += [str(result[f"n_alleles_diff_{thresh}"]) for thresh in [1, 5, 10, 15, 20, 25, 30]]
                line += [str(result[f"n_alleles_diff_norm_{i}"]) for i in [1, 2, 3, 4, 5, 6]]
                assert len(header) == len(line), "Mismatch between header and output columns"
                out_file.write("\t".join(line) + "\n")

        end_time = datetime.now()
        logger.info(f"MSI analysis completed in {end_time - start_time}.")
