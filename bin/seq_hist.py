#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt
from statistics import mean, median

#this script limits the length of reads in histogram to 2201. Can be changed by using the bins parameter.
parser = argparse.ArgumentParser(
    description="Plot read length histogram and compute stats from seqkit output"
)
parser.add_argument("-i", "--input", required=True, help="seqkit stats-like file")
parser.add_argument("-o", "--outfile", required=True, help="Output image path")
parser.add_argument("-s", "--sample_name", required=True, help="Sample name")
parser.add_argument(
    "-t", "--threshold", type=int, default=500,
    help="Length threshold for 'total reads > X bp' (default: 500)"
)
args = parser.parse_args()


def read_lengths_from_seqkit(file):
    lengths = []
    with open(file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split()
            try:
                length = int(parts[1])
                lengths.append(length)
            except:
                continue
    return lengths


def calculate_nxx(lengths, fraction):
    """Generic NXX (e.g., N50=0.5, N75=0.75)"""
    lengths_sorted = sorted(lengths, reverse=True)
    total = sum(lengths_sorted)
    cutoff = total * fraction

    cumsum = 0
    for l in lengths_sorted:
        cumsum += l
        if cumsum >= cutoff:
            return l
    return 0


def main():
    lengths = read_lengths_from_seqkit(args.input)

    if not lengths:
        print("No valid reads found.")
        return

    total_reads = len(lengths)
    total_bases = sum(lengths)
    mean_len = mean(lengths)
    median_len = median(lengths)

    min_len = min(lengths)
    max_len = max(lengths)

    n50_len = calculate_nxx(lengths, 0.5)
    n75_len = calculate_nxx(lengths, 0.75)

    threshold = args.threshold
    bases_above_thresh = sum(l for l in lengths if l > threshold)
    reads_above_thresh = sum(1 for l in lengths if l > threshold)

    # ---- PRINT ----
    print(f"Total reads: {total_reads}")
    print(f"Total bases: {total_bases}")
    print(f"Mean length: {mean_len:.2f}")
    print(f"Median length: {median_len}")
    print(f"Min length: {min_len}")
    print(f"Max length: {max_len}")
    print(f"N50: {n50_len}")
    print(f"N75: {n75_len}")
    print(f"Reads > {threshold} bp: {reads_above_thresh}")
    print(f"Bases in reads > {threshold} bp: {bases_above_thresh}")

    # ---- PLOT ----
    bins = list(range(0, 2201, 100))

    plt.figure(figsize=(12, 7))
    plt.hist(lengths, bins=bins, edgecolor="black")
    plt.title(f"{args.sample_name} Read Length Distribution")
    plt.xlabel("Read Length (bp)")
    plt.ylabel("Read Count")
    plt.grid(axis="y", linestyle="--", alpha=0.7)

    stats_text = (
        f"Total Reads: {total_reads}\n"
        f"Total Bases: {total_bases}\n"
        f"Mean length: {mean_len:.2f}\n"
        f"Median length: {median_len}\n"
        f"Min length: {min_len} | Max length: {max_len}\n"
        f"N50: {n50_len}\n"
        f"N75: {n75_len}\n"
        f"Reads > {threshold}bp: {reads_above_thresh}\n"
		f"Bases in reads > {threshold}bp: {bases_above_thresh}"
    )

    plt.text(
        0.95, 0.95, stats_text,
        transform=plt.gca().transAxes,
        ha="right", va="top",
        bbox=dict(facecolor="white", alpha=0.7)
    )

    plt.tight_layout()
    plt.savefig(args.outfile)

    print(f"Plot saved to {args.outfile}")


if __name__ == "__main__":
    main()

