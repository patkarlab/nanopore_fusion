#!/usr/bin/env python3

from statistics import mean, median
import gzip
import matplotlib
import matplotlib.pyplot as plt
import argparse
import os

parser = argparse.ArgumentParser(description="Extract read lengths and plot histogram from FASTQ. As this script is specifically for targeted data, the read lengths evaluated are from 0-2000 basepairs. This can be changed by changing the bins values")
parser.add_argument("-s", "--sample_name", required=True, help="Sample name")
parser.add_argument("-f", "--fastq", required=True, help="Path to FASTQ file")
parser.add_argument("-o", "--outfile", required=True, help="Output image path")
args = parser.parse_args()

def read_lengths_and_gc_from_fastq(fastq_file):
    lengths, gc_counts, bases_counts = [], [], []
    opener = gzip.open if fastq_file.endswith(".gz") else open
    with opener(fastq_file, "rt") as f:
        for i, line in enumerate(f):
            if i % 4 == 1:
                seq = line.strip()
                if seq:  # avoid empty lines
                    lengths.append(len(seq))
                    gc_counts.append(seq.count("G") + seq.count("C") + seq.count("g") + seq.count("c"))
                    bases_counts.append(len(seq))
    return lengths, gc_counts, bases_counts

def calculate_n50(lengths):
    lengths_sorted = sorted(lengths, reverse=True)
    total_bases = sum(lengths_sorted)
    half_total = total_bases / 2
    cum_sum = 0
    for l in lengths_sorted:
        cum_sum += l
        if cum_sum >= half_total:
            return l
    return 0

def summarize_read_lengths():
    print(f"Reading: {args.fastq}")
    lengths, gc_counts, bases_counts = read_lengths_and_gc_from_fastq(args.fastq)

    if not lengths:
        print(" No reads found — file may be empty or malformed.")
        plt.figure()
        plt.text(0.5, 0.5, "No reads found", ha='center', va='center', fontsize=12)
        plt.savefig(args.outfile)
        print(f"Empty plot saved to {args.outfile}")
        return

    print(f"Total reads: {len(lengths)}")
    mean_len, median_len, n50_len = mean(lengths), median(lengths), calculate_n50(lengths)
    avg_gc = (sum(gc_counts) / sum(bases_counts)) * 100 if bases_counts else 0

    plt.figure(figsize=(12,7))
    bins = [0,100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000]
    plt.hist(lengths, bins=bins, color="skyblue", edgecolor="black")
    plt.title(f"{args.sample_name} Read Length Distribution")
    plt.xlabel("Read Length (bases)")
    plt.ylabel("Read Count")
    #plt.yscale("log")
    plt.grid(axis="y", linestyle="--", alpha=0.7)

    stats_text = (
        f"Mean: {mean_len:.2f}\nMedian: {median_len}\nN50: {n50_len}\n"
        f"Avg GC%: {avg_gc:.2f}\nReads: {len(lengths)}"
    )
    plt.text(0.95, 0.95, stats_text, transform=plt.gca().transAxes,
             fontsize=10, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.7))
    
    plt.tight_layout()
    plt.savefig(args.outfile)
    print(f" Plot saved to {os.path.abspath(args.outfile)}")

if __name__ == "__main__":
    summarize_read_lengths()
