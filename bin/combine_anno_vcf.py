#!/usr/bin/env python3

import pandas as pd
import argparse
import matplotlib.pyplot as plt

# Argument parsing
parser = argparse.ArgumentParser(
    description="Merge annotated VCF tables into a multi-sheet Excel file with common variants"
)
parser.add_argument("-s", "--sample_id", required=True,
                    help="Sample ID for labeling outputs")
parser.add_argument("-l", "--longshot", required=True,
                    help="Longshot annotated VCF (tabular) file")
parser.add_argument("-n", "--nanocaller", required=True,
                    help="NanoCaller annotated VCF (tabular) file")
parser.add_argument("-c", "--clairsTO", required=True,
                    help="ClairsTO annotated VCF (tabular) file")
parser.add_argument("-r", "--clair3_rna", required=True,
                    help="Clair3-RNA annotated VCF (tabular) file")
parser.add_argument("-g", "--longcallr", required=True,
                    help="LongcallR annotated VCF (tabular) file")
parser.add_argument("-o", "--output", required=True,
                    help="Output Excel file")
args = parser.parse_args()

# Columns to extract
required_cols = [
    "Chr",
    "Start",
    "End",
    "Ref",
    "Alt",
    "AAChange.refGene",
    "cosmic70",
    "CLNSIG",
    "AF_popmax"
]

key_cols = ["Chr", "Start", "End", "Ref", "Alt"]

# VAF extraction
def extract_vaf(otherinfo):
    try:
        fields = otherinfo.split("\t")
        format_keys = fields[-2].split(":")
        sample_vals = fields[-1].split(":")
        fmt = dict(zip(format_keys, sample_vals))

        if "AF" in fmt:
            return round(float(fmt["AF"]) * 100, 2)

        if "VF" in fmt:
            return round(float(fmt["VF"]) * 100, 2)

        if "AD" in fmt:
            ad_vals = fmt["AD"].split(",")
            if len(ad_vals) >= 2:
                ref = float(ad_vals[0])
                alt = float(ad_vals[1])
                if (ref + alt) > 0:
                    return round((alt / (ref + alt)) * 100, 2)

    except Exception:
        return None

    return None

# AD extraction
def extract_ad_counts(otherinfo):
    try:
        fields = otherinfo.split("\t")
        format_keys = fields[-2].split(":")
        sample_vals = fields[-1].split(":")
        fmt = dict(zip(format_keys, sample_vals))

        if "AD" in fmt:
            ad_vals = fmt["AD"].split(",")
            if len(ad_vals) >= 2:
                return ad_vals[0], ad_vals[1]

    except Exception:
        pass

    return None, None

# Process file
def process_file(file):
    df = pd.read_csv(file, dtype=str, low_memory=False)
    df.columns = df.columns.str.strip()

    col_map = {
        "#Chr": "Chr",
        "Chrom": "Chr",
        "CLIN_SIG": "CLNSIG",
        "ClinSig": "CLNSIG",
        "COSMIC": "cosmic70"
    }
    df.rename(columns=col_map, inplace=True)

    missing = [col for col in required_cols if col not in df.columns]
    if missing:
        print(f"Warning: Missing columns in {file}: {missing}")

    # VAF + AD
    if "Otherinfo" in df.columns:
        df["VAF_%"] = df["Otherinfo"].apply(extract_vaf)

        counts = df["Otherinfo"].apply(extract_ad_counts)
        df["Ref_count"] = [c[0] for c in counts]
        df["Alt_count"] = [c[1] for c in counts]

    # Column order
    cols_present = [
        "Chr", "Start", "End", "Ref", "Alt",
        "VAF_%", "Ref_count", "Alt_count",
        "AAChange.refGene", "cosmic70", "CLNSIG", "AF_popmax"
    ]

    cols_present = [col for col in cols_present if col in df.columns]
    df = df[cols_present]

    df = df.drop_duplicates(subset=key_cols)

    return df

# -----------------------------
# Natural chromosome sorting function
# -----------------------------
def sort_df(df):
    def chr_sort_key(chrom):
        chrom = str(chrom).replace("chr", "").upper()

        if chrom.isdigit():
            return int(chrom)
        elif chrom == "X":
            return 23
        elif chrom == "Y":
            return 24
        elif chrom in ["M", "MT"]:
            return 25
        else:
            return 100  # unknown contigs

    df["chr_order"] = df["Chr"].apply(chr_sort_key)

    df = df.sort_values(
        by=["chr_order", "Start"],
        ascending=[True, True]
    )

    return df.drop(columns=["chr_order"])

# Process inputs
longshot_df   = process_file(args.longshot)
nanocaller_df = process_file(args.nanocaller)
clairsTO_df  = process_file(args.clairsTO)
clair3rna_df  = process_file(args.clair3_rna)
longcallr_df  = process_file(args.longcallr)

# Find common variants (inner join) - ALL 5 callers
common_df = clair3rna_df.merge(nanocaller_df[key_cols], on=key_cols, how="inner") \
                      .merge(clairsTO_df[key_cols], on=key_cols, how="inner") \
                      .merge(longshot_df[key_cols], on=key_cols, how="inner") \
                      .merge(longcallr_df[key_cols], on=key_cols, how="inner")

# Keep only one copy of annotation columns (from first dataset)
common_df = common_df[key_cols + [
    col for col in required_cols if col not in key_cols
] + [col for col in ["VAF_%", "Ref_count", "Alt_count"] if col in common_df.columns]]

# -----------------------------
# Build CONSENSUS table
# -----------------------------
all_variants = pd.concat([
    longshot_df.assign(Caller="Longshot"),
    nanocaller_df.assign(Caller="NanoCaller"),
    clairsTO_df.assign(Caller="ClairsTO"),
    clair3rna_df.assign(Caller="Clair3RNA"),
    longcallr_df.assign(Caller="LongcallR")
], ignore_index=True)

# Group by variant
consensus_df = all_variants.groupby(key_cols).agg({
    "Caller": lambda x: ",".join(sorted(set(x))),
    "VAF_%": lambda x: round(pd.to_numeric(x, errors="coerce").mean(), 2)
}).reset_index()

# Add support count
consensus_df["Support"] = consensus_df["Caller"].apply(lambda x: len(x.split(",")))

# Reorder columns
consensus_df = consensus_df[key_cols + ["Support", "Caller", "VAF_%"]]
consensus_df = consensus_df.rename(columns={"VAF_%": "Mean_VAF%"})

# -----------------------------
# Support distribution (1–5 callers)
# -----------------------------
support_dist = consensus_df["Support"].value_counts().sort_index()

# Ensure all categories (1–5) are present
for i in range(1, 6):
    if i not in support_dist:
        support_dist.loc[i] = 0

support_dist = support_dist.sort_index().reset_index()
support_dist.columns = ["Num_Callers", "Variant_Count"]

#bar chart
plt.figure()
plt.bar(support_dist["Num_Callers"], support_dist["Variant_Count"])
plt.xlabel("Number of Callers Supporting Variant")
plt.ylabel("Number of Variants")
plt.title(f"{args.sample_id} - Variant Support Distribution")
plt.xticks([1, 2, 3, 4, 5])

plt.savefig(f"{args.sample_id}_support_distribution.png")
plt.close()

# -----------------------------
# Apply sorting to all dataframes
# -----------------------------
longshot_df   = sort_df(longshot_df)
nanocaller_df = sort_df(nanocaller_df)
clairsTO_df   = sort_df(clairsTO_df)
clair3rna_df  = sort_df(clair3rna_df)
longcallr_df  = sort_df(longcallr_df)
common_df     = sort_df(common_df)
consensus_df  = sort_df(consensus_df)

# Write Excel

with pd.ExcelWriter(args.output) as writer:
    consensus_df.to_excel(writer, sheet_name="Consensus", index=False)
    common_df.to_excel(writer, sheet_name="Common", index=False)
    support_dist.to_excel(writer, sheet_name="Support_Distribution", index=False)
    clair3rna_df.to_excel(writer, sheet_name="Clair3RNA", index=False)
    longcallr_df.to_excel(writer, sheet_name="LongcallR", index=False)
    longshot_df.to_excel(writer, sheet_name="Longshot", index=False)
    nanocaller_df.to_excel(writer, sheet_name="NanoCaller", index=False)
    clairsTO_df.to_excel(writer, sheet_name="ClairsTO", index=False)


print(f"Excel file created with Common variants: {args.output}")
