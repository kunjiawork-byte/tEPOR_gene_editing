#!/usr/bin/env python3
# Strict-set spectrum only; not the later 51-variant FINAL spectrum.
# Original computation is retained, including assumptions of nonempty SNV sets.

import csv
import os
from pathlib import Path
from collections import Counter

OUT = Path(os.environ.get("OUT", "/path/to/WESsample/mutect2_reanalysis"))

SAMPLES = ["3XBE", "tEPORCas9"]

# Complement bases
COMP = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C"
}

def normalize_snv(ref, alt):
    """
    Convert SNV to pyrimidine-centered representation.

    Examples:
      C>T -> C>T
      G>A -> C>T
      A>G -> T>C
      T>C -> T>C
    """

    ref = ref.upper()
    alt = alt.upper()

    if ref in ("C", "T"):
        return f"{ref}>{alt}"

    if ref in ("A", "G"):
        return f"{COMP[ref]}>{COMP[alt]}"

    return "OTHER"


all_snv_rows = []

raw_counts = {}
norm_counts = {}

for sample in SAMPLES:

    infile = OUT / sample / f"{sample}.PASS.DP30.ALT5.tsv"

    sample_rows = []

    with open(infile) as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:

            if row["TYPE"] != "SNV":
                continue

            ref = row["REF"].upper()
            alt = row["ALT"].upper()

            raw = f"{ref}>{alt}"
            normalized = normalize_snv(ref, alt)

            row["Sample"] = sample
            row["RAW_SUBSTITUTION"] = raw
            row["NORMALIZED_SUBSTITUTION"] = normalized

            # CBE-type transition:
            # C>T or complementary G>A
            row["CBE_CT_GA"] = (
                "YES"
                if raw in ("C>T", "G>A")
                else "NO"
            )

            sample_rows.append(row)
            all_snv_rows.append(row)

    raw_counts[sample] = Counter(
        r["RAW_SUBSTITUTION"]
        for r in sample_rows
    )

    norm_counts[sample] = Counter(
        r["NORMALIZED_SUBSTITUTION"]
        for r in sample_rows
    )


# --------------------------------------------------
# Write per-variant SNV table
# --------------------------------------------------

outfile = OUT / "SNV_substitution_classification.tsv"

fieldnames = [
    "Sample",
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "DP",
    "REF_READS",
    "ALT_READS",
    "VAF",
    "RAW_SUBSTITUTION",
    "NORMALIZED_SUBSTITUTION",
    "CBE_CT_GA"
]

with open(outfile, "w", newline="") as f:

    writer = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=fieldnames,
        extrasaction="ignore"
    )

    writer.writeheader()

    for row in all_snv_rows:
        writer.writerow(row)


# --------------------------------------------------
# 12-class raw spectrum
# --------------------------------------------------

RAW_CLASSES = [
    "A>C", "A>G", "A>T",
    "C>A", "C>G", "C>T",
    "G>A", "G>C", "G>T",
    "T>A", "T>C", "T>G"
]

raw_out = OUT / "SNV_12class_spectrum.tsv"

with open(raw_out, "w", newline="") as f:

    writer = csv.writer(f, delimiter="\t")

    writer.writerow([
        "Substitution",
        "3XBE_count",
        "3XBE_percent",
        "tEPORCas9_count",
        "tEPORCas9_percent"
    ])

    for mut in RAW_CLASSES:

        n3 = raw_counts["3XBE"][mut]
        nt = raw_counts["tEPORCas9"][mut]

        total3 = sum(raw_counts["3XBE"].values())
        totalt = sum(raw_counts["tEPORCas9"].values())

        writer.writerow([
            mut,
            n3,
            100 * n3 / total3 if total3 else 0,
            nt,
            100 * nt / totalt if totalt else 0
        ])


# --------------------------------------------------
# 6-class normalized spectrum
# --------------------------------------------------

NORM_CLASSES = [
    "C>A",
    "C>G",
    "C>T",
    "T>A",
    "T>C",
    "T>G"
]

norm_out = OUT / "SNV_6class_spectrum.tsv"

with open(norm_out, "w", newline="") as f:

    writer = csv.writer(f, delimiter="\t")

    writer.writerow([
        "Substitution",
        "3XBE_count",
        "3XBE_percent",
        "tEPORCas9_count",
        "tEPORCas9_percent"
    ])

    for mut in NORM_CLASSES:

        n3 = norm_counts["3XBE"][mut]
        nt = norm_counts["tEPORCas9"][mut]

        total3 = sum(norm_counts["3XBE"].values())
        totalt = sum(norm_counts["tEPORCas9"].values())

        writer.writerow([
            mut,
            n3,
            100 * n3 / total3 if total3 else 0,
            nt,
            100 * nt / totalt if totalt else 0
        ])


# --------------------------------------------------
# C>T / G>A summary
# --------------------------------------------------

summary_out = OUT / "CBE_CT_GA_summary.tsv"

with open(summary_out, "w", newline="") as f:

    writer = csv.writer(f, delimiter="\t")

    writer.writerow([
        "Sample",
        "Total_SNV",
        "CT_GA",
        "Other_SNV",
        "CT_GA_percent"
    ])

    for sample in SAMPLES:

        rows = [
            r for r in all_snv_rows
            if r["Sample"] == sample
        ]

        cbe = sum(
            r["CBE_CT_GA"] == "YES"
            for r in rows
        )

        other = len(rows) - cbe

        writer.writerow([
            sample,
            len(rows),
            cbe,
            other,
            100 * cbe / len(rows)
        ])


# --------------------------------------------------
# Print results
# --------------------------------------------------

print("\n==============================")
print("RAW 12-class substitution")
print("==============================")

for sample in SAMPLES:

    total = sum(raw_counts[sample].values())

    print(f"\n{sample} (n={total})")

    for mut in RAW_CLASSES:

        n = raw_counts[sample][mut]

        if n > 0:
            print(
                f"{mut:4s}  "
                f"{n:2d}  "
                f"({100*n/total:5.1f}%)"
            )


print("\n==============================")
print("Normalized 6-class spectrum")
print("==============================")

for sample in SAMPLES:

    total = sum(norm_counts[sample].values())

    print(f"\n{sample} (n={total})")

    for mut in NORM_CLASSES:

        n = norm_counts[sample][mut]

        print(
            f"{mut:4s}  "
            f"{n:2d}  "
            f"({100*n/total:5.1f}%)"
        )


print("\n==============================")
print("CBE-type C>T / G>A")
print("==============================")

for sample in SAMPLES:

    rows = [
        r for r in all_snv_rows
        if r["Sample"] == sample
    ]

    cbe = sum(
        r["CBE_CT_GA"] == "YES"
        for r in rows
    )

    print(
        f"{sample}: "
        f"{cbe}/{len(rows)} "
        f"({100*cbe/len(rows):.1f}%)"
    )


print("\nOutput files:")
print(outfile)
print(raw_out)
print(norm_out)
print(summary_out)
