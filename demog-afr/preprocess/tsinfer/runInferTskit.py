#!/usr/bin/env python

import os
import sys
from InferTskit import InferTskit


def main():

    if len(sys.argv) != 2:
        print("Usage: python runInferTskit.py <chromosome>")
        sys.exit(1)

    chrom = sys.argv[1]

    anc_folder = "homo_sapiens_ancestor_GRCh38"
    vcf_folder = "subsetAfricansphased_processed"
    output_dir = "tsinfer_results"

    vcf_file = os.path.join(
        vcf_folder,
        f"subset_chr{chrom}.phased.vcf.gz"
    )

    fasta_file = os.path.join(
        anc_folder,
        f"homo_sapiens_ancestor_{chrom}.fa"
    )

    # checks for file existence
    if not os.path.isfile(vcf_file):
        raise FileNotFoundError(
            f"Missing VCF: {vcf_file}"
        )

    if not os.path.isfile(fasta_file):
        raise FileNotFoundError(
            f"Missing FASTA: {fasta_file}"
        )

    print(f"Running chromosome {chrom}")
    print(f"VCF: {vcf_file}")
    print(f"FASTA: {fasta_file}")


    runner = InferTskit(
        vcf_file=vcf_file,
        output_prefix=f"subset_chr{chrom}",
        output_dir=output_dir,

        ancestral_mode="fasta",
        fasta_file=fasta_file,

        apply_tsdate=True,
        tsdate_time_unit="generations",
        tsdate_mutation_rate=1.2e-8,

        remove_zarr=True
    )

    runner.run()


if __name__ == "__main__":
    main()