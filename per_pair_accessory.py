import os
import argparse
import pandas as pd
import numpy as np


# Script for extracting user specified accessory genes between core-core segments
# Run order_annotated_gene_presence_absence.sh first in same directory
# order_annotated_gene_presence_absence.sh will also generate a sequence_id.txt file listing all genomes

# If there is readthrough to core genes, these will be automatically removed


def extract(shortlist, start, finish, gpa, core, output):
    # Extract accessory genes between core segments from ordered and annotated gene_presence_absence files
    # Check first gene is not Sequence_break which is non-unique
    if start != "Sequence_break":
        if finish == "Sequence_break":
            for sequence in shortlist.keys():
                # Read in segment
                segment = pd.read_csv(gpa + "/" + sequence + "_gene_presence_absence.csv").fillna("")
                segment = segment.rename(columns={"Non-unique Gene name": "Non-unique_name", "Annotation": "Annotation",
                                                  "No. isolates": "No_isolates", "No. sequences": "No_sequences"})
                index = segment.index[segment["Gene"] == start][0]
                # Need to add one as pandas slice excludes end index
                n_acc = shortlist[sequence]["Core_region_accessory_count"]
                index_fow = index + n_acc + 1
                index_rev = index - n_acc
                # Check we will slice in right direction - sometimes flipped
                for_acc = len(segment.iloc[(index + 1):index_fow, 6][segment.iloc[(index+1):index_fow, 6] < core])
                rev_acc = len(segment.iloc[index_rev:index, 6][segment.iloc[index_rev:index, 6] < core])
                if for_acc >= rev_acc:
                    segment = segment.iloc[(index + 1):index_fow, :]
                else:
                    segment = segment.iloc[index_rev:index, :]
                # Remove readthrough into core genes
                segment = segment[segment["No_isolates"] < core]
                # Write segments to csv
                gff = segment.columns[8]
                segment.to_csv(output + "/" + gff + ".filtered.temp", index=False)
        # If both genes listed
        else:
            for sequence in shortlist.keys():
                segment = pd.read_csv(gpa + "/" + sequence + "_gene_presence_absence.csv").fillna("")
                segment = segment.rename(columns={"Non-unique Gene name": "Non-unique_name", "Annotation": "Annotation",
                                                  "No. isolates": "No_isolates", "No. sequences": "No_sequences"})
                index_start = segment.index[segment["Gene"] == start][0]
                index_finish = segment.index[segment["Gene"] == finish][0]
                # Check which direction to slice in
                if index_start < index_finish:
                    segment = segment.iloc[(index_start + 1):index_finish, :]
                else:
                    segment = segment.iloc[(index_finish + 1):index_start, :]
                # Remove readthrough into core genes
                segment = segment[segment["No_isolates"] < core]
                # Write segments to csv
                gff = segment.columns[8]
                segment.to_csv(output + "/" + gff + ".filtered.temp", index=False)

    # In the case of first gene being "Sequence_break", index from second gene
    else:
        for sequence in shortlist.keys():
            # Read in gene_presence_absence file
            segment = pd.read_csv(gpa + "/" + sequence + "_gene_presence_absence.csv").fillna("")
            segment = segment.rename(columns={"Non-unique Gene name": "Non-unique_name", "Annotation": "Annotation",
                                              "No. isolates": "No_isolates", "No. sequences": "No_sequences"})
            index = segment.index[segment["Gene"] == finish][0]
            # Need to add one as pandas slice excludes end index
            n_acc = shortlist[sequence]["Core_region_accessory_count"]
            index_fow = index + n_acc + 1
            index_rev = index - n_acc
            # Check we will slice in right direction - sometimes flipped
            for_acc = len(segment.iloc[(index + 1):index_fow, 6][segment.iloc[(index + 1):index_fow, 6] < core])
            rev_acc = len(segment.iloc[index_rev:index, 6][segment.iloc[index_rev:index, 6] < core])
            if rev_acc >= for_acc:
                segment = segment.iloc[index_rev:index, :]
            else:
                segment = segment.iloc[(index + 1):index_fow, :]
            # Remove readthrough into core genes
            segment = segment[segment["No_isolates"] < core]
            # Write segments to csv
            gff = segment.columns[8]
            segment.to_csv(output + "/" + gff + ".filtered.temp", index=False)


def main():
    description = "Pulls accessory segment between specified core-core pair"
    parser = argparse.ArgumentParser(description=description)

    io_opts = parser.add_argument_group("Input/output")

    io_opts.add_argument(
        "-f",
        "--first",
        dest="f_gene",
        help="First gene in core-core pair",
        required=True,
        type=str,
    )

    io_opts.add_argument(
        "-s",
        "--second",
        dest="s_gene",
        help="Second gene in core-core pair",
        required=True,
        type=str,
    )

    io_opts.add_argument(
        "-l",
        "--low_freq",
        dest="low_freq",
        help="low_frequency_gene_placement.tsv",
        required=True,
        type=os.path.abspath,
    )

    io_opts.add_argument(
        "-r",
        "--rec",
        dest="rec",
        help="recombinase_rules.tsv",
        required=True,
        type=os.path.abspath,
    )

    io_opts.add_argument(
        "-o",
        "--output",
        dest="output_dir",
        help="Output directory",
        required=True,
        type=os.path.abspath,
    )

    io_opts.add_argument(
        "-i",
        "--seq_id",
        dest="seq_id",
        help="sequence_id.txt",
        required=True,
        type=os.path.abspath,
    )

    io_opts.add_argument(
        "-g",
        "--gpa",
        dest="gpa",
        help="Path to ordered_gene_presence_absence folder",
        required=True,
        type=os.path.abspath,
    )

    io_opts.add_argument(
        "-c",
        "--core",
        dest="core",
        help="Core threshold expressed as fraction up to 1",
        required=False,
        type=float,
        default=0.99
    )

    io_opts.add_argument(
        "-m",
        "--min",
        dest="min",
        help="Minimum accessory content",
        required=False,
        type=int,
        default=1
    )

    args = parser.parse_args()

    # Make dedicated directory for outputs
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)

    # Make directory for core-core combination
    output = os.path.join(args.output_dir, args.f_gene + "-" + args.s_gene)
    if not os.path.exists(output):
        os.mkdir(output)

    # read in low_frequency_gene_placement.tsv
    low_freq = pd.read_csv(args.low_freq, sep="\t")
    low_freq = low_freq[(low_freq["Core_gene_1"] == args.f_gene) &
                        (low_freq["Core_gene_2"] == args.s_gene) &
                        (low_freq["Core_region_accessory_count"] >= args.min)]

    # Create an ordered sequence list
    # This will be merged later to create a master list
    seq_id = pd.read_csv(args.seq_id, header=None, names=["Gff"])
    # Left merge low_freq with seq_id to maintain order of Gffs
    joined = seq_id.merge(low_freq, how="left", on="Gff")

    # Change to integer
    joined[["Core_region_accessory_count", "Core_region_size"]] = \
        joined[["Core_region_accessory_count", "Core_region_size"]].fillna(0).astype(int)

    # Output accessory gene count
    acc_joined = joined[["Gff", "Core_region_accessory_count"]].astype(str).replace("0", "")
    acc_joined = acc_joined.rename({"Core_region_accessory_count": (args.f_gene + "-" + args.s_gene)}, axis=1)
    acc_joined.to_csv(output + "/" + args.f_gene + "-" + args.s_gene + "_joined_count.tsv", sep="\t", index=False)
    dist_joined = joined[["Gff", "Core_region_size"]].astype(str).replace("0", "")
    dist_joined = dist_joined.rename({"Core_region_size": (args.f_gene + "-" + args.s_gene)}, axis=1)
    dist_joined.to_csv(output + "/" + args.f_gene + "-" + args.s_gene + "_joined_distance.tsv", sep="\t", index=False)

    # Calculate number of core genes and core gene cutoff
    core = np.floor(len(seq_id.index) * args.core)

    # Sequences to examine
    shortlist = low_freq[["Gff", "Core_region_accessory_count"]].set_index("Gff").to_dict("index")

    # Extract accessory genes between core segments from ordered and annotated gene_presence_absence files
    extract(shortlist, args.f_gene, args.s_gene, args.gpa, core, output)


if __name__ == "__main__":
    main()
