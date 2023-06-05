import os
import argparse
import pandas as pd
import numpy as np
from tqdm import tqdm


# Script for extracting coreless segment to look for missed phage/ICE
# By default looks for segments with >=10 accessory genes as more interested in missed phage/ICE here rather than IS
# Run order_annotated_gene_presence_absence.sh first in same directory
# order_annotated_gene_presence_absence.sh will also generate a sequence_id.txt file listing all genomes


def parse_gff(sequence, gff_dir, coreless_contigs):
    # Adapted from post_run_gff_output.py
    # Read in and split off FASTA portion
    raw_file = open(os.path.join(gff_dir, sequence + ".gff"), 'r')
    # Read in and split off FASTA portion
    lines = raw_file.read().replace(',', '')
    split = lines.split('##FASTA')[0]
    # Remove headers
    body = []
    contig_length = {}
    for line in split.splitlines():
        if "##sequence-region" in line:
            line_split = line.split(" ")
            # Capture contig and contig length
            # Gives a dictionary with contig number as the key and contig length as value
            if line_split[1] in coreless_contigs:
                contig_length[line_split[1]] = [line_split[3]]
        elif "##" not in line:
            body.append(line)
    # Find locus tags in contig
    # Gives dictionary with contig number as the key and a list of locus tags as the value
    parsed_locus = {}
    for gff_line in body:
        initial_split = gff_line.split("\t")
        if initial_split[0] in coreless_contigs:
            attribute_split = " ".join(initial_split[8:]).split(";")
            for attribute in attribute_split:
                if "locus_tag=" in attribute:
                    locus_tag = attribute.split("=")[1]
                    if initial_split[0] in parsed_locus.keys():
                        parsed_locus[initial_split[0]].append(locus_tag)
                    else:
                        parsed_locus[initial_split[0]] = [locus_tag]

    return parsed_locus, contig_length


def extract_coreless(shortlist, coreless, gff_dir, gpa_dir, core, output_dir):
    # Extract accessory genes from coreless segment using annotated gene_presence_absence files
    # Extract locus tags from post-Panaroo corrected gffs

    # Create empty dataframe for contig length
    length_all = pd.DataFrame(columns=["Gff", "Contig", "distance", "dummy_contig"])
    coreless["dummy_contig"] = ""

    for sequence in tqdm(shortlist.keys()):
        segment = pd.read_csv(gpa_dir + "/" + sequence + "_gene_presence_absence.csv").fillna("")
        segment = segment.rename(columns={"Non-unique Gene name": "Non-unique_name", "Annotation": "Annotation",
                                          "No. isolates": "No_isolates", "No. sequences": "No_sequences"})
        # Temp column containing only the first locus_tag for CDS which were found to belong to the same gene
        unique_tag = [x.split(";")[0] for x in segment[sequence].to_list()]
        segment = segment.assign(locus_unique=unique_tag)

        # Find locus tags
        accessory_locus, contig_length = parse_gff(sequence, gff_dir, shortlist[sequence])

        # Make dataframe with columns "Gff", "Contig" and "segment_length"
        length_df = pd.DataFrame(data=contig_length)
        length_df = length_df.assign(Gff=sequence)
        length_df = pd.melt(length_df, id_vars="Gff", value_vars=length_df.columns,
                            var_name="Contig", value_name="distance")
        length_df["dummy_contig"] = ""

        i = 0
        # Find index in coreless dataframe where sequence appears
        first_index = np.where(coreless["Gff"] == sequence)[0][0]
        # Get accessory segments
        for contig, locus_tags in accessory_locus.items():
            contig_segment = segment[segment.loc[:, "locus_unique"].isin(locus_tags)]
            # Remove temporary "locus_unique" column
            contig_segment = contig_segment.drop("locus_unique", axis=1)
            # Remove readthrough into core genes
            contig_segment = contig_segment[contig_segment["No_isolates"] < core]
            # Write into subdirectory with generic contig name - this is important for summarise_mge_count.py
            subdir = os.path.join(output_dir, "coreless_" + str(i + 1))
            if not os.path.exists(subdir):
                os.mkdir(subdir)
            contig_segment.to_csv(os.path.join(subdir, sequence + ".filtered.temp"), index=False)
            # Add dummy variable to length_df and coreless df
            coreless.iloc[first_index + i, -1] = "coreless_" + str(i + 1)
            length_df.iloc[i, -1] = "coreless_" + str(i + 1)
            i += 1

        # Merge with overall record
        length_all = pd.concat([length_all, length_df])

    # Merge contig_length with coreless
    coreless = pd.merge(coreless, length_all, how="left", on=["Gff", "Contig", "dummy_contig"]).fillna(0)

    return coreless


def main():
    description = "Extracts coreless segments with minimum number of genes (default 10)"
    parser = argparse.ArgumentParser(description=description)

    io_opts = parser.add_argument_group("Input/output")

    io_opts.add_argument(
        "-a",
        "--acc",
        dest="acc_only",
        help="coreless_contig_accessory_gene_content.tsv",
        required=True,
        type=os.path.abspath,
    )

    io_opts.add_argument(
        "-f",
        "--gff",
        dest="gff",
        help="Path to directory containing post-Panaroo corrected gff",
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
        help="Core threshold expressed as fraction up to 1 (default 0.99)",
        required=False,
        type=float,
        default=0.99
    )

    io_opts.add_argument(
        "-m",
        "--min",
        dest="min",
        help="Minimum number of genes in accessory fragment cutoff (default 10)",
        required=False,
        type=int,
        default=10
    )

    args = parser.parse_args()

    # Make dedicated directory for outputs
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)

    # Make directory for core-core combination
    output = os.path.join(args.output_dir, "coreless_contigs")
    if not os.path.exists(output):
        os.mkdir(output)

    # Read in coreless_contig_accessory_gene_content.tsv
    coreless = pd.read_csv(args.acc_only, sep="\t")
    # Filter for accessory genes >=min (default 10) to look for missed phage/ICE elements
    # Chosen 10 for balance between finding missed elements and including potentially poorly assembled contigs or
    # very small fragmented MGEs
    # If interested in small plasmids, consider reducing cutoff
    coreless = coreless[coreless["Accessory_count"] >= args.min]

    # Create dictionary with key as genome name and value as list of contigs of interest
    shortlist = coreless[["Gff", "Contig"]].groupby("Gff")["Contig"].apply(list).to_dict()

    # Read in sequence list
    seq_id = pd.read_csv(args.seq_id, header=None, names=["Gff"])

    # Calculate number of core genes and core gene cutoff
    core = np.floor(len(seq_id.index) * args.core)

    # Extract accessory only contigs with >= 10 accessory genes to look for missed phage/ICE
    coreless = extract_coreless(shortlist, coreless, args.gff, args.gpa, core, output)
    coreless["distance"] = coreless["distance"].astype(int)
    coreless["Accessory_count"] = coreless["Accessory_count"].astype(int)

    # Save this output in case want to trace back manually
    coreless.to_csv(os.path.join(output, "coreless_segments.csv"), index=False)

    # Pivot after joining contig length (as surrogate for acc distance)
    # Pivot coreless dataframe wider for merging with sequence list
    coreless_count = coreless.pivot(
        index="Gff", columns="dummy_contig", values="Accessory_count").fillna(0).\
        astype(int).replace(0, "").astype(object)
    coreless_distance = coreless.pivot(
        index="Gff", columns="dummy_contig", values="distance").fillna(0).astype(int).replace(0, "").astype(object)

    # Create an ordered sequence list
    # This will be merged later to create a master list
    # Left merge low_freq with seq_id to maintain order of Gffs
    acc_joined = seq_id.merge(coreless_count, how="left", on="Gff")
    dist_joined = seq_id.merge(coreless_distance, how="left", on="Gff")

    # Output accessory gene count
    acc_joined.to_csv(os.path.join(output, "coreless_joined_count.tsv"), sep="\t", index=False)
    dist_joined.to_csv(os.path.join(output, "coreless_joined_distance.tsv"), sep="\t", index=False)


if __name__ == "__main__":
    main()
