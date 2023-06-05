import os
import argparse
import pandas as pd
from tqdm import tqdm


def summarise_mge(input_file, summary_df, ignore):

    mge_types = ["IS", "Phage", "Phage_like", "ICE", "ME", "Integron", "Degraded", "Hotspot"]

    core_pair = os.path.dirname(input_file).split("/")[-1]
    segment = pd.read_csv(input_file)
    gff = list(segment.columns)[8]
    mge_count = {}

    # check df not empty
    if not segment.empty:
        # record frequencies of each nested mge type
        for mge in mge_types:
            n = segment[mge].unique()[0]
            if n > 0:
                mge_count[mge] = n

        # if no mge detected
        if len(mge_count) == 0:
            # if ignoring classifications with sequence break and no MGE detected
            if ignore and "Sequence_break" in core_pair:
                summary_df.loc[gff, core_pair] = "Unclassified"
            else:
                summary_df.loc[gff, core_pair] = "Non_MGE"
        # fill gff/core-core cell with element type and count e.g., IS&2;Phage&1
        else:
            value = []
            for x in mge_count.keys():
                value.append(x + "&" + str(mge_count[x]))
            summary_df.loc[gff, core_pair] = ';'.join(value)
    else:
        return summary_df

    return summary_df


def main():
    description = "Summarises mge count per sequence and core-core pair"
    parser = argparse.ArgumentParser(description=description)

    io_opts = parser.add_argument_group("Input/output")

    io_opts.add_argument(
        "-i",
        "--input",
        dest="input_file",
        help="List of classified MGE csvs for each sequence between a single core-core pair",
        required=True,
        type=os.path.abspath,
    )

    io_opts.add_argument(
        "-s",
        "--sequences",
        dest="seq_id",
        help="File listing sequence/gff ids",
        required=True,
        type=os.path.abspath,
    )

    io_opts.add_argument(
        "-o",
        "--output",
        dest="output_dir",
        help="Output directory",
        required=False,
        type=os.path.abspath,
    )

    parser.add_argument(
        "-b",
        action="store_true",
        help="Ignore gene classification if no MGE found and sequence break in core pair",
    )

    args = parser.parse_args()

    # create list of input csvs paths
    csvs = []
    with open(args.input_file) as f:
        for line in f:
            line = line.strip()
            csvs.append(line)

    # create dataframe with sequence id/gff as index
    summary_df = pd.read_csv(args.seq_id, header=None, names=["Gff"])

    # add core-core pair columns
    core_pair = list(dict.fromkeys([os.path.dirname(x).split("/")[-1] for x in csvs]))
    summary_df = pd.concat([summary_df, pd.DataFrame(columns=core_pair)], axis=1).fillna("")

    summary_df.set_index("Gff", inplace=True)

    print("Summarsing MGE counts...")
    for file in tqdm(csvs):
        summarise_mge(file, summary_df, args.b)

    if args.output_dir is None:
        summary_df.to_csv("summary_mge_count.csv")
    else:
        summary_df.to_csv(args.output_dir + "/" + "summary_mge_count.csv")


if __name__ == "__main__":
    main()
