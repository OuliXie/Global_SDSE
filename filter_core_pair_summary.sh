#!/bin/bash

# Script for filtering core_pair_summary by max accessory content
# Then calls per_pair_accessory.sh to generate csv with accessory genes between each core pair with accessory >= cutoff

# Requires csvtk, pandas and numpy
# Has been tested with conda install of Panaroo v1.2.10, python v3.7.12, numpy v1.21.6, pandas v1.3.4 and csvtk v0.23.0 

# Gene names cannot contain hyphen "-". Hyphen can only exist between core pair
# order_annotated_gene_presence_absence.sh should be run first in the same directory
# order_annotated_gene_presence_absence.sh will also generate gene_presence_absence files and sequence_id.txt 
# sequence_id.txt file lists all genomes included in the analysis and is required for this script

usage() { echo "usage $(basename $0) 
    [-n max_acc_cutoff] 
    [-i path/to/corekaburra/output/folder] 
    [-r path/to/recombinase_rules.tsv]
    [-g path/to/postpanaroo_gffs/folder]
    [-c core threshold (default 0.99)]
    [-m min accessory content (default 1)]
    [-b ignore gene classifications from sequence breaks]
    [-a min accessory genes in coreless fragments (default 10)]" 1>&2; exit 1; }

while getopts ':n:i:r:g:c:m:ba:h' OPTION; do
    case "${OPTION}" in
        n)
            n=${OPTARG}
            ;;
        i)
            i=${OPTARG}
            ;;
        r)
            r=${OPTARG}
            ;;
        g)
            g=${OPTARG}
            ;;
        c)
            core=${OPTARG}
            ;;
        m)
            min=${OPTARG}
            ;;
        b)
            seq_ig=1
            ;;
        a)
            acc=${OPTARG}
            ;;

        h)
            usage
            ;;
        *)
            usage
            ;;
    esac
done

if [ -z "${min}" ] ; then
    min=1
fi

if [ -z "${core}" ] ; then
    core="0.99"
fi

if [ -z "${acc}" ] ; then
    acc=10
fi

if [ -z "${n}" ] || [ -z "${i}" ] || [ -z "${r}" ] || [ -z "${g}" ] ; then
    usage
fi

shift $(($OPTIND -1))

core_pair="${i}/core_pair_summary.csv"
low_freq="${i}/low_frequency_gene_placement.tsv"

# Find all core pairs with max_acc > cutoff
csvtk filter -f "max_acc>=${n}" ${core_pair} | csvtk del-header > filtered_core_pair.csv

# Add progress bar
prog=0
total=$(wc -l filtered_core_pair.csv | grep -o '[0-9]\+')
bar=##################################################

while IFS=$',' read pair n occ_1 occ_2 co_occ minim max mean median min_acc max_acc mean_acc median_acc
do
    ((prog++))
    perc=$( expr ${prog} \* 100 / ${total} )
    bar_prog=$( expr ${perc} / 2 )
    gene_1=$(echo ${pair} | cut -d"-" -f1)
    gene_2=$(echo ${pair} | cut -d"-" -f2)
    output="per_pair_output/${gene_1}-${gene_2}"
    python $(dirname $0)/per_pair_accessory.py -f ${gene_1} -s ${gene_2} -l ${low_freq} -r ${r} -c ${core} \
      -m ${min} -o per_pair_output -i sequence_id.txt -g ordered_gene_presence_absence
    # Check if segment sequence break gene classifications is on
    if [ $seq_ig == 1 ] ; then
        if [ ${gene_1} == Sequence_break ] || [ ${gene_2} == Sequence_break ] ; then
            python $(dirname ${0})/mge_type.py -i ${output}/*.filtered.temp -r ${r} -o ${output} -b
        else
            python $(dirname ${0})/mge_type.py -i ${output}/*.filtered.temp -r ${r} -o ${output}
        fi
    else
        python $(dirname ${0})/mge_type.py -i ${output}/*.filtered.temp -r ${r} -o ${output}
    fi
    # clean temporary files
    rm ${output}/*.temp
    echo -ne "\r${perc}% ${bar:0:${bar_prog}}"
done < filtered_core_pair.csv
echo -e "\nAnalysing accessory only segments..."

# Analyse coreless segments
python $(dirname $0)/coreless_extract.py -a ${i}/coreless_contig_accessory_gene_content.tsv \
    -f ${g} -o per_pair_output -i sequence_id.txt -g ordered_gene_presence_absence \
    -c ${core} -m ${acc}
cp per_pair_output/coreless_contigs/coreless_segments.csv .
# Find subdirectories of per_pair_output/coreless_contigs
ls -l per_pair_output/coreless_contigs | grep -e "^d" | \
    awk '{print $NF}' > per_pair_output/coreless_contigs/coreless_list.txt
# Classify segments
while read coreless_dir
do
    if [ $seq_ig == 1 ] ; then
        python $(dirname $0)/mge_type.py -i per_pair_output/coreless_contigs/${coreless_dir}/*.filtered.temp -r ${r} \
        -o per_pair_output/coreless_contigs/${coreless_dir} -b
    else
        python $(dirname $0)/mge_type.py -i per_pair_output/coreless_contigs/${coreless_dir}/*.filtered.temp -r ${r} \
        -o per_pair_output/coreless_contigs/${coreless_dir}
    fi
    # clean temporary files
    rm per_pair_output/coreless_contigs/${coreless_dir}/*.temp
done < per_pair_output/coreless_contigs/coreless_list.txt

echo -e "\nGenerating summary outputs..."

# Create summary tsv with accessory segment gene counts for each sequence and core-core pair
csvtk -t join -L -f Gff per_pair_output/*/*_joined_count.tsv > summary_acc_count.tsv

# Create summary tsv with accessory segment lengths for each sequence and core-core pair
csvtk -t join -L -f Gff per_pair_output/*/*_joined_distance.tsv > summary_dist_acc_summary.tsv

# Create summary csv with MGE count from resolved recombinases (includes counting nested MGEs)
# List all *_mge_full.csv files
find per_pair_output -name "*_mge_full.csv" -type f > list_mge_full.txt
if [ $seq_ig == 1 ] ; then
    python $(dirname $0)/summarise_mge_count.py -i list_mge_full.txt -s sequence_id.txt -b
else
    python $(dirname $0)/summarise_mge_count.py -i list_mge_full.txt -s sequence_id.txt
fi
rm list_mge_full.txt

# Create summary csv with summarised resolved recombinases (collapses nested MGEs without counts of components)
sed 's/&[0-9]//g' summary_mge_count.csv > summary_mge_classes.csv

# Create gene_presence_absence style csv with classification of each accessory to MGE or non-MGE element type
find per_pair_output -name "*_mge_summarised.csv" -type f > list_mge_summarised.txt
python $(dirname $0)/summarise_gene_type.py -i list_mge_summarised.txt \
    -g ordered_gene_presence_absence/annotated_gene_presence_absence.csv -c ${core}
rm list_mge_summarised.txt

echo "Done!"
