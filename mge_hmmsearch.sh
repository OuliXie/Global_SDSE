#!/bin/bash

# Script for running hmmsearch against pangenome genes
# Looks for recombinase subfamilies, phage structural proteins and essential type IV secretion system proteins
# Uses eggnog-mapper for annotating phage structural proteins - consider running this separately as takes >30 mins - commented out by default
# Can place *.emapper.annotations to directory this script is run from
# Expect the gene_presence_absence_roary.csv file from Panaroo if using Panaroo output
# Expect a folder containing relevant HMMs with subfolders recombinase/ and T4SS/ containing their respective HMMs
# Pan_genome_reference.fa should be translated to protein sequence - see get_pangenome_protein.py 
# Requires csvtk and hmmer3

usage() { echo "usage $(basename $0) [-g path/to/pan_genome_reference.fa] [-i path/to/gene_presence_absence.csv] [-p path/to/hmm/folder]" 1>&2; exit 1; }


while getopts ':g:i:p:h' OPTION; do
    case "${OPTION}" in
        g)
            g=${OPTARG}
            ;;
        i)
            i=${OPTARG}
            ;;
        p)
            p=${OPTARG}
            ;;
        h)
            usage
            ;;
        *)
            usage
            ;;
    esac
done
shift "$(($OPTIND -1))"

if [ -z "${g}" ] || [ -z "${i}" ] || [ -z "${p}" ]; then
    usage
fi

# Make temp directory
temp_dir="$(mktemp -p . -d temp.XXXXXX)"

# List recombinase hmms
ls ${p}/recombinase | grep .hmm > ${temp_dir}/recombinase_hmms.tmp

# Run hmmersearch for each recombinase profile HMM against representative pangenome
# Save hmmersearch results in separate directory for inspection later
mkdir hmm_results
while read hmms
do
    hmmsearch --cut_ga --tblout hmm_results/${hmms%.hmm}.out ${p}/recombinase/${hmms} ${g} > /dev/null
    # Check that there has been a hit
    if grep -qvE "^#" hmm_results/${hmms%.hmm}.out
    then
    # Organise output from hmmsearch taking gene name, recombinase hit and bit score
    grep -vE "^#" hmm_results/${hmms%.hmm}.out | sed 's/ \+/ /g' | cut -d" " -f1,3,6 \
        | sed 's/ /,/g' | csvtk add-header -n Gene,recombinase,score > ${temp_dir}/${hmms%.hmm}.recombinase.clean
    fi
done < ${temp_dir}/recombinase_hmms.tmp

# Gives merged list of gene/recombinase combinations but may have duplicates if multiple hmmer hits
csvtk join -O -f "Gene,recombinase,score" ${temp_dir}/*.recombinase.clean > ${temp_dir}/recombinase.merged

# Take only the highest bit scoring recombinase for each gene
seq_uniq=$(csvtk del-header ${temp_dir}/recombinase.merged | csvtk cut -f 1 | sort -u)
for seq in ${seq_uniq}
# Use non-csvtk commands so don't have to deal with header
do grep -E "^${seq}," ${temp_dir}/recombinase.merged | sort -t"," -k3 -r -g | head -1 >> ${temp_dir}/recombinase.uniq
done
# Add back header and remove bit score column
csvtk add-header -n Gene,recombinase,score ${temp_dir}/recombinase.uniq | \
    csvtk cut -f Gene,recombinase > ${temp_dir}/recombinase.uniq.clean

# Merge recombinase genes with gene_presence_absence.csv
csvtk join -L -f "Gene" ${i} ${temp_dir}/recombinase.uniq.clean > ${temp_dir}/gene_presence_absence_recombinase.csv


### Find T4SS proteins ###
# List T4SS hmms
ls ${p}/T4SS | grep .hmm > ${temp_dir}/T4SS_hmms.tmp

# Run hmmersearch for each T4SS profile HMM (uses E-value threshold of 0.001 as no gathering threshold defined)
# Save hmmsersearch results for inspection later
while read hmms
do
    hmmsearch -E 0.001 --tblout hmm_results/${hmms%.hmm}.out ${p}/T4SS/${hmms} ${g} > /dev/null
    if grep -qvE "^#" hmm_results/${hmms%.hmm}.out
    then
    # Organise output from hmmsearch taking gene name, recombinase hit and bit score
    grep -vE "^#" hmm_results/${hmms%.hmm}.out | sed 's/ \+/ /g' | cut -d" " -f1,3,6 \
        | sed 's/ /,/g' | csvtk add-header -n Gene,T4SS,score > ${temp_dir}/${hmms%.hmm}.T4SS.clean
    fi
done < ${temp_dir}/T4SS_hmms.tmp

# Gives merged list of gene/recombinase combinations but may have duplicates if multiple hmmer hits
csvtk join -O -f "Gene,T4SS,score" ${temp_dir}/*.T4SS.clean > ${temp_dir}/T4SS.merged

# Take only the highest bit scoring T4SS protein for each gene
seq_uniq=$(csvtk del-header ${temp_dir}/T4SS.merged | csvtk cut -f 1 | sort -u)
for seq in ${seq_uniq}
# Use non-csvtk commands so don't have to deal with header
do grep -E "^${seq}," ${temp_dir}/T4SS.merged | sort -t"," -k3 -r -g | head -1 >> ${temp_dir}/T4SS.uniq
done
# Add back header and remove bit score column
csvtk add-header -n Gene,T4SS,score ${temp_dir}/T4SS.uniq | \
    csvtk cut -f Gene,T4SS > ${temp_dir}/T4SS.uniq.clean

# Merge T4SS hmmersearch results with gene_presence_absence file
csvtk join -L -f "Gene" ${temp_dir}/gene_presence_absence_recombinase.csv \
    ${temp_dir}/T4SS.uniq.clean > ${temp_dir}/gene_presence_absence_T4SS.csv

# Find phage structural proteins from eggNOG annotations
# eggnog-mapper may take >30 mins for searching ~6000 genes in SDSE and S. pyo genome so is commented out
# emapper.py -i ${g} -o sdse
# consider running eggnog-mapper separately and move *.emapper.annotations to directory this script is run from
# Uses a word search from COG descriptions and pfam annotations
cut -f1,8,21 *.emapper.annotations | \
    awk 'tolower($0) ~ /phage|capsid protein|portal protein|protein hol|sipho_tail/ \
    && tolower($0) !~ /domain|macrophage|integrase|primase|replisome|abortive|antirepressor|collagen|translation|replication/' | sed 's/,//g' | sed 's/\t/,/' | sed 's/,.*/,1/' \
    | csvtk add-header -n Gene,phage > ${temp_dir}/phage_annotations.csv

# Left merge with gene_presence_absence file
csvtk join -L -f "Gene" ${temp_dir}/gene_presence_absence_T4SS.csv ${temp_dir}/phage_annotations.csv > annotated_gene_presence_absence.csv

# Clean temporary directory
rm -r ${temp_dir}

printf "Script has completed! \nUse annotated_gene_presence_absence.csv with filter_core_pair_summary.sh\n"
