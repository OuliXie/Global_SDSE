#!/usr/bin/bash

# Convenience script for downloading hmms from proMGE, pfam and MacSyFinder

mkdir hmms
mkdir hmms/recombinase
mkdir hmms/T4SS

# Get proMGE hmms
wget https://promge.embl.de/data/hmms.tgz
tar -xf hmms.tgz
# Move to correct directory and clean
mv hmm/*.hmm hmms/recombinase
rm -r hmm hmms.tgz

# Get pfam hmms
pfambase="https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/"
pfamhmm="?annotation=hmm"

while IFS=$'\t' read -r pfam name
do
	wget -O ${name}.hmm.gz ${pfambase}${pfam}${pfamhmm}
	gzip -d ${name}.hmm.gz
	mv ${name}.hmm hmms/recombinase
done < pfam_hmms.txt

# Get hmms from MacSyFinder
while IFS=$'\t' read -r path name
do
	wget -O ${name} ${path}
	mv ${name} hmms/T4SS
done < txss_hmms.txt
