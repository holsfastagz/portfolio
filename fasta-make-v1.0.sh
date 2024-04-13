#!/bin/bash

# FASTA Make
# Makes a FASTA file from the BLAST outputs of Pax6.
# Last Updated: 2024-04-12
# Version: 1.1

# Records the species name to memory.
species=$(pwd | awk -F'/' '{print $5}')

# Records the species abbreviation.
# This is the same as the species name for many species.
abbr=$(ls ../raw/ | grep "L001_R1" | awk -F'_' '{print $1}')

# Records the sample number, e.g. "S1."
sample_no=$(ls ../raw/ | grep "L001_R1" | awk -F'_' '{print $2}')

# Records transcript IDs from BLAST search to memory.
script_ids=$(awk -F'\t' '{print $2}' *Pax6*)

for i in ${script_ids}
do
	blastdbcmd -db ${species} -entry ${i} -out .${i}.fasta
	printf "\n" >> .${i}.fasta
done

blastdbcmd -db ${species} -entry ${script_ids[1]} -out ${species}_pax6_top.fasta
cat .*.fasta > ${species}_pax6_cat.fasta
