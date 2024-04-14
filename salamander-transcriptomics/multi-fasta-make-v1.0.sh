#!/bin/bash

# Multi FASTA Make
# Makes a multi FASTA file from the BLAST outputs of Pax6.
# This must be run in the species' blast directory.
# By Holsen B. Moore
# Last Updated: 2024-04-14
# Version: 1.0

# Records the species name to memory.
species=$(pwd | awk -F'/' '{print $5}')

# Records transcript IDs from BLAST search to memory.
script_ids=$(awk -F'\t' '{print $2}' *Pax6*)

# Iterates over each transcript ID.
for i in ${script_ids}
do
	# Creates invisible FASTA file of transcript from BLAST database.
	blastdbcmd -db ${species} -entry ${i} -out .${i}.fasta

	# Adds newline to end of invisible FASTA file.
	# This helps with readability.
	printf "\n" >> .${i}.fasta
done

# Concatenates all invisible FASTA files into a multi FASTA file
cat .*.fasta > ${species}_pax6_multi.fasta

# Removes all invisible FASTA files.
rm .*.fasta
