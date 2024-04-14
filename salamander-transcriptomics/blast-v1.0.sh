#!/bin/bash

# BLAST v1.0
# Let's have a BLAST!
# Last Updated: 2024-04-12

# Records the species name to memory.
species=$(pwd | awk -F'/' '{print $5}')

# Records the species abbreviation.
# This is the same as the species name for many species.
abbr=$(ls ../raw/ | grep "L001_R1" | awk -F'_' '{print $1}')

# Records the sample number, e.g. "S1."
sample_no=$(ls ../raw/ | grep "L001_R1" | awk -F'_' '{print $2}')

# Create a BLAST database of the assemblies in a new directory called blast.
makeblastdb -in ../trinity/${abbr}_transcript.fasta -parse_seqids -dbtype nucl -out ${species}
wait

# Protein sequences can be found from NCBI https://www.ncbi.nlm.nih.gov/nuccore/U36574.1?report=fasta.

# Run new data base against sequence (inside the blast directory) using code below.
tblastn -query ../../homo_sequences/Pax6_Cynops.fasta -db ${species} -out ${species}_Pax6_Opsin.txt -outfmt 6 -evalue 1e-20
wait

# This should not take long. Open output file to see matching transcripts.

# Connecting Transcripts

# This process requires a BLAST output of an assembled transcriptome, a mapping of reads to the same assembled transcriptome, and a mapping of reads to the reference transcriptome

# The reference map should be in a directory called trim_map, the BLAST output should be in a directory called blast, and the assembled transcriptome map should be in a directory called script_map.

# Copy Connecting Scripts Version 1.0 into the blast directory.
cp /stor/work/Hillis/programs/connecting-scripts-v1.0.sh ./

# Run Connecting Scripts Version 1.0.
bash connecting-scripts-v1.0.sh 

# Output will be saved by the program in a file like so: ${species}_pax6_scripts.txt
# You do __NOT__ have to make this file yourself.
