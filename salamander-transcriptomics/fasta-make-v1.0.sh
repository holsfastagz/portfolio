#!/bin/bash

# FASTA Make
# Makes a FASTA file from the top Pax6 BLAST output.
# This can be run anywhere, and it will create a FASTA file for all the necessary species.
# Last Updated: 2024-04-14
# Version: 1.0 


species_list=( "aquatica" "pterophila_Comal" "pterophila_PC" "wakei" "pricei" "pronatura" "sosorum" "wallacei" "chisholmensis" "cirrigera" "hillisi" "sosorum" "tonkawae_TT" "rathbuni" "neotenes" "naufragia" "NewBraunfels" "latitans_HCC" "latitans_HCSNA" "pronatura_trans" "winkleri" "nana" "waterlooensis" "troglodytes" "tonkawae_BC" )

for i in ${species_list[@]}
do
	# Checks to see if blast directory exists.
	# If not, it will display a message and move on.
	if [ ! -f /stor/work/Hillis/${i}/blast/*Pax6* ]
	then
		printf "No BLAST output found for $i.\n"
		continue
	fi

	# Checks to see if a FASTA Make output exists.
	# If not, it will display a message and move on.
	if [ -f /stor/work/Hillis/${i}/blast/${i}_pax6.fasta ] 
	then
		printf "FASTA Make already completed for $i.\n"
		continue
	fi

	# Records top transcript ID from BLAST search to memory.
	script_id=$(head -n 1 /stor/work/Hillis/${i}/blast/*Pax6* | awk -F'\t' '{print $2}')

	# Writes BLAST output sequence to fasta file.
	blastdbcmd -db /stor/work/Hillis/${i}/blast/${i} -entry ${script_id} -out /stor/work/Hillis/${i}/blast/${i}_pax6.fasta
	
	# Copies fasta file to pax6_scripts directory.
	cp /stor/work/Hillis/${i}/blast/${i}_pax6.fasta /stor/work/Hillis/pax6_scripts

	# Displays a nice message.
	printf "FASTA Make completed for $i"
done
