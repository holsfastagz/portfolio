#!/bin/bash

# Connecting Scripts
# Version 1.0
# Salamander Transcriptomics
# Biodiversity Discovery
# March 12, 2024
# Holsen B. Moore

echo 'Program is running. This may take a while.'

# Creates the temporary directory and temporary file for purged Trinity IDs.
mkdir .cs_temp 
touch .cs_temp/trinity-ids.tmp

# Gathers Trinity IDs from BLAST output file and writes them to memory.
redundant_trinity_ids=$(awk -F $'\t' '{print $2}' *Pax6_Opsin.txt)

# Removes duplicates from Trinity IDs and writes to file.
for i in ${redundant_trinity_ids}
do
	grep ${i} .cs_temp/trinity-ids.tmp || 
	echo ${i} >> .cs_temp/trinity-ids.tmp 
done

# Writes content of .cs_temp/trinity-ids.tmp to memory.
purged_trinity_ids=$(cat .cs_temp/trinity-ids.tmp)

# Uses grep to grab lines including Trinity IDs and writes them to a file.
touch .cs_temp/read_ids.sam.tmp
grep_for_reads=""
for i in ${purged_trinity_ids}
do
	grep_for_reads+="${i}|"
done
grep -E "${grep_for_reads%|}" ../script_map/*.sam >> .cs_temp/read_ids.sam.tmp

# Writes read IDs to file corresponding to their Trinity IDs.
for i in ${purged_trinity_ids}
do
	touch .cs_temp/${i}
	redundant_read_ids=$(awk -F $'\t' "BEGIN {OFS = \"~\"} /${i}/ {print \$1, \$5}" .cs_temp/read_ids.sam.tmp)
	for j in ${redundant_read_ids}
	do
		grep "${j}" .cs_temp/${i} ||
		echo ${j} >> .cs_temp/${i}
	done
done

# Removes extraneous "@SQ~" sequence from Trinity ID files.
for i in ${purged_trinity_ids}
do
	sed -i "/@SQ~/d" .cs_temp/${i}
done

# Creates temporary output file.
output=".output.tmp"
touch ${output}

# Uses grep to grab lines including Read IDs and writes them to a file.
touch .cs_temp/script_ids.sam.tmp
grep_for_scripts=""
for i in ${purged_trinity_ids}
do
	local_read_ids=$(cat .cs_temp/${i})
	for j in ${local_read_ids}
	do
		grep_for_scripts+="${j%~*}|"
	done
done
grep -E "${grep_for_scripts%|}" ../trim_map/*.sam >> .cs_temp/script_ids.sam.tmp

# Writes Trinity IDs to file. Uses awk to find transcript IDs from file. Writes Read IDs and read quality, and transcript IDs and read quality to output file.
for i in ${purged_trinity_ids}
do
	printf "$i\n" >> ${output}
	local_read_ids=$(cat .cs_temp/${i})
	for j in ${local_read_ids}
	do
		printf "\t$j\n" >> ${output}
		mkdir .cs_temp/${i}-${j}
		awk -F $'\t' "BEGIN {OFS = \"~\"} /${j%~*}/ {print \$3, \$5}" .cs_temp/script_ids.sam.tmp >> .cs_temp/${i}-${j}/${j}_redundant.tmp
		redundant_script_ids=$(cat .cs_temp/${i}-${j}/${j}_redundant.tmp)
		touch .cs_temp/${i}-${j}/${j}_purged.tmp
		for k in ${redundant_script_ids}
		do
			grep "${k}" .cs_temp/${i}-${j}/${j}_purged.tmp ||
			echo ${k} >> .cs_temp/${i}-${j}/${j}_purged.tmp
		done
		purged_script_ids=$(cat .cs_temp/${i}-${j}/${j}_purged.tmp)
		for k in ${purged_script_ids}
		do
			printf "\t\t$k\n" >> ${output}
		done
	done
done	

# Formatting and background stuff.
sed -i "s/~/\t/" ${output}
rm -r .cs_temp

# Finds species name by isolating species directory name.
parent_dir=$(echo ${PWD%/*})
species=$(echo ${parent_dir#*Hillis/})

# Moves temporary output file to permanent output file including species name.
mv ${output} ${species}_pax6_scripts.txt

echo 'Program has finished!'
