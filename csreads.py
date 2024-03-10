#!/usr/bin/python3

# Connecting Scripts (Reads)
# Read ID duplicate removal and AWK command generator.
# Holsen B. Moore
# 09 March 2024
# v1

import regex as re

# Open the contents of raw_output.txt and write into memory.
with open("raw_output.txt", "r") as file:
    redundant_ids = file.readlines()

# Remove redundant items from the list and output to terminal.
purged_ids = list(set(redundant_ids))
for i in purged_ids:
    lone_id = re.sub(r'\n',"",i)
    print(lone_id)

print("")

# Format the IDs into a new AWK command to search for transcript IDs.
script_awk = []
for i in purged_ids:
    lone_id = re.sub(r'\n',"",i)
    lone_id = re.sub(r' *',"",lone_id)
    script_awk.append("awk -F $'\\t' '/" + lone_id + "/ {print $3, $5}' ../trim_map/*.sam >> *pax6_scripts.txt")

# Output the generated AWK commands to terminal.
for i in script_awk:
    print(i)
