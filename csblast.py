#!/usr/bin/python3

# Connecting Scripts (BLAST output)
# AWK command generator.
# Holsen B. Moore
# 09 March 2024
# v1

import regex as re

# Read trinity_ids.txt and write to memory.
with open("trinity_ids.txt","r") as file:
    redundant_ids = file.readlines()

# Remove duplicate IDs and output to terminal.
purged_ids = list(set(redundant_ids))
for i in purged_ids:
    lone_id = re.sub(r'\n', "",i)
    print(lone_id)

print("")

# Output AWK commands to terminal to search for read IDs.
for i in purged_ids:
    lone_id = re.sub(r'\n',"",i)
    print("awk -F $'\\t' '/" + lone_id + "/ {print $1, $5}' ../script_map/*.sam")
