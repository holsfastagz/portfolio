#!/usr/bin/python3

# Connecting Scripts (BLAST output)
# AWK command generator.
# Holsen B. Moore
# 09 March 2024
# v1

import regex as re

# Read trinity_ids.txt and write to memory.
with open("trinity_ids.txt","r") as file:
    trinity_ids = file.readlines()

# Output AWK commands to terminal to search for read IDs.
for i in trinity_ids:
    lone_id = re.sub(r'\n',"",i)
    print("awk -F $'\\t' '/" + lone_id + "/ {print $1, $5}' ../script_map/*.sam")
