#!/bin/bash

#### Line deletion program for Trinity failed commands - Biodiversity Discovery Spring 2024.
#### By Holsen B. Moore :)
#### Last Updated 2024-02-15

#### Use this program in /trinity/trinity_out_dir if you have failed commands you need to delete.

#### Global variables:

dir_exists=1      # Used to record if ../original_recursive_trinity.cmds existed before the use of this program.
backup_exists=1   # Used to record if the recursive_trinity.cmds backup existed before the use of this program.
command_ids=()    # Used in Step 2 to record command IDs from FailedCommands.
progress=0        # Used in Step 4 to display progress.
periods=""        # Used in Step 4 to display progress.
original_length=0 # Used to record the original length of recursive_trinity.cmds in lines.
final_length=0    # Used to record the final length of recursive_trinity.cmds in lines.

#### Step 0: Check if the correct files exist in the current directory. ####

[ ! -f 'FailedCommands' ] && [ ! -f 'recursive_trinity.cmds' ] &&
echo -e "\nThe necessary files do not exist in this directory. Perhaps you are in the wrong directory or haven't run Trinity yet.\n:)" && exit 1

#### Step 1: Create a backup of the recursive_trinity.cmds file in case we make a mistake, provided a backup doesn't already exist. ####

[ ! -d '../original_recursive_trinity.cmds' ] &&
mkdir ../original_recursive_trinity.cmds && dir_exists=0

[ ! -f '../original_recursive_trinity.cmds/recursive_trinity.cmds' ] &&
cp recursive_trinity.cmds ../original_recursive_trinity.cmds && backup_exists=0

[ $backup_exists -eq 0 ] &&
echo -e '\n[Step 1]: Backup of recursive_trinity.cmds created in ../original_recursive_trinity.cmds.' ||
echo -e '\n[Step 1]: Backup of recursive_trinity.cmds already exists in ../original_recursive_trinity.cmds.'

#### Step 2: Gather the command IDs from the FailedCommands file. ####

command_ids=( $(sed -n 's@\.trinity\.reads\.fa\.out.*@\\.@
                 s@.*/c@c@p' FailedCommands) )
echo -e "\n[Step 2]: Gathered Command IDs from FailedCommands.\nThe following ${#command_ids[@]} commands were gathered:\n
${command_ids[@]%??}\n"

#### Step 3: Delete failed commands from recursive_trinity.cmds following confirmation. ####

original_length=$(wc -l recursive_trinity.cmds)

read -p "Do you want to delete these commands? (y/N): " del
case $del in
    [yY])
    echo -e "\nDeleting commands. Please wait..."
    ;;
    *)
    echo "The program has been cancelled. Nothing has been deleted."
    [ $backup_exists -eq 0 ] && rm ../original_recursive_trinity.cmds/recursive_trinity.cmds &&
    [ $dir_exists -eq 0 ] && rmdir ../original_recursive_trinity.cmds
    exit 0
    ;;
esac

for n in ${command_ids[@]}
do
    (( ++progress ))
    periods+="."
    sed -i "/$n/d" recursive_trinity.cmds
    echo -ne "    "$(( 100 * $progress / ${#command_ids[@]} ))"% Done $periods\r"
done

final_length=$(wc -l recursive_trinity.cmds)

echo -e "\n\n[Step 3]: Deleted ${#command_ids[@]} failed commands from recursive_trinity.cmds.\nIf all is well, you may now restart Trinity :)\n
File Length (Lines)\nBefore: $original_length\n After: $final_length\n"

