#!/bin/bash

#### CMDSDEL: Line deletion program for Trinity failed commands - Biodiversity Discovery Spring 2024. ####
#### By Holsen B. Moore :) ####
#### Last Updated 2024-01-11 ####

#### Use this program in /trinity/trinity_out_dir if you have failed commands you need to delete. ###

#### Variables: ####

dirq=1    # Used to record if ../original_recursive_trinity.cmds existed before the use of this program. #
recq=1    # Used to record if the recursive_trinity.cmds backup existed before the use of this program. #
idar=()   # Used in Step 2 to record command IDs from FailedCommands. #
tic=0     # Used in Step 4 to display progress. #
per=""    # Used in Step 4 to display progress. #
orlen=0   # Used to record the original length of recursive_trinity.cmds in lines. #
delen=0   # Used to record the final length of recursive_trinity.cmds in lines. #

#### Step 0: Check if the correct files exist in the current directory. ####

[ ! -f 'FailedCommands' ] && [ ! -f 'recursive_trinity.cmds' ] &&
echo "
The necessary files do not exist in this directory. Perhaps you are in the wrong directory or haven't run Trinity yet.
:)
" && exit 1

#### Step 1: Create a backup of the recursive_trinity.cmds file in case we make a mistake, provided a backup doesn't already exist. ####

[ ! -d '../original_recursive_trinity.cmds' ] &&
mkdir ../original_recursive_trinity.cmds && dirq=0

[ ! -f '../original_recursive_trinity.cmds/recursive_trinity.cmds' ] &&
cp recursive_trinity.cmds ../original_recursive_trinity.cmds && recq=0

[ $recq -eq 0 ] &&
echo '
[Step 1]: Backup of recursive_trinity.cmds created in ../original_recursive_trinity.cmds.' ||
echo '
[Step 1]: Backup of recursive_trinity.cmds already exists in ../original_recursive_trinity.cmds.'

#### Step 2: Gather the command IDs from the FailedCommands file. ####

idar=( $(sed -n 's@\.trinity\.reads\.fa\.out.*@\\.@
                 s@.*/c@c@p' FailedCommands) )
echo "
[Step 2]: Gathered Command IDs from FailedCommands.
The following ${#idar[@]} commands were gathered:

${idar[@]%??}
"

#### Step 3: Delete failed commands from recursive_trinity.cmds following confirmation. ####

orlen=$(wc -l recursive_trinity.cmds)

read -p "Do you want to delete these commands? (y/N): " del
case $del in
    [yY])
    echo "
Deleting commands. Please wait..."
    ;;
    *)
    echo "The program has been cancelled. Nothing has been deleted."
    [ $recq -eq 0 ] && rm ../original_recursive_trinity.cmds/recursive_trinity.cmds &&
    [ $dirq -eq 0 ] && rmdir ../original_recursive_trinity.cmds
    exit 0
    ;;
esac

for n in ${idar[@]}
do
    sed -i "/$n/d" recursive_trinity.cmds
    ((++tic))
    per+="."
    echo $(( 100 * $tic / ${#idar[@]} ))"% Done
$per"
done

delen=$(wc -l recursive_trinity.cmds)

echo "
[Step 3]: Deleted ${#idar[@]} from recursive_trinity.cmds.
If all is well, you may now restart Trinity :)

File Length (Lines)
Before: $orlen
 After: $delen 
"
