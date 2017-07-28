#!/bin/bash
# This script is for jacob's convenience - personal to his machine
# This opens a new terminal and runs felix on GaAs_short on 2 parrallel cores
# saving the terminal output to terminal_log.txt
gnome-terminal --working-directory='/home/jacob/a/Felix/samples/GaAs_short' \
-e "bash -c \"mpirun -n 2 /home/jacob/a/Felix/src/felixrefine \
| tee terminal_log.txt\""


