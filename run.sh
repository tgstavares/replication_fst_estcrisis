#!/bin/bash

#make smm_estimation
#make main stats_show
#make main_exporters
#make greek_crisis_exp
make smm_estimation_par_exporters

echo
echo "  ####> Now the program..."
echo

#./smm_estimation
#./main && ./stats_show
#./stats_show
#./greek_crisis_exp
mpirun -n 81 --machinefile /home/ubuntu/.hosts_file.txt ./smm_estimation_par

echo "  ####> make clean before sync..."
echo
make clean
