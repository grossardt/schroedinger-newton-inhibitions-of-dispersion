#!/bin/bash
#
# This code is part of Schroedinger Newton Inhibitions of Dispersion
# (C) Copyright Andre Grossardt 2010-2023
# https://github.com/grossardt/schroedinger-newton-inhibitions-of-dispersion
# 
# This code is licensed under the MIT License (see LICENSE.txt for details)
#
# This script compiles and runs the program.
# Can be called with 2 or 3 arguments:
# ./run.sh s                : starts a new calculation
# ./run.sh c path time      : continues from given path and timestep
# ./run.sh safe s           : compile in safe mode with fewer optimization
# ./run.sh safe c path time : same but continues calculation
# 
# Edit line 22 to set the path to the output files.
#
if [[ "$1" != "safe" ]]
then echo "Warning: you are in fast mode. Use 'run.sh safe' for safe mode."
fi
tms=`date +%s`
out="/PATH/TO/OUTPUT/FILES/sne${tms}"
if [[ "$1" == "safe" ]]
then gcc -lm -ffast-math -funroll-loops -march=core2 -o $out sne.c wf.c helpers.c
else gcc -O3 -DCHECK_OFF -lm -ffast-math -funroll-loops -march=core2 -o $out sne.c wf.c helpers.c
fi
date
if [[ "$1" == "safe" ]]
then $out $2 $3 $4
else $out $1 $2 $3
fi
rm $out
