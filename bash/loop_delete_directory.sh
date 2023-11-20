#!/bin/bash -x
# Written by Anders Lillevik Thorsen, October 2021
# Loops over list of folders and deletes them. Takes text file as input.

set -euo pipefail

# Sets up relevant paths and files
outputdir=/data/OBIC/Freesurfer # Where you want the Freesurfer data to end 

cd ${outputdir}
echo 'current directory is '`pwd` # Prints out current directory to be sure

for subject in `cat ${1}` # Loops through list of all subject
	do
	rm -R ${subject}
	done
