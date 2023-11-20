#!/bin/bash -x
# Written by Anders Lillevik Thorsen, September 2021
# This scripts runs Freesurfer's recon-all with the cw256 flag for a given subject

set -euo pipefail

# Sets up relevant paths and files
outputdir=/data/OBIC/Freesurfer # Where you want the Freesurfer data to end up
subjectdir=/data/OBIC/BIDS # Where the original niftis are

cd ${outputdir}
echo 'current directory is '`pwd` # Prints out current directory to be sure

export SUBJECTS_DIR=${outputdir} # Sets directory where Freesurfer will output files
echo 'SUBJECTS_DIR is '$SUBJECTS_DIR

recon-all -s ${1} -i ${subjectdir}/${1}/anat/${1}_T1w.nii.gz -all -cw256 # Inputs subject IDs to recon-all
