#!/bin/bash -x
# Written by Anders Lillevik Thorsen, September 2021
# This script checks if the stats file from Freesurfer is present, and outputs a list where the ffile is missing. This can ba used to identify failed Freesurfer segmentatnons. The script also checks if the failed segmentation is due to a field of view above 256.

set -euo pipefail

# Sets up relevant paths
outputdir=/data/OBIC/Freesurfer # Where you want the Freesurfer data to end up

	cd ${outputdir}
	echo 'current directory is '`pwd` # Prints out current directory to be sure

# Creates relevant list files that are later used for looping over subjects or logging results
	ls -1 | grep "sub-*" > all.txt # Creates txt file with list of all subjects

	rm missing_asegs_stats.txt -f # Deletes old list to avoid appending to previous version
	echo "subject_ID number_of_slices" >> missing_asegs_stats.txt # appends headers to list


# Loops over subjects and finds those missing aseg.stats file
	for subject in `cat all.txt` # Loops through list of all subject
		do
		
			if [ ! -f ${subject}/stats/aseg.stats ]; then # Checks if aseg.stats file exists, if not then run code below
				nslices=$(mri_info ${subject}/mri/orig.mgz --nslices) # Calculates number of slices in original T1-weighted volume
				echo "${subject} ${nslices}" >> missing_asegs_stats.txt # If file does not exist, put subjectID and number of slices into list
			fi

		done


