#!/bin/bash -x
# Written by Anders Lillevik Thorsen, September 2021
# This scripts runs Freesurfer in paralell on three test batches of 4, 6 and 8 subjects in parallel. The aim is to determine how much RAM is used, and what the reasonable number of maximum parallel subjects are.

set -euo pipefail

# Available resources: 20 CPU 30Gb RAM
# 8 GB RAM per subject = 3.75 subjects in parallel
# 6 GB RAM per subject = 5 subjects in parallel
# 4 GB RAM per subject = 7.5 subjects in parallel
# 3 GB RAM per subject = 10 subjects in parallel
# 2.5 GB RAM per subject = 12 subjects in parallel
# 4. September 2021. Current max observed RAM usage per subject is  about 2.1 GB. Might be possible to run up to 12 subject at once, will try 10 to be safe.
# Batch of 4 took 5 hours 20 min
# Batch of 6 about 5 hours 20 min
# Batch of 8 took about 5 hours 20 min

# Sets up relevant paths and files
outputdir=/data/OBIC/Freesurfer # Where you want the Freesurfer data to end up
subjectdir=/data/OBIC/BIDS # Where the original niftis are

cd ${outputdir}
echo 'current directory is '`pwd` # Prints out current directory to be sure

export SUBJECTS_DIR=${outputdir} # Sets directory where Freesurfer will output files
echo 'SUBJECTS_DIR is '$SUBJECTS_DIR

	for i in {4..8..2} # Starts at listn 4, moves up 2 each iteration, ends at listn 8
	do
		echo 'Number of paralell jobs is '${i}
		echo 'Current subjects are '`cat listn${i}.txt` # Looks for subject IDs in files called listn.txt
		
		cat listn${i}.txt | parallel --jobs ${i} recon-all -s {} -i ${subjectdir}/{}/anat/{}_T1w.nii.gz -all # Inputs subject IDs to parallel
		#break # Used for debugging, makes sure that only the first iteration of the loop runs
	done


