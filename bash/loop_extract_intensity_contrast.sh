#!/bin/bash
# Written by Anders Lillevik Thorsen, December 2021
# Loops over list of subjects and runs the script for calculating GWC for each subject/line.

cd /data/OBIC/Freesurfer

for subject in `cat GWC_include20Dec2021.txt`
	do
		bash /data/OBIC/scripts/extract_intensity_contrast_bash.sh ${subject}
	done
