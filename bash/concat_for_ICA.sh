#!/bin/bash
# Written by Anders Lillevik Thorsen, December 2021
# Concatenates GWC surfaces to one file for each hemisphere, which can then be inputted into ICA/ICASSO script

export SUBJECTS_DIR=/data/OBIC/Freesurfer
cd $SUBJECTS_DIR


	
for hemi in lh rh
	do
	mri_concat --f GWC_include20Dec2021_${hemi}.txt --o ICA/${hemi}.nu.w-g.avg_fsaverage.mgh
done
