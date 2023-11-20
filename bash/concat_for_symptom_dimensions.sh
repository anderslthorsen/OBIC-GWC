#!/bin/bash
# Written by Anders Lillevik Thorsen, December 2021
# Concatenates GWC surfaces to one file for each hemisphere, which can then be inputted into ICA/ICASSO script

export SUBJECTS_DIR=/data/OBIC/Freesurfer/
cd $SUBJECTS_DIR

for hemi in lh rh
	do
	mri_concat --f /data/OBIC/Freesurfer/GWC/symptom_dim_files_${hemi}.txt --o /data/OBIC/Freesurfer/PALM/symptom_dims/symptom_dims_${hemi}.nu.w-g.avg_fsaverage10.mgh
done
