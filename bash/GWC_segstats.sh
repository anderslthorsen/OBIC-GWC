#!/bin/bash
# Written by Anders Lillevik Thorsen, March 2022
# This script extracts the gray-white matter contrast from cortical regions in Freesurfer's aparc regions (Desikan-Killiany atlas)

set -euo pipefail

# Sets up relevant paths and files
data_dir=/data/OBIC/Freesurfer # Where you want the Freesurfer data to end up
output_dir=/data/OBIC/segstats # Where the output ends up

cd ${data_dir}
echo 'current directory is '`pwd` # Prints out current directory to be sure

export SUBJECTS_DIR=${data_dir} # Set subject directory for Freesurfers
echo 'SUBJECTS_DIR is '$SUBJECTS_DIR

# Define test subject
#subject=sub-01subject001

# Loop over subjects using a line seperated list of subject IDs
for subject in `cat GWC/GWC_include20Dec2021.txt`; do

	# Loop over hemispheres
	for hemi in lh rh; do

		# Run mri_segstats for cortical thickness (using subject aparc)
		mri_segstats --annot ${subject} ${hemi} aparc --i ${SUBJECTS_DIR}/${subject}/surf/${hemi}.thickness --sum ${SUBJECTS_DIR}/${subject}/stats/${hemi}.thickness.DK
	
		# Run mri_segstats for subcortical volume (using subject aseg)
		mri_segstats --annot ${subject} ${hemi} aseg --i ${SUBJECTS_DIR}/${subject}/surf/${hemi}.thickness --sum ${SUBJECTS_DIR}/${subject}/stats/${hemi}.thickness.DK

		# Run mri_segstats for GWC (using subject aparc)
		mri_segstats --annot $subject ${hemi} aparc --i ${SUBJECTS_DIR}/${subject}/surf/${hemi}.nu.w-g.avg.mgh --sum ${SUBJECTS_DIR}/${subject}/stats/${hemi}.GWC.DK
		
		# Run mri_segstats for WM intensity (using subject aparc)
		mri_segstats --annot $subject ${hemi} aparc --i ${SUBJECTS_DIR}/${subject}/surf/${hemi}.nu.avg.wm.mgh --sum ${SUBJECTS_DIR}/${subject}/stats/${hemi}.WM.DK
		
		# Run mri_segstats for GM intensity (using subject aparc)
		mri_segstats --annot $subject ${hemi} aparc --i ${SUBJECTS_DIR}/${subject}/surf/${hemi}.nu.avg.gm.mgh --sum ${SUBJECTS_DIR}/${subject}/stats/${hemi}.GM.DK
	done
done

cd ${output_dir}

# Run aseg2stats to get a table of the results

# Cortical thickness
asegstats2table --subjectsfile=/data/OBIC/Freesurfer/GWC/GWC_include20Dec2021.txt --stats=lh.thickness.DK -m mean --tablefile N848_thickness.lh.DK.csv -v
asegstats2table --subjectsfile=/data/OBIC/Freesurfer/GWC/GWC_include20Dec2021.txt --stats=rh.thickness.DK -m mean --tablefile N848_thickness.rh.DK.csv -v

# Subcortical volume
asegstats2table --subjectsfile=/data/OBIC/Freesurfer/GWC/GWC_include20Dec2021.txt --tablefile N848_subcortical_volume.bilateral.DK.csv -v

# GWC
asegstats2table --subjectsfile=/data/OBIC/Freesurfer/GWC/GWC_include20Dec2021.txt --stats=lh.GWC.DK -m mean --tablefile N848_GWC.lh.DK.csv -v
asegstats2table --subjectsfile=/data/OBIC/Freesurfer/GWC/GWC_include20Dec2021.txt --stats=rh.GWC.DK -m mean --tablefile N848_GWC.rh.DK.csv -v

# GM signal intensity
asegstats2table --subjectsfile=/data/OBIC/Freesurfer/GWC/GWC_include20Dec2021.txt --stats=lh.GM.DK -m mean --tablefile N848_GM.lh.DK.csv -v
asegstats2table --subjectsfile=/data/OBIC/Freesurfer/GWC/GWC_include20Dec2021.txt --stats=rh.GM.DK -m mean --tablefile N848_GM.rh.DK.csv -v

# WM signal intensity
asegstats2table --subjectsfile=/data/OBIC/Freesurfer/GWC/GWC_include20Dec2021.txt --stats=lh.WM.DK -m mean --tablefile N848_WM.lh.DK.csv -v
asegstats2table --subjectsfile=/data/OBIC/Freesurfer/GWC/GWC_include20Dec2021.txt --stats=rh.WM.DK -m mean --tablefile N848_WM.rh.DK.csv -v








