#!/bin/bash
# Written by Anders Lillevik Thorsen, May 2022
# Concatenates GWC surfaces for gray and white matter to one file for each hemisphere

export SUBJECTS_DIR=/data/OBIC/Freesurfer/
cd $SUBJECTS_DIR

dcm=nu
target=fsaverage
fwhm=10

echo "SUBJECTS_DIR is $SUBJECTS_DIR"

#for subject in `cat /data/OBIC/Freesurfer/GWC/GWC_include20Dec2021.txt`
#do#
#		echo "Subject is ${s}"
#		cd $SUBJECTS_DIR
#		cd ${subject}
#
#		#smooth GM and WM intensity
#		for meas in avg.gm avg.wm
#		do
#			for hemi in lh rh
#			do
#
#				mri_surf2surf \
#					--hemi $hemi \
#					--s $target \
#					--sval $SUBJECTS_DIR/${subject}/surf/${hemi}.${dcm}.${meas}_fsaverage.mgh \
#					--fwhm $fwhm \
#					--cortex \
#					--tval $SUBJECTS_DIR/${subject}/surf/${hemi}.${dcm}.${meas}_fsaverage.${fwhm}.mgh
#			done
#		done
#done

cd $SUBJECTS_DIR

echo "current directory is $(pwd)"	

#for hemi in lh rh
#do
#	mri_concat --f /data/OBIC/Freesurfer/GWC/GWC_include20Dec2021_${hemi}_GM.txt --o ICA/${hemi}.nu.gm.avg_fsaverage.mgh
#	mri_concat --f /data/OBIC/Freesurfer/GWC/GWC_include20Dec2021_${hemi}_WM.txt --o ICA/${hemi}.nu.wm.avg_fsaverage.mgh
#done

# Calculate mean intensity for GM and WM per hemisphere
mri_concat \
--i /data/OBIC/Freesurfer/ICA/lh.nu.gm.avg_fsaverage.mgh \
--o /data/OBIC/Freesurfer/estimates/GM_lh_mean.mgz \
--mean

mri_concat \
--i /data/OBIC/Freesurfer/ICA/rh.nu.gm.avg_fsaverage.mgh \
--o /data/OBIC/Freesurfer/estimates/GM_rh_mean.mgz \
--mean

mri_concat \
--i /data/OBIC/Freesurfer/ICA/lh.nu.wm.avg_fsaverage.mgh \
--o /data/OBIC/Freesurfer/estimates/WM_lh_mean.mgz \
--mean

mri_concat \
--i /data/OBIC/Freesurfer/ICA/rh.nu.wm.avg_fsaverage.mgh \
--o /data/OBIC/Freesurfer/estimates/WM_rh_mean.mgz \
--mean
