#!/bin/bash
# By Anders Lillevik Thorsen, August 2021
# This script converts all OBIC T1-weighted images to compressed nifti (nii.gz)

headdir=/data/OBIC/raw_data # Defines starting directory

# For each directory (OCD and HC)

for i in "T1_HC_n388" "T1_OCD_n446"
do
	echo $i
	cd $headdir
	pwd
	
		for file in `cat ${i}.txt`
		do		
		new_filename=`echo $file | sed 's/....$//'` # Removes .img suffix from filename
		echo $new_filename		
		
		mri_convert -ot nii ${i}/${file} ${i}/${new_filename}.nii.gz # Uses Freesurfer's tool to convert image to nii.gz

		done
done
