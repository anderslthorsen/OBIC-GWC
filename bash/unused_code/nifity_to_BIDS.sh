#!/bin/bash
set -euo pipefail
# By Anders Lillevik Thorsen, August 2021
# This script takes OBIC T1-weigthed nii.gz files and creates a minimal BIDS structure

headdir=/data/OBIC # Defines starting directory

cd ${headdir}/raw_data

# For each directory (OCD and HC)

for i in "T1_HC_n388" "T1_OCD_n446"
do

		for file in `cat ${i}_nii.txt`
		do
		#mkdir -p ${headdir}/BIDS/sub-${file}
			#mkdir -p ${headdir}/BIDS/sub-${file}/anat
				#echo ${headdir}/raw_data/${i}/${file}.nii.gz ${headdir}/BIDS/sub-${file}/anat/				
				#cp ${headdir}/raw_data/${i}/${file}.nii.gz ${headdir}/BIDS/sub-${file}/anat
				#mv ${headdir}/BIDS/sub-${file}/anat/${file}.nii.gz ${headdir}/BIDS/sub-${file}/anat/sub-${file}.nii.gz
				mv ${headdir}/BIDS/sub-${file}/anat/sub-${file}.nii.gz ${headdir}/BIDS/sub-${file}/anat/sub-${file}_T1w.nii.gz
		done
done
