#!/bin/bash
set -euo pipefail
# By Anders Lillevik Thorsen, August 2021
# This script takes Bergen T1-weigthed nii.gz files and copies them into the OBIC 


Bergendir=/data/ENIGMA/data
headdir=/data/OBIC/BIDS # Defines starting directory

cd $Bergendir

for file in `cat Bergen.txt`
do

	#mkdir -p ${headdir}/sub-09-${file}
echo ${headdir}/sub-09-${file}
mkdir -p ${headdir}/sub-09-${file}/anat
cp ${Bergendir}/${file}/T1_${file}_T0.nii.gz ${headdir}/sub-09-${file}/anat/sub-09-${file}_T1w.nii.gz
echo ${Bergendir}/${file}/T1_${file}_T0.nii.gz ${headdir}/sub-09-${file}/anat/sub-09-${file}_T1w.nii.gz

 #T1_00001_T0.nii.gz
	
done
