#!/bin/bash
set -euo pipefail
# By Anders Lillevik Thorsen, August 2021
# This script changes the sub formatting to be BIDS compatible. 

headdir=/data/OBIC/BIDS # Defines starting directory

cd ${headdir}

# For each directory


for sub in `cat list.txt`
#for sub in "sub-01-001"
do
#echo ${headdir}/${sub}/anat
cd ${headdir}/${sub}/anat

mv "${sub}-T1w.nii.gz" "$(echo "${sub}-T1w.nii.gz" | sed 's/.\{11\}$//')"
mv ${sub} ${sub}_T1w.nii.gz

#for i in *_*;do mv $i ${i//"_"/"-"};done


# filename = sub-01-001-T1w.nii.gz
# sub = sub-01-001

done

