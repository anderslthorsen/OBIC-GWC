#!/bin/bash
headdir=/data/OBIC/BIDS

cd ${headdir}

for i in `cat list.txt`
	do
	cp -n ${headdir}/${i}/anat/${i}_T1w.nii.gz ${headdir}/${i}/anat/${i}_T1w_preorientation.nii.gz
	fslreorient2std ${headdir}/${i}/anat/${i}_T1w.nii.gz
	done
