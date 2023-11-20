#!/bin/bash
# Written by Anders Lillevik Thorsen, December 2021
# Loops over list of subjects to regenerate surfaces after gcut to remove dura from brain mask.

export SUBJECTS_DIR=/data/OBIC/Freesurfer
cd /data/OBIC/Freesurfer

for subject in `cat regen_surfaces.txt`
	do
		recon-all -autorecon-pial -subjid ${subject}
	done
